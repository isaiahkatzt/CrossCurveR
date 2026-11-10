#############################################
#         Shock States and Features         #
#############################################

st_factor_names <- function(factor_dim) {
  if (factor_dim == 2) {
    return(c("L", "S"))
  }
  if (factor_dim == 3) {
    return(c("L", "S", "C"))
  }
  stop("Factor matrices must have either 2 or 3 columns.")
}

.st_factor_matrix <- function(x, model_factors = st_factor_names(ncol(as.matrix(x)))) {
  x <- as.matrix(x)
  if (ncol(x) != length(model_factors)) {
    stop("All factor matrices must have the same number of columns.", call. = FALSE)
  }
  labels <- colnames(x)
  if (is.null(labels) || identical(labels, model_factors)) return(x)
  identities <- c(L = "L", S = "S", C = "C", l1 = "L", l2 = "S", l3 = "C")[labels]
  if (anyNA(identities) || anyDuplicated(identities) || !setequal(identities, model_factors)) {
    stop("Factor columns must uniquely identify ", paste(model_factors, collapse = ", "), ".", call. = FALSE)
  }
  x <- x[, match(model_factors, identities), drop = FALSE]
  colnames(x) <- model_factors
  x
}

st_ns_factors <- function(bl_yields, curves) {
  stats::setNames(lapply(curves, function(curve_name) {
    factor_mat <- .st_factor_matrix(bl_yields$betas[[curve_name]]$NSbetas)
    colnames(factor_mat) <- st_factor_names(ncol(factor_mat))
    factor_mat
  }), curves)
}

st_check_yields <- function(yields, curves, reference) {
  yields <- .cc_named(yields, curves, "`yields`")
  .cc_reference(reference, curves)
  yields
}

st_sc_state <- function(yields, lambdas, cutoffs, curves, mats) {
  yields <- .cc_yields(yields, mats, curves)
  ytime <- yields[[1]][[1]]
  single_curve <- .bl_prepare_single_curve_state(yields, lambdas, cutoffs, curves, mats)

  return(list(
    curves = curves,
    mats = mats,
    time = ytime,
    mds = single_curve$mds,
    lambdas = single_curve$lambdas,
    yields = yields,
    bl_yields = single_curve$bl_yields,
    phi_hat = single_curve$phi_hat,
    W = single_curve$W,
    H = single_curve$H,
    ns_factor = st_ns_factors(single_curve$bl_yields, curves),
    input_spec = .cc_input_spec(yields, mats, curves,
      decay = list(selected = single_curve$lambdas, grid = lambdas,
                   parameterization = "inverse months"))
  ))
}

st_xt_ecm <- function(bl_yields, curves, reference, curve_ecm, X_normalize = TRUE,
                      X_trunc = FALSE) {
  cc_betas <- blmc_VAR(bl_yields, reference, curves)

  cc_cspread <- stats::setNames(lapply(names(cc_betas), function(curve_name) {
    ecm_obj <- curve_ecm[[curve_name]]
    if (is.list(ecm_obj) && is.numeric(ecm_obj$rank) &&
        length(ecm_obj$rank) == 1 && ecm_obj$rank > 0) {
      blmc_cspread(
        cc_betas = cc_betas[[curve_name]],
        Gmatrix = ecm_obj$Gmatrix,
        curve_name = curve_name,
        normalize = X_normalize
      )
    }
  }), names(cc_betas))

  Xt <- blcc_build_Xt(cc_cspread, curve_ecm, trunc = X_trunc)

  return(list(
    Xt = Xt,
    cc_ci_features = list(
      cc_cspread = cc_cspread,
      cc_ecm = curve_ecm,
      cc_beta = cc_betas
    )
  ))
}

st_xt_state <- function(bl_yields, curves, reference, reuse_ECM = TRUE, baseline_ecm = NULL,
                        ECM_estim = "ML", ECM_type = "eigen", ECM_alpha = 0.1,
                        ECM_fallback = c("error", "baseline"),
                        X_normalize = TRUE, X_trunc = FALSE,
                        Xt_smooth = FALSE, smoother = c("ns", "bs", "rm", "hp", "henderson"),
                        knot_count = 5, k_count = 9, k_pass = 3,
                        hp_lambda = 1600, henderson_k = 13, ytime) {
  ECM_fallback <- match.arg(ECM_fallback)
  baseline_xt <- function() st_xt_ecm(
    bl_yields = bl_yields, curves = curves, reference = reference,
    curve_ecm = baseline_ecm, X_normalize = X_normalize, X_trunc = X_trunc
  )

  if (reuse_ECM) {
    if (is.null(baseline_ecm)) {
      stop("`baseline_ecm` must be supplied when `reuse_ECM = TRUE`.")
    }
    xt_state <- baseline_xt()
  } else {
    cc_ci_features <- NULL
    xt_state <- tryCatch(
      {
        cc_ci_features <- blcc_cointegration(
          blsc_list = bl_yields, curves = curves, reference = reference,
          estim = ECM_estim, type = ECM_type, alpha = ECM_alpha, normalize = X_normalize)
        list(Xt = blcc_build_Xt(cc_ci_features$cc_cspread, cc_ci_features$cc_ecm, trunc = X_trunc),
             cc_ci_features = cc_ci_features)
      },
      error = function(e) {
        if (ECM_fallback != "baseline" || is.null(baseline_ecm)) {
          stop(e)
        }

        warning(
          "Re-estimated ECM did not produce usable Xt features; falling back to baseline ECM.",
          call. = FALSE
        )

        fallback_state <- baseline_xt()
        fallback_state$cc_ci_features$refit_cc_ci_features <- cc_ci_features
        fallback_state$cc_ci_features$ECM_fallback <- "baseline"
        fallback_state$cc_ci_features$ECM_fallback_reason <- conditionMessage(e)
        fallback_state
      }
    )
  }

  if (Xt_smooth) {
    xt_state$Xt <- blce_smooth_Xt(
      Xt = xt_state$Xt,
      ytime = ytime,
      smoother = smoother,
      knot_count = knot_count,
      k_count = k_count,
      k_pass = k_pass,
      hp_lambda = hp_lambda,
      henderson_k = henderson_k
    )
  }

  xt_state$Xt <- .cc_features(xt_state$Xt, ytime)
  xt_state$feature_spec <- .bl_feature_spec(
    xt_state$cc_ci_features, X_normalize, X_trunc, Xt_smooth, smoother,
    knot_count, k_count, k_pass, hp_lambda, henderson_k)

  xt_state
}

st_shock_state <- function(single_curve_state, shocked_factors) {
  curves <- single_curve_state$curves
  factor_names <- st_factor_names(ncol(as.matrix(shocked_factors[[1]])))

  shocked_curve_fits <- stats::setNames(lapply(curves, function(curve_name) {
    curve_phi <- single_curve_state$bl_yields$phi[[curve_name]]
    baseline_nsfit <- single_curve_state$bl_yields$nsfit[[curve_name]]
    shocked_beta <- .st_factor_matrix(shocked_factors[[curve_name]], factor_names)
    colnames(shocked_beta) <- factor_names

    shocked_nsfit <- shocked_beta %*% curve_phi$phi
    # Shocked curve panels are measured against the independent DNS fit path.
    shocked_matrix_yield <- shocked_nsfit

    list(
      phi = curve_phi,
      betas = list(
        NSbetas = shocked_beta,
        m_yield = shocked_matrix_yield
      ),
      nsfit = list(
        NSfit = shocked_nsfit,
        W = baseline_nsfit$W
      ),
      H = single_curve_state$bl_yields$H[[curve_name]]
    )
  }), curves)

  shocked_bl_yields <- purrr::transpose(shocked_curve_fits)

  shocked_yields <- stats::setNames(lapply(curves, function(curve_name) {
    curve_time <- single_curve_state$yields[[curve_name]][[1]]
    curve_cols <- colnames(single_curve_state$yields[[curve_name]])[-1]
    curve_df <- data.frame(
      time = curve_time,
      as.data.frame(shocked_bl_yields$betas[[curve_name]]$m_yield)
    )
    colnames(curve_df) <- c("time", curve_cols)
    curve_df
  }), curves)

  return(list(
    yields = shocked_yields,
    bl_yields = shocked_bl_yields,
    phi_hat = single_curve_state$phi_hat,
    W = blcc_build_Wj(shocked_bl_yields$nsfit, single_curve_state$mats),
    ns_factor = stats::setNames(lapply(curves, function(curve_name) {
      shocked_bl_yields$betas[[curve_name]]$NSbetas
    }), curves)
  ))
}

st_factor_shock <- function(curves, curve_dns_factor,
                            shock_curves, shock_factors, shock_magnitude,
                            shock_type = c("additive", "multiplicative"),
                            shock_window = NULL,
                            shock_window_position = c("first", "last")) {
  shock_type <- match.arg(shock_type)
  shock_window_position <- match.arg(shock_window_position)

  apply_shock <- switch(shock_type, additive = st_additive, multiplicative = st_multiplicative)
  apply_shock(
    curves = curves,
    curve_dns_factor = curve_dns_factor,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_window = shock_window,
    shock_window_position = shock_window_position
  )
}

st_factor_profile <- function(curves, curve_dns_factor,
                              shock_curves, shock_factors,
                              shock_magnitudes,
                              shock_type = c("additive", "multiplicative")) {
  shock_type <- match.arg(shock_type)

  if (!is.list(curve_dns_factor) || length(curve_dns_factor) == 0 ||
      is.null(names(curve_dns_factor))) {
    stop("`curve_dns_factor` must be a non-empty named list.")
  }
  if (length(shock_factors) != 1) {
    stop("Profile shocks currently support one shocked factor at a time.")
  }
  if (!all(shock_curves %in% curves)) {
    stop("`shock_curves` must be a subset of `curves`.")
  }

  factor_dim <- ncol(as.matrix(curve_dns_factor[[1]]))
  model_factors <- st_factor_names(factor_dim)
  if (!shock_factors %in% model_factors) {
    stop("`shock_factors` must be compatible with the model factor dimension.")
  }

  factor_idx <- match(shock_factors, model_factors)
  shock_magnitudes <- as.numeric(shock_magnitudes)
  if (any(!is.finite(shock_magnitudes))) {
    stop("`shock_magnitudes` must contain only finite values.", call. = FALSE)
  }

  shocked_curve_dns_factor <- stats::setNames(lapply(names(curve_dns_factor), function(curve_name) {
    curve_factor_mat <- .st_factor_matrix(curve_dns_factor[[curve_name]], model_factors)
    if (nrow(curve_factor_mat) != length(shock_magnitudes)) {
      stop("Shock profile length must match the number of factor rows.")
    }
    shocked_mat <- curve_factor_mat
    if (curve_name %in% shock_curves) {
      if (shock_type == "additive") {
        shocked_mat[, factor_idx] <- curve_factor_mat[, factor_idx] + shock_magnitudes
      } else {
        shocked_mat[, factor_idx] <- curve_factor_mat[, factor_idx] * (1 + shock_magnitudes)
      }
    }

    colnames(shocked_mat) <- model_factors
    shocked_mat
  }), names(curve_dns_factor))

  list(
    shock_profile = shock_magnitudes,
    shocked_factors = shocked_curve_dns_factor
  )
}

st_dns_shock <- function(single_curve_state,
                         shock_curves, shock_factors, shock_magnitude,
                         shock_type = c("additive", "multiplicative"),
                         shock_window = NULL,
                         shock_window_position = c("first", "last")) {
  shock_fit <- st_factor_shock(
    curves = single_curve_state$curves,
    curve_dns_factor = single_curve_state$ns_factor,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window,
    shock_window_position = shock_window_position
  )

  shocked_state <- st_shock_state(
    single_curve_state = single_curve_state,
    shocked_factors = shock_fit$shocked_factors
  )

  return(list(
    yields = shocked_state$yields,
    ns_factor = shocked_state$ns_factor,
    shock_matrix = shock_fit$shock_matrix,
    shocked_factors = shock_fit$shocked_factors,
    baseline_state = single_curve_state
  ))
}

st_mean_profile <- function(curve_factor, shock_fraction_profile,
                            shock_factor = "L", curves = names(curve_factor)) {
  if (!is.list(curve_factor) || length(curve_factor) == 0 || is.null(names(curve_factor))) {
    stop("`curve_factor` must be a non-empty named list.")
  }

  if (!is.character(shock_factor) || length(shock_factor) != 1 || is.na(shock_factor) ||
      !shock_factor %in% c("L", "S", "C")) {
    stop("`shock_factor` must be one of L, S or C.")
  }

  curve_factor <- curve_factor[curves]
  factor_panel <- do.call(cbind, lapply(curve_factor, function(factor_mat) {
    factor_mat <- as.matrix(factor_mat)
    factor_mat <- .st_factor_matrix(factor_mat,
      if (ncol(factor_mat) == 1L) shock_factor else st_factor_names(ncol(factor_mat)))
    if (!shock_factor %in% colnames(factor_mat)) {
      stop("All factor matrices must contain `shock_factor`.")
    }
    factor_mat[, shock_factor]
  }))

  shock_fraction_profile <- as.numeric(shock_fraction_profile)
  if (any(!is.finite(shock_fraction_profile))) {
    stop("`shock_fraction_profile` must contain only finite values.", call. = FALSE)
  }
  if (nrow(factor_panel) != length(shock_fraction_profile)) {
    stop("Shock fraction profile length must match the number of factor rows.")
  }

  as.numeric(rowMeans(factor_panel) * shock_fraction_profile)
}

st_align_Xt <- function(Xt_new, Xt_reference) {
  Xt_new <- .cc_features(Xt_new)
  Xt_reference <- .cc_features(Xt_reference)
  if (nrow(Xt_new) != nrow(Xt_reference)) {
    stop("Scenario and baseline Xt must have the same observation count.", call. = FALSE)
  }
  if (ncol(Xt_new) != ncol(Xt_reference)) {
    stop(
      "Re-estimated Xt is incompatible with the reused baseline BS0. ",
      "Set `reuse_BS0 = FALSE` or `reuse_ECM = TRUE`."
    )
  }

  ref_names <- colnames(Xt_reference)
  new_names <- colnames(Xt_new)

  if (!is.null(ref_names) && !is.null(new_names)) {
    if (!setequal(ref_names, new_names)) {
      stop(
        "Re-estimated Xt columns do not match the reused baseline BS0 design. ",
        "Set `reuse_BS0 = FALSE` or `reuse_ECM = TRUE`."
      )
    }
    Xt_new <- Xt_new[, ref_names, drop = FALSE]
  }

  Xt_new
}

st_project_fit <- function(yields, time, phi_hat, Xt, BS0, curves, mats, W) {
  yields <- .cc_yields(yields, mats, curves)
  if (!identical(.cc_index(time), .cc_index(yields[[1]][[1]]))) {
    stop("Fit time index must match the yield observation index.", call. = FALSE)
  }
  Xt <- .cc_features(Xt, time)
  Y <- vecY(yields, curves, mats)
  P <- bln_build_P(curves, mats)
  pY <- PY_full(Y$tY, P)

  sigma_JT <- tSigma_optim(BS0 = BS0, X = Xt)
  full_sigma_t <- sqrt_inv_build(sigma_JT$tSigma)

  sqrt_sigma_t <- full_sigma_t$sqrt
  sqrt_inv_sigma_t <- full_sigma_t$inverse

  projected <- .bl_project_rescale(pY$py, P, phi_hat, sqrt_sigma_t, sqrt_inv_sigma_t,
                                   time = time, curves = curves, mats = mats)

  return(list(
    tenor = projected$tenor,
    curve = projected$curve,
    sigma_JT = sigma_JT,
    Xt = Xt,
    W = W,
    input_spec = .cc_input_spec(yields, mats, curves, Xt)
  ))
}

st_resolve_config <- function(value, fallback) {
  if (is.null(value)) {
    fallback
  } else {
    value
  }
}

st_check_baseline <- function(baseline_state, curves = NULL, mats = NULL, reference = NULL) {
  required_names <- c(
    "curves", "mats", "reference", "single_curve_state",
    "xt_state", "BS0", "baseline_fit", "config"
  )

  if (!is.list(baseline_state) || !all(required_names %in% names(baseline_state))) {
    stop(
      "`baseline_state` must be created by `mc_fit_shock_baseline()` and contain ",
      paste(required_names, collapse = ", "),
      "."
    )
  }

  if (!is.null(curves) && !identical(as.character(curves), as.character(baseline_state$curves))) {
    stop("Supplied `curves` do not match the curves stored in `baseline_state`.")
  }

  if (!is.null(mats) && !identical(as.numeric(mats), as.numeric(baseline_state$mats))) {
    stop("Supplied `mats` do not match the maturities stored in `baseline_state`.")
  }

  if (!is.null(reference) && !identical(reference, baseline_state$reference)) {
    stop("Supplied `reference` does not match the reference curve stored in `baseline_state`.")
  }

  observed <- .cc_yields(baseline_state$single_curve_state$yields,
                         baseline_state$mats, baseline_state$curves, baseline_state$reference)
  time <- observed[[1]][[1]]
  if (!identical(.cc_index(time), .cc_index(baseline_state$single_curve_state$time))) {
    stop("Cached single-curve time index does not match its yields.", call. = FALSE)
  }
  X <- .cc_features(baseline_state$xt_state$Xt, time)
  .st_check_cached_coefficients(baseline_state$BS0, X, baseline_state$curves, baseline_state$mats)
  spec <- baseline_state$input_spec
  if (!is.null(spec) &&
      (!identical(spec$curves, baseline_state$curves) ||
       !identical(as.numeric(spec$mats), as.numeric(baseline_state$mats)) ||
       !identical(.cc_index(spec$time), .cc_index(time)) ||
       !identical(spec$features$columns, colnames(X)) || spec$features$count != ncol(X))) {
    stop("Cached baseline inputs no longer match its stored input specification.", call. = FALSE)
  }
  stored_features <- spec$features$preprocessing
  if (!is.null(stored_features) && !isTRUE(all.equal(
      stored_features, baseline_state$xt_state$feature_spec, tolerance = 0))) {
    stop("Cached feature transformations no longer match the stored input specification.", call. = FALSE)
  }

  baseline_state
}

#############################################
#           Compatibility Aliases           #
#############################################

st_extract_ns_factor <- st_ns_factors
st_validate_yields_input <- st_check_yields
st_prepare_single_curve_state <- st_sc_state
st_build_xt_from_ecm <- st_xt_ecm
st_build_xt_state <- st_xt_state
st_build_shocked_curve_state <- st_shock_state
st_apply_factor_shock <- st_factor_shock
st_apply_factor_shock_profile <- st_factor_profile
st_dns_shock_from_single_curve_state <- st_dns_shock
st_average_level_shock_profile <- st_mean_profile
st_align_Xt_to_baseline <- st_align_Xt
st_finalize_shocked_fit <- st_project_fit
st_validate_baseline_state <- st_check_baseline
