#############################################
#          Shocked Cross-Curve Fit          #
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

st_extract_ns_factor <- function(bl_yields, curves) {
  setNames(lapply(curves, function(curve_name) {
    factor_mat <- as.matrix(bl_yields$betas[[curve_name]]$NSbetas)
    colnames(factor_mat) <- st_factor_names(ncol(factor_mat))
    factor_mat
  }), curves)
}

st_validate_yields_input <- function(yields, curves, reference) {
  if (!is.character(curves) || length(curves) == 0 || anyNA(curves) || any(!nzchar(curves))) {
    stop("`curves` must be a non-empty character vector.")
  }

  if (!is.list(yields) || length(yields) != length(curves)) {
    stop("`yields` must be a list with the same length as `curves`.")
  }

  if (is.null(names(yields)) || anyNA(names(yields)) || any(!nzchar(names(yields)))) {
    stop("`yields` must be a named list.")
  }

  if (!all(curves %in% names(yields))) {
    stop("All `curves` must be present in the names of `yields`.")
  }

  if (!reference %in% curves) {
    stop("`reference` must be one of the supplied `curves`.")
  }

  yields[curves]
}

st_prepare_single_curve_state <- function(yields, lambdas, cutoffs, curves, mats) {
  mat_str <- sapply(mats, numeric_to_matname)
  mds <- as.vector(t(outer(curves, mat_str, paste, sep = ".")))
  ytime <- yields[[1]][[1]]
  curve_lambdas <- vapply(yields, function(curve_yield) {
    lambda_grid_search(lambdas, curve_yield, mats)
  }, numeric(1))

  curve_fits <- setNames(lapply(seq_along(curves), function(i) {
    bl_single_curve(curve_lambdas[[i]], yields[[i]], cutoffs, mats)
  }), curves)

  bl_yields <- purrr::transpose(curve_fits)
  phi_hat <- bln_phi_hat(bl_yields$phi, mds)
  W <- blcc_build_Wj(bl_yields$nsfit, mats)
  H <- bdiag(bl_yields$H)

  return(list(
    curves = curves,
    mats = mats,
    time = ytime,
    mds = mds,
    lambdas = curve_lambdas,
    yields = yields,
    bl_yields = bl_yields,
    phi_hat = phi_hat,
    W = W,
    H = H,
    ns_factor = st_extract_ns_factor(bl_yields, curves)
  ))
}

st_build_xt_from_ecm <- function(bl_yields, curves, reference, curve_ecm, X_normalize = TRUE,
                                 X_trunc = FALSE) {
  cc_betas <- blmc_VAR(bl_yields, reference, curves)

  cc_cspread <- setNames(lapply(names(cc_betas), function(curve_name) {
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

st_build_xt_state <- function(bl_yields, curves, reference, reuse_ECM = TRUE, baseline_ecm = NULL,
                              ECM_estim = "ML", ECM_type = "eigen", ECM_alpha = 0.1,
                              X_normalize = TRUE, X_trunc = FALSE,
                              Xt_smooth = FALSE, smoother = c("ns", "bs", "rm", "hp", "henderson"),
                              knot_count = 5, k_count = 9, k_pass = 3,
                              hp_lambda = 1600, henderson_k = 13, ytime) {
  if (reuse_ECM) {
    if (is.null(baseline_ecm)) {
      stop("`baseline_ecm` must be supplied when `reuse_ECM = TRUE`.")
    }
    xt_state <- st_build_xt_from_ecm(
      bl_yields = bl_yields,
      curves = curves,
      reference = reference,
      curve_ecm = baseline_ecm,
      X_normalize = X_normalize,
      X_trunc = X_trunc
    )
  } else {
    cc_ci_features <- blcc_cointegration(
      blsc_list = bl_yields,
      curves = curves,
      reference = reference,
      estim = ECM_estim,
      type = ECM_type,
      alpha = ECM_alpha,
      normalize = X_normalize
    )
    xt_state <- list(
      Xt = blcc_build_Xt(cc_ci_features$cc_cspread, cc_ci_features$cc_ecm, trunc = X_trunc),
      cc_ci_features = cc_ci_features
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

  xt_state
}

st_build_shocked_curve_state <- function(single_curve_state, shocked_factors) {
  curves <- single_curve_state$curves
  factor_names <- st_factor_names(ncol(as.matrix(shocked_factors[[1]])))

  shocked_curve_fits <- setNames(lapply(curves, function(curve_name) {
    curve_phi <- single_curve_state$bl_yields$phi[[curve_name]]
    baseline_nsfit <- single_curve_state$bl_yields$nsfit[[curve_name]]
    shocked_beta <- as.matrix(shocked_factors[[curve_name]])
    colnames(shocked_beta) <- factor_names

    shocked_nsfit <- shocked_beta %*% curve_phi$phi
    shocked_matrix_yield <- shocked_nsfit - baseline_nsfit$W

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

  shocked_yields <- setNames(lapply(curves, function(curve_name) {
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
    ns_factor = setNames(lapply(curves, function(curve_name) {
      shocked_bl_yields$betas[[curve_name]]$NSbetas
    }), curves)
  ))
}

st_apply_factor_shock <- function(curves, curve_dns_factor,
                                  shock_curves, shock_factors, shock_magnitude,
                                  shock_type = c("additive", "multiplicative"),
                                  shock_window = NULL) {
  shock_type <- match.arg(shock_type)

  switch(
    shock_type,
    additive = st_additive(
      curves = curves,
      curve_dns_factor = curve_dns_factor,
      shock_curves = shock_curves,
      shock_factors = shock_factors,
      shock_magnitude = shock_magnitude,
      shock_window = shock_window
    ),
    multiplicative = st_multiplicative(
      curves = curves,
      curve_dns_factor = curve_dns_factor,
      shock_curves = shock_curves,
      shock_factors = shock_factors,
      shock_magnitude = shock_magnitude,
      shock_window = shock_window
    )
  )
}

st_dns_shock_from_single_curve_state <- function(single_curve_state,
                                                 shock_curves, shock_factors, shock_magnitude,
                                                 shock_type = c("additive", "multiplicative"),
                                                 shock_window = NULL) {
  shock_fit <- st_apply_factor_shock(
    curves = single_curve_state$curves,
    curve_dns_factor = single_curve_state$ns_factor,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window
  )

  shocked_state <- st_build_shocked_curve_state(
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

st_build_shocked_dly_single_curve <- function(dly_single_curve, shocked_factors,
                                              curves = names(dly_single_curve)) {
  factor_names <- st_factor_names(ncol(as.matrix(shocked_factors[[1]])))

  setNames(lapply(curves, function(curve_name) {
    baseline_curve <- dly_single_curve[[curve_name]]
    shocked_beta <- as.matrix(shocked_factors[[curve_name]])
    colnames(shocked_beta) <- factor_names
    shocked_nsfit <- shocked_beta %*% baseline_curve$phi$phi

    baseline_curve$betas$NSbetas <- shocked_beta
    baseline_curve$betas$m_yield <- shocked_nsfit
    baseline_curve$nsfit$NSfit <- shocked_nsfit
    baseline_curve
  }), curves)
}

st_dly_shock_from_fit <- function(dly_fit_obj,
                                  shock_curves, shock_factors, shock_magnitude,
                                  shock_type = c("additive", "multiplicative"),
                                  shock_window = NULL) {
  shock_fit <- st_apply_factor_shock(
    curves = dly_fit_obj$curves,
    curve_dns_factor = dly_pull_beta(dly_fit_obj$single_curve, curves = dly_fit_obj$curves),
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window
  )

  shocked_single_curve <- st_build_shocked_dly_single_curve(
    dly_single_curve = dly_fit_obj$single_curve,
    shocked_factors = shock_fit$shocked_factors,
    curves = dly_fit_obj$curves
  )

  shocked_fit <- dly_fit_from_single_curve(
    dly_sc = shocked_single_curve,
    curves = dly_fit_obj$curves,
    reference = dly_fit_obj$reference,
    center = dly_fit_obj$center,
    scale. = dly_fit_obj$scale.,
    mats = dly_fit_obj$mats
  )

  return(list(
    fit = shocked_fit,
    shock_matrix = shock_fit$shock_matrix,
    shocked_factors = shock_fit$shocked_factors,
    baseline_fit = dly_fit_obj
  ))
}

st_align_Xt_to_baseline <- function(Xt_new, Xt_reference) {
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

st_finalize_shocked_fit <- function(yields, time, phi_hat, Xt, BS0, curves, mats, W) {
  Y <- vecY(yields, curves, mats)
  P <- bln_build_P(curves, mats)
  pY <- PY_full(Y$tY, P)

  sigma_JT <- tSigma_optim(BS0 = BS0, X = Xt)
  full_sigma_t <- sqrt_inv_build(sigma_JT$tSigma)

  sqrt_sigma_t <- full_sigma_t$sqrt
  sqrt_inv_sigma_t <- full_sigma_t$inverse

  yb <- bln_YB(pY$py, sqrt_inv_sigma_t)
  phib <- bln_phi_BT(phi_hat, P, sqrt_inv_sigma_t)
  phic <- bln_phi_CT(phib, P)
  yc <- PY_full(yb, t(P))

  ycut <- cf_curve_cut(yc$py, curves, mats)
  phicut <- cf_phi_cut(phic, curves, mats)
  normalized <- bln_normalized_fit(ycut, phicut, curves, mats)

  shocked_tenor <- blce_rescale_yield(
    normalized$n_yield,
    sqrt_sigma_t,
    P,
    "tenor",
    time = time,
    curves = curves,
    mats = mats
  )
  shocked_curve <- blce_rescale_yield(
    normalized$n_yield,
    sqrt_sigma_t,
    P,
    "curve",
    time = time,
    curves = curves,
    mats = mats
  )

  return(list(
    tenor = shocked_tenor,
    curve = shocked_curve,
    sigma_JT = sigma_JT,
    Xt = Xt,
    W = W
  ))
}

st_resolve_config <- function(value, fallback) {
  if (is.null(value)) {
    fallback
  } else {
    value
  }
}

st_validate_baseline_state <- function(baseline_state, curves = NULL, mats = NULL, reference = NULL) {
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

  baseline_state
}

mc_fit_shock_baseline <- function(yields, lambdas, cutoffs, reference,
                                  ECM_estim = "ML", ECM_type = "eigen", ECM_alpha = 0.1,
                                  X_normalize = TRUE, X_trunc = FALSE,
                                  CR_algo = zero_mean_covreg_em, CR_init = "adaptive",
                                  Xt_smooth = FALSE, smoother = c("ns", "bs", "rm", "hp", "henderson"),
                                  knot_count = 5, k_count = 9, k_pass = 3,
                                  hp_lambda = 1600, henderson_k = 13,
                                  CR_maxiter = 1000, CR_tol = 1e-8,
                                  CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = 0,
                                  CR_verb = FALSE, CR_term = TRUE,
                                  curves, mats) {
  yields <- st_validate_yields_input(yields, curves, reference)

  single_curve_state <- st_prepare_single_curve_state(
    yields = yields,
    lambdas = lambdas,
    cutoffs = cutoffs,
    curves = curves,
    mats = mats
  )

  baseline_xt_state <- st_build_xt_state(
    bl_yields = single_curve_state$bl_yields,
    curves = curves,
    reference = reference,
    reuse_ECM = FALSE,
    baseline_ecm = NULL,
    ECM_estim = ECM_estim,
    ECM_type = ECM_type,
    ECM_alpha = ECM_alpha,
    X_normalize = X_normalize,
    X_trunc = X_trunc,
    Xt_smooth = Xt_smooth,
    smoother = smoother,
    knot_count = knot_count,
    k_count = k_count,
    k_pass = k_pass,
    hp_lambda = hp_lambda,
    henderson_k = henderson_k,
    ytime = single_curve_state$time
  )

  baseline_BS0 <- blcc_fe(
    Xt = baseline_xt_state$Xt,
    Wt = single_curve_state$W,
    mats = mats,
    covreg = CR_algo,
    init = CR_init,
    max_iter = CR_maxiter,
    tol = CR_tol,
    S0 = CR_S0init,
    B = CR_Binit,
    S0_shrink_diag = CR_S0_shrink_diag,
    verb = CR_verb,
    term = CR_term
  )

  baseline_fit_state <- st_finalize_shocked_fit(
    yields = single_curve_state$yields,
    time = single_curve_state$time,
    phi_hat = single_curve_state$phi_hat,
    Xt = baseline_xt_state$Xt,
    BS0 = baseline_BS0,
    curves = curves,
    mats = mats,
    W = single_curve_state$W
  )

  structure(list(
    curves = curves,
    mats = mats,
    reference = reference,
    yields = yields,
    single_curve_state = single_curve_state,
    xt_state = baseline_xt_state,
    BS0 = baseline_BS0,
    baseline_fit = baseline_fit_state,
    config = list(
      lambdas = lambdas,
      cutoffs = cutoffs,
      ECM_estim = ECM_estim,
      ECM_type = ECM_type,
      ECM_alpha = ECM_alpha,
      X_normalize = X_normalize,
      X_trunc = X_trunc,
      CR_algo = CR_algo,
      CR_init = CR_init,
      Xt_smooth = Xt_smooth,
      smoother = smoother,
      knot_count = knot_count,
      k_count = k_count,
      k_pass = k_pass,
      hp_lambda = hp_lambda,
      henderson_k = henderson_k,
      CR_maxiter = CR_maxiter,
      CR_tol = CR_tol,
      CR_Binit = CR_Binit,
      CR_S0init = CR_S0init,
      CR_S0_shrink_diag = CR_S0_shrink_diag,
      CR_verb = CR_verb,
      CR_term = CR_term
    )
  ), class = "mc_shock_baseline")
}

mc_fit_shock_from_baseline <- function(baseline_state,
                                       shock_curves, shock_factors, shock_magnitude,
                                       shock_type = c("additive", "multiplicative"),
                                       shock_window = NULL,
                                       reuse_BS0 = TRUE, reuse_ECM = TRUE,
                                       ECM_estim = NULL, ECM_type = NULL, ECM_alpha = NULL,
                                       X_normalize = NULL, X_trunc = NULL,
                                       CR_algo = NULL, CR_init = NULL,
                                       Xt_smooth = NULL, smoother = NULL,
                                       knot_count = NULL, k_count = NULL, k_pass = NULL,
                                       hp_lambda = NULL, henderson_k = NULL,
                                       CR_maxiter = NULL, CR_tol = NULL,
                                       CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = NULL,
                                       CR_verb = NULL, CR_term = NULL) {
  shock_type <- match.arg(shock_type)
  baseline_state <- st_validate_baseline_state(baseline_state)
  config <- baseline_state$config
  single_curve_state <- baseline_state$single_curve_state

  ECM_estim <- st_resolve_config(ECM_estim, config$ECM_estim)
  ECM_type <- st_resolve_config(ECM_type, config$ECM_type)
  ECM_alpha <- st_resolve_config(ECM_alpha, config$ECM_alpha)
  X_normalize <- st_resolve_config(X_normalize, config$X_normalize)
  X_trunc <- st_resolve_config(X_trunc, config$X_trunc)
  CR_algo <- st_resolve_config(CR_algo, config$CR_algo)
  CR_init <- st_resolve_config(CR_init, config$CR_init)
  Xt_smooth <- st_resolve_config(Xt_smooth, config$Xt_smooth)
  smoother <- st_resolve_config(smoother, config$smoother)
  knot_count <- st_resolve_config(knot_count, config$knot_count)
  k_count <- st_resolve_config(k_count, config$k_count)
  k_pass <- st_resolve_config(k_pass, config$k_pass)
  hp_lambda <- st_resolve_config(hp_lambda, config$hp_lambda)
  henderson_k <- st_resolve_config(henderson_k, config$henderson_k)
  CR_maxiter <- st_resolve_config(CR_maxiter, config$CR_maxiter)
  CR_tol <- st_resolve_config(CR_tol, config$CR_tol)
  CR_Binit <- st_resolve_config(CR_Binit, config$CR_Binit)
  CR_S0init <- st_resolve_config(CR_S0init, config$CR_S0init)
  CR_S0_shrink_diag <- st_resolve_config(CR_S0_shrink_diag, config$CR_S0_shrink_diag)
  CR_verb <- st_resolve_config(CR_verb, config$CR_verb)
  CR_term <- st_resolve_config(CR_term, config$CR_term)

  shock_fit <- st_apply_factor_shock(
    curves = baseline_state$curves,
    curve_dns_factor = single_curve_state$ns_factor,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window
  )

  shocked_state <- st_build_shocked_curve_state(
    single_curve_state = single_curve_state,
    shocked_factors = shock_fit$shocked_factors
  )

  shocked_xt_state <- st_build_xt_state(
    bl_yields = shocked_state$bl_yields,
    curves = baseline_state$curves,
    reference = baseline_state$reference,
    reuse_ECM = reuse_ECM,
    baseline_ecm = if (reuse_ECM) baseline_state$xt_state$cc_ci_features$cc_ecm else NULL,
    ECM_estim = ECM_estim,
    ECM_type = ECM_type,
    ECM_alpha = ECM_alpha,
    X_normalize = X_normalize,
    X_trunc = X_trunc,
    Xt_smooth = Xt_smooth,
    smoother = smoother,
    knot_count = knot_count,
    k_count = k_count,
    k_pass = k_pass,
    hp_lambda = hp_lambda,
    henderson_k = henderson_k,
    ytime = single_curve_state$time
  )

  if (reuse_BS0) {
    shocked_Xt <- st_align_Xt_to_baseline(shocked_xt_state$Xt, baseline_state$xt_state$Xt)
    shocked_BS0 <- baseline_state$BS0
  } else {
    shocked_Xt <- shocked_xt_state$Xt
    shocked_BS0 <- blcc_fe(
      Xt = shocked_Xt,
      Wt = shocked_state$W,
      mats = baseline_state$mats,
      covreg = CR_algo,
      init = CR_init,
      max_iter = CR_maxiter,
      tol = CR_tol,
      S0 = CR_S0init,
      B = CR_Binit,
      S0_shrink_diag = CR_S0_shrink_diag,
      verb = CR_verb,
      term = CR_term
    )
  }

  shocked_fit_state <- st_finalize_shocked_fit(
    yields = shocked_state$yields,
    time = single_curve_state$time,
    phi_hat = shocked_state$phi_hat,
    Xt = shocked_Xt,
    BS0 = shocked_BS0,
    curves = baseline_state$curves,
    mats = baseline_state$mats,
    W = shocked_state$W
  )

  return(list(
    tenor = baseline_state$baseline_fit$tenor,
    curve = baseline_state$baseline_fit$curve,
    sigma_JT = shocked_fit_state$sigma_JT,
    Xt = shocked_fit_state$Xt,
    W = shocked_fit_state$W,
    ns_factor = shocked_state$ns_factor,
    shock_matrix = shock_fit$shock_matrix,
    shocked_factors = shock_fit$shocked_factors,
    shocked_tenor = shocked_fit_state$tenor,
    shocked_curve = shocked_fit_state$curve,
    baseline_state = baseline_state
  ))
}

mc_fit_shock <- function(yields = NULL, lambdas = NULL, cutoffs = NULL, reference = NULL,
                         shock_curves, shock_factors, shock_magnitude,
                         shock_type = c("additive", "multiplicative"),
                         shock_window = NULL,
                         reuse_BS0 = TRUE, reuse_ECM = TRUE,
                         ECM_estim = "ML", ECM_type = "eigen", ECM_alpha = 0.1,
                         X_normalize = TRUE, X_trunc = FALSE,
                         CR_algo = zero_mean_covreg_em, CR_init = "adaptive",
                         Xt_smooth = FALSE, smoother = c("ns", "bs", "rm", "hp", "henderson"),
                         knot_count = 5, k_count = 9, k_pass = 3,
                         hp_lambda = 1600, henderson_k = 13,
                         CR_maxiter = 1000, CR_tol = 1e-8,
                         CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = 0,
                         CR_verb = FALSE, CR_term = TRUE,
                         curves = NULL, mats = NULL, baseline_state = NULL) {
  if (is.null(baseline_state)) {
    if (is.null(yields) || is.null(lambdas) || is.null(cutoffs) ||
        is.null(reference) || is.null(curves) || is.null(mats)) {
      stop(
        "When `baseline_state` is not supplied, `yields`, `lambdas`, `cutoffs`, ",
        "`reference`, `curves`, and `mats` must be provided."
      )
    }

    baseline_state <- mc_fit_shock_baseline(
      yields = yields,
      lambdas = lambdas,
      cutoffs = cutoffs,
      reference = reference,
      ECM_estim = ECM_estim,
      ECM_type = ECM_type,
      ECM_alpha = ECM_alpha,
      X_normalize = X_normalize,
      X_trunc = X_trunc,
      CR_algo = CR_algo,
      CR_init = CR_init,
      Xt_smooth = Xt_smooth,
      smoother = smoother,
      knot_count = knot_count,
      k_count = k_count,
      k_pass = k_pass,
      hp_lambda = hp_lambda,
      henderson_k = henderson_k,
      CR_maxiter = CR_maxiter,
      CR_tol = CR_tol,
      CR_Binit = CR_Binit,
      CR_S0init = CR_S0init,
      CR_S0_shrink_diag = CR_S0_shrink_diag,
      CR_verb = CR_verb,
      CR_term = CR_term,
      curves = curves,
      mats = mats
    )
  } else {
    baseline_state <- st_validate_baseline_state(
      baseline_state = baseline_state,
      curves = curves,
      mats = mats,
      reference = reference
    )
  }

  mc_fit_shock_from_baseline(
    baseline_state = baseline_state,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window,
    reuse_BS0 = reuse_BS0,
    reuse_ECM = reuse_ECM,
    ECM_estim = ECM_estim,
    ECM_type = ECM_type,
    ECM_alpha = ECM_alpha,
    X_normalize = X_normalize,
    X_trunc = X_trunc,
    CR_algo = CR_algo,
    CR_init = CR_init,
    Xt_smooth = Xt_smooth,
    smoother = smoother,
    knot_count = knot_count,
    k_count = k_count,
    k_pass = k_pass,
    hp_lambda = hp_lambda,
    henderson_k = henderson_k,
    CR_maxiter = CR_maxiter,
    CR_tol = CR_tol,
    CR_Binit = CR_Binit,
    CR_S0init = CR_S0init,
    CR_S0_shrink_diag = CR_S0_shrink_diag,
    CR_verb = CR_verb,
    CR_term = CR_term
  )
}
