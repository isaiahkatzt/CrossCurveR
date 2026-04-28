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
                              ECM_fallback = c("error", "baseline"),
                              X_normalize = TRUE, X_trunc = FALSE,
                              Xt_smooth = FALSE, smoother = c("ns", "bs", "rm", "hp", "henderson"),
                              knot_count = 5, k_count = 9, k_pass = 3,
                              hp_lambda = 1600, henderson_k = 13, ytime) {
  ECM_fallback <- match.arg(ECM_fallback)

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
    xt_state <- tryCatch(
      list(
        Xt = blcc_build_Xt(cc_ci_features$cc_cspread, cc_ci_features$cc_ecm, trunc = X_trunc),
        cc_ci_features = cc_ci_features
      ),
      error = function(e) {
        if (ECM_fallback != "baseline" || is.null(baseline_ecm)) {
          stop(e)
        }

        warning(
          "Re-estimated ECM did not produce usable Xt features; falling back to baseline ECM.",
          call. = FALSE
        )

        fallback_state <- st_build_xt_from_ecm(
          bl_yields = bl_yields,
          curves = curves,
          reference = reference,
          curve_ecm = baseline_ecm,
          X_normalize = X_normalize,
          X_trunc = X_trunc
        )
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
                                  shock_window = NULL,
                                  shock_window_position = c("first", "last")) {
  shock_type <- match.arg(shock_type)
  shock_window_position <- match.arg(shock_window_position)

  switch(
    shock_type,
    additive = st_additive(
      curves = curves,
      curve_dns_factor = curve_dns_factor,
      shock_curves = shock_curves,
      shock_factors = shock_factors,
      shock_magnitude = shock_magnitude,
      shock_window = shock_window,
      shock_window_position = shock_window_position
    ),
    multiplicative = st_multiplicative(
      curves = curves,
      curve_dns_factor = curve_dns_factor,
      shock_curves = shock_curves,
      shock_factors = shock_factors,
      shock_magnitude = shock_magnitude,
      shock_window = shock_window,
      shock_window_position = shock_window_position
    )
  )
}

st_dns_shock_from_single_curve_state <- function(single_curve_state,
                                                 shock_curves, shock_factors, shock_magnitude,
                                                 shock_type = c("additive", "multiplicative"),
                                                 shock_window = NULL,
                                                 shock_window_position = c("first", "last")) {
  shock_fit <- st_apply_factor_shock(
    curves = single_curve_state$curves,
    curve_dns_factor = single_curve_state$ns_factor,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window,
    shock_window_position = shock_window_position
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

st_dly_curve_yields_from_single_curve <- function(dly_single_curve,
                                                  curves = names(dly_single_curve)) {
  setNames(lapply(curves, function(curve_name) {
    curve_state <- dly_single_curve[[curve_name]]
    curve_df <- data.frame(
      time = curve_state$time,
      as.data.frame(curve_state$betas$m_yield),
      check.names = FALSE
    )
    colnames(curve_df) <- c("time", curve_state$yield_names)
    curve_df
  }), curves)
}

st_dly_shock_refit_core <- function(dly_fit_obj,
                                    shock_curves, shock_factors, shock_magnitude,
                                    shock_type, shock_window = NULL,
                                    shock_window_position = c("first", "last"),
                                    lambdas = NULL) {
  lambda_grid <- lambdas
  if (is.null(lambda_grid) && !is.null(dly_fit_obj$lambda_grid)) {
    lambda_grid <- dly_fit_obj$lambda_grid
  }

  if (is.null(lambda_grid)) {
    stop("`lambdas` must be supplied unless `dly_fit_obj` stores `lambda_grid`.")
  }

  shock_fit <- st_apply_factor_shock(
    curves = dly_fit_obj$curves,
    curve_dns_factor = dly_pull_beta(dly_fit_obj$single_curve, curves = dly_fit_obj$curves),
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window,
    shock_window_position = shock_window_position
  )

  shocked_single_curve <- st_build_shocked_dly_single_curve(
    dly_single_curve = dly_fit_obj$single_curve,
    shocked_factors = shock_fit$shocked_factors,
    curves = dly_fit_obj$curves
  )
  shocked_yields <- st_dly_curve_yields_from_single_curve(
    dly_single_curve = shocked_single_curve,
    curves = dly_fit_obj$curves
  )

  shocked_fit <- dly_fit(
    yields = shocked_yields,
    mats = dly_fit_obj$mats,
    lambdas = lambda_grid,
    curves = dly_fit_obj$curves,
    reference = dly_fit_obj$reference,
    center = dly_fit_obj$center,
    scale. = dly_fit_obj$scale.
  )

  return(list(
    fit = shocked_fit,
    yields = shocked_yields,
    shock_matrix = shock_fit$shock_matrix,
    shocked_factors = shock_fit$shocked_factors,
    shocked_single_curve = shocked_single_curve,
    baseline_fit = dly_fit_obj
  ))
}

st_dly_shock_refit_from_fit <- function(dly_fit_obj,
                                        shock_curves, shock_factors, shock_magnitude,
                                        shock_type = c("additive", "multiplicative"),
                                        shock_window = NULL,
                                        shock_window_position = c("first", "last"),
                                        lambdas = NULL) {
  shock_type <- match.arg(shock_type)
  shock_window_position <- match.arg(shock_window_position)

  st_dly_shock_from_fit(
    dly_fit_obj = dly_fit_obj,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window,
    shock_window_position = shock_window_position,
    refit_dly = TRUE,
    lambdas = lambdas
  )
}

st_dly_project_global_factor <- function(factor_draw, global_factor_obj) {
  factor_draw <- as.matrix(factor_draw)
  pca_fit <- global_factor_obj$pca

  if (is.null(pca_fit) || is.null(pca_fit$rotation)) {
    stop("Baseline DLY fit must include a PCA object for frozen shock propagation.")
  }

  expected_cols <- rownames(pca_fit$rotation)

  if (is.null(colnames(factor_draw)) || !setequal(colnames(factor_draw), expected_cols)) {
    stop("Shocked factor draw columns must match the baseline DLY PCA variables.")
  }

  factor_draw <- factor_draw[, expected_cols, drop = FALSE]
  projection <- stats::predict(pca_fit, newdata = factor_draw)
  as.numeric(projection[, 1])
}

st_dly_predict_country_factor <- function(factor_model, global_factor, curves) {
  intercept <- setNames(as.numeric(factor_model$intercept[curves]), curves)
  loading <- setNames(as.numeric(factor_model$loading[curves]), curves)

  fitted <- outer(as.numeric(global_factor), loading)
  fitted <- sweep(fitted, 2, intercept, FUN = "+")
  colnames(fitted) <- curves

  template <- factor_model$fitted
  if (!is.null(template) && !is.null(rownames(template))) {
    rownames(fitted) <- rownames(template)
  } else if (!is.null(factor_model$residual) && !is.null(rownames(factor_model$residual))) {
    rownames(fitted) <- rownames(factor_model$residual)
  }

  fitted
}

st_dly_shock_residual_states <- function(raw_level_factor, raw_slope_factor,
                                         fitted_level_factor, fitted_slope_factor,
                                         curves) {
  setNames(lapply(curves, function(curve_name) {
    out <- cbind(
      L = raw_level_factor[, curve_name] - fitted_level_factor[, curve_name],
      S = raw_slope_factor[, curve_name] - fitted_slope_factor[, curve_name]
    )
    rownames(out) <- rownames(fitted_level_factor)
    out
  }), curves)
}

st_dly_shock_global_path <- function(dly_fit_obj, shocked_global_level, shocked_global_slope,
                                     shock_rows) {
  global_path <- as.matrix(dly_fit_obj$dynamics$global$propagated)
  override_state <- cbind(L = shocked_global_level, S = shocked_global_slope)
  rownames(override_state) <- rownames(global_path)

  if (1 %in% shock_rows) {
    global_path[1, ] <- override_state[1, ]
  }

  if (nrow(global_path) >= 2) {
    for (t in 2:nrow(global_path)) {
      if (t %in% shock_rows) {
        global_path[t, ] <- override_state[t, ]
      } else {
        global_path[t, ] <- dly_transition_var1(
          prev_state = global_path[t - 1, ],
          intercept = dly_fit_obj$dynamics$global$intercept,
          phi = dly_fit_obj$dynamics$global$phi
        )
      }
    }
  }

  global_path
}

st_dly_shock_residual_path <- function(dly_fit_obj, global_path, shock_residual_states,
                                       shock_rows) {
  curves <- dly_fit_obj$curves

  setNames(lapply(curves, function(curve_name) {
    residual_model <- dly_fit_obj$dynamics$residual[[curve_name]]
    residual_path <- as.matrix(residual_model$propagated)
    override_state <- as.matrix(shock_residual_states[[curve_name]])

    if (1 %in% shock_rows) {
      residual_path[1, ] <- override_state[1, ]
    }

    if (nrow(residual_path) >= 2) {
      for (t in 2:nrow(residual_path)) {
        if (t %in% shock_rows) {
          residual_path[t, ] <- override_state[t, ]
        } else {
          residual_path[t, ] <- dly_transition_country_residual(
            prev_residual = residual_path[t - 1, ],
            prev_global = global_path[t - 1, ],
            intercept = residual_model$intercept,
            phi = residual_model$phi,
            kappa = residual_model$kappa
          )
        }
      }
    }

    residual_path
  }), curves)
}

st_dly_shock_from_fit <- function(dly_fit_obj,
                                  shock_curves, shock_factors, shock_magnitude,
                                  shock_type = c("additive", "multiplicative"),
                                  shock_window = NULL,
                                  shock_window_position = c("first", "last"),
                                  refit_dly = FALSE,
                                  lambdas = NULL) {
  shock_type <- match.arg(shock_type)
  shock_window_position <- match.arg(shock_window_position)

  if (isTRUE(refit_dly)) {
    return(st_dly_shock_refit_core(
      dly_fit_obj = dly_fit_obj,
      shock_curves = shock_curves,
      shock_factors = shock_factors,
      shock_magnitude = shock_magnitude,
      shock_type = shock_type,
      shock_window = shock_window,
      shock_window_position = shock_window_position,
      lambdas = lambdas
    ))
  }

  shock_fit <- st_apply_factor_shock(
    curves = dly_fit_obj$curves,
    curve_dns_factor = dly_pull_beta(dly_fit_obj$single_curve, curves = dly_fit_obj$curves),
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window,
    shock_window_position = shock_window_position
  )

  shocked_single_curve <- st_build_shocked_dly_single_curve(
    dly_single_curve = dly_fit_obj$single_curve,
    shocked_factors = shock_fit$shocked_factors,
    curves = dly_fit_obj$curves
  )

  raw_level_factor <- dly_group_factor(shocked_single_curve, curves = dly_fit_obj$curves, latent = "L")
  raw_slope_factor <- dly_group_factor(shocked_single_curve, curves = dly_fit_obj$curves, latent = "S")

  shocked_global_level <- st_dly_project_global_factor(
    factor_draw = raw_level_factor,
    global_factor_obj = dly_fit_obj$global_factors$level
  )
  shocked_global_slope <- st_dly_project_global_factor(
    factor_draw = raw_slope_factor,
    global_factor_obj = dly_fit_obj$global_factors$slope
  )

  shocked_level_fitted <- st_dly_predict_country_factor(
    factor_model = dly_fit_obj$factor_models$level,
    global_factor = shocked_global_level,
    curves = dly_fit_obj$curves
  )
  shocked_slope_fitted <- st_dly_predict_country_factor(
    factor_model = dly_fit_obj$factor_models$slope,
    global_factor = shocked_global_slope,
    curves = dly_fit_obj$curves
  )

  baseline_level_residual <- as.matrix(
    dly_fit_obj$factor_models$level$residual[, dly_fit_obj$curves, drop = FALSE]
  )
  baseline_slope_residual <- as.matrix(
    dly_fit_obj$factor_models$slope$residual[, dly_fit_obj$curves, drop = FALSE]
  )

  shock_rows <- st_resolve_shock_rows(
    row_count = nrow(raw_level_factor),
    shock_window = shock_window,
    shock_window_position = shock_window_position
  )

  shock_residual_states <- st_dly_shock_residual_states(
    raw_level_factor = raw_level_factor,
    raw_slope_factor = raw_slope_factor,
    fitted_level_factor = shocked_level_fitted,
    fitted_slope_factor = shocked_slope_fitted,
    curves = dly_fit_obj$curves
  )
  propagated_global_state <- st_dly_shock_global_path(
    dly_fit_obj = dly_fit_obj,
    shocked_global_level = shocked_global_level,
    shocked_global_slope = shocked_global_slope,
    shock_rows = shock_rows
  )
  propagated_residual_state <- st_dly_shock_residual_path(
    dly_fit_obj = dly_fit_obj,
    global_path = propagated_global_state,
    shock_residual_states = shock_residual_states,
    shock_rows = shock_rows
  )
  propagated_global_factor <- dly_country_global_factors(
    level_global_factor = propagated_global_state[, "L"],
    slope_global_factor = propagated_global_state[, "S"],
    factor_models = dly_fit_obj$factor_models,
    curves = dly_fit_obj$curves
  )
  propagated_residual_factor <- dly_unpack_country_state(
    country_state_list = propagated_residual_state,
    curves = dly_fit_obj$curves
  )
  propagated_level_factor <- propagated_global_factor$level + propagated_residual_factor$level
  propagated_slope_factor <- propagated_global_factor$slope + propagated_residual_factor$slope

  raw_betas <- dly_bind_betas(raw_level_factor, raw_slope_factor, curves = dly_fit_obj$curves)
  global_betas <- dly_bind_betas(shocked_level_fitted, shocked_slope_fitted, curves = dly_fit_obj$curves)
  idio_betas <- dly_bind_betas(baseline_level_residual, baseline_slope_residual, curves = dly_fit_obj$curves)
  propagated_betas <- dly_bind_betas(
    propagated_level_factor,
    propagated_slope_factor,
    curves = dly_fit_obj$curves
  )

  raw_yields <- dly_curve_yields(
    shocked_single_curve,
    raw_level_factor,
    raw_slope_factor,
    curves = dly_fit_obj$curves
  )
  global_yields <- dly_curve_yields(
    shocked_single_curve,
    shocked_level_fitted,
    shocked_slope_fitted,
    curves = dly_fit_obj$curves
  )
  idiosyncratic_yields <- dly_curve_yields(
    shocked_single_curve,
    baseline_level_residual,
    baseline_slope_residual,
    curves = dly_fit_obj$curves
  )
  propagated_yields <- dly_curve_yields(
    shocked_single_curve,
    propagated_level_factor,
    propagated_slope_factor,
    curves = dly_fit_obj$curves
  )

  shocked_fit <- list(
    curves = dly_fit_obj$curves,
    mats = dly_fit_obj$mats,
    reference = dly_fit_obj$reference,
    center = dly_fit_obj$center,
    scale. = dly_fit_obj$scale.,
    lambdas = dly_fit_obj$lambdas,
    single_curve = shocked_single_curve,
    country_factors = list(
      level = raw_level_factor,
      slope = raw_slope_factor
    ),
    global_factors = list(
      level = list(
        factor = shocked_global_level,
        loading = dly_fit_obj$global_factors$level$loading,
        pca = dly_fit_obj$global_factors$level$pca,
        explained_variance = dly_fit_obj$global_factors$level$explained_variance
      ),
      slope = list(
        factor = shocked_global_slope,
        loading = dly_fit_obj$global_factors$slope$loading,
        pca = dly_fit_obj$global_factors$slope$pca,
        explained_variance = dly_fit_obj$global_factors$slope$explained_variance
      )
    ),
    factor_models = list(
      level = list(
        intercept = dly_fit_obj$factor_models$level$intercept,
        loading = dly_fit_obj$factor_models$level$loading,
        fitted = shocked_level_fitted,
        residual = baseline_level_residual,
        r_squared = dly_fit_obj$factor_models$level$r_squared
      ),
      slope = list(
        intercept = dly_fit_obj$factor_models$slope$intercept,
        loading = dly_fit_obj$factor_models$slope$loading,
        fitted = shocked_slope_fitted,
        residual = baseline_slope_residual,
        r_squared = dly_fit_obj$factor_models$slope$r_squared
      )
    ),
    country_betas = list(
      raw = raw_betas,
      global = global_betas,
      idiosyncratic = idio_betas,
      dislocation_propagated = propagated_betas
    ),
    yields = list(
      raw = raw_yields,
      global = global_yields,
      idiosyncratic = idiosyncratic_yields,
      dislocation_propagated = propagated_yields
    ),
    residuals = list(
      raw = dly_curve_residuals(raw_yields, raw_yields, dly_fit_obj$curves),
      global = dly_curve_residuals(raw_yields, global_yields, dly_fit_obj$curves),
      dislocation_propagated = dly_curve_residuals(
        raw_yields,
        propagated_yields,
        dly_fit_obj$curves
      )
    ),
    dynamics = list(
      global = c(
        dly_fit_obj$dynamics$global[names(dly_fit_obj$dynamics$global) %in% c("intercept", "phi", "fitted", "residual")],
        list(
          actual = as.matrix(dly_fit_obj$dynamics$global$actual),
          propagated = propagated_global_state
        )
      ),
      residual = setNames(lapply(dly_fit_obj$curves, function(curve_name) {
        c(
          dly_fit_obj$dynamics$residual[[curve_name]][names(dly_fit_obj$dynamics$residual[[curve_name]]) %in% c("intercept", "phi", "kappa", "fitted", "residual")],
          list(
            actual = as.matrix(dly_fit_obj$dynamics$residual[[curve_name]]$actual),
            propagated = propagated_residual_state[[curve_name]]
          )
        )
      }), dly_fit_obj$curves)
    )
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

st_curve_shock_delta <- function(baseline_yields, shocked_yields, shock_curves) {
  shock_curves <- intersect(as.character(shock_curves), names(shocked_yields))

  setNames(lapply(shock_curves, function(curve_name) {
    baseline_curve <- baseline_yields[[curve_name]] %>%
      dplyr::mutate(time = as.Date(time)) %>%
      dplyr::arrange(time)
    shocked_curve <- shocked_yields[[curve_name]] %>%
      dplyr::mutate(time = as.Date(time)) %>%
      dplyr::arrange(time)

    if (!identical(as.Date(baseline_curve$time), as.Date(shocked_curve$time))) {
      stop("Baseline and shocked curve histories must align to compute anchored shock deltas.")
    }

    delta_curve <- data.frame(
      time = baseline_curve$time,
      shocked_curve[, -1, drop = FALSE] - baseline_curve[, -1, drop = FALSE],
      check.names = FALSE
    )
    colnames(delta_curve) <- colnames(baseline_curve)
    delta_curve
  }), shock_curves)
}

st_anchor_curve_output <- function(curve_panel, baseline_curve_panel, shock_delta_yields,
                                   shock_curves, curves, mats, anchor_rows = NULL) {
  shock_curves <- intersect(as.character(shock_curves), as.character(curves))

  if (length(shock_curves) == 0) {
    return(curve_panel)
  }

  mat_str <- sapply(mats, numeric_to_matname)
  anchored_panel <- curve_panel
  panel_time <- as.Date(curve_panel$time)
  panel_rows <- seq_len(nrow(curve_panel))

  if (is.null(anchor_rows)) {
    anchor_rows <- panel_rows
  } else {
    if (!is.numeric(anchor_rows) || anyNA(anchor_rows) ||
        any(anchor_rows %% 1 != 0)) {
      stop("`anchor_rows` must be NULL or an integer row index vector.")
    }
    anchor_rows <- as.integer(anchor_rows)
    anchor_rows <- anchor_rows[anchor_rows >= 1 & anchor_rows <= nrow(curve_panel)]
  }

  if (length(anchor_rows) == 0) {
    return(anchored_panel)
  }

  for (curve_name in shock_curves) {
    shock_delta_curve <- shock_delta_yields[[curve_name]]
    if (is.null(shock_delta_curve)) {
      stop("Each anchored shock curve must be present in `shock_delta_yields`.")
    }

    baseline_curve <- baseline_curve_panel %>%
      dplyr::mutate(time = as.Date(time)) %>%
      dplyr::arrange(time)
    shock_delta_curve <- shock_delta_curve %>%
      dplyr::mutate(time = as.Date(time)) %>%
      dplyr::arrange(time)

    row_idx <- match(panel_time, baseline_curve$time)
    if (anyNA(row_idx)) {
      stop("Baseline curve panel must align with the finalized shocked output time index.")
    }

    delta_idx <- match(panel_time, shock_delta_curve$time)
    if (anyNA(delta_idx)) {
      stop("Shock deltas must align with the finalized shocked output time index.")
    }

    curve_cols <- paste0(curve_name, ".", mat_str)
    anchored_panel[anchor_rows, curve_cols] <-
      baseline_curve[row_idx[anchor_rows], curve_cols, drop = FALSE] +
      shock_delta_curve[delta_idx[anchor_rows], colnames(shock_delta_curve)[-1], drop = FALSE]
  }

  anchored_panel
}

st_tenor_panel_from_curve_panel <- function(curve_panel, curves, mats) {
  P <- bln_build_P(curves, mats)
  tenor_panel <- PY_full(curve_panel[, -1, drop = FALSE], P)$py
  data.frame(time = as.Date(curve_panel$time), tenor_panel, check.names = FALSE)
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
                                  CR_check_every = 1, CR_use_loglik = TRUE, CR_store_path = TRUE,
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
    check_every = CR_check_every,
    use_loglik = CR_use_loglik,
    store_path = CR_store_path,
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
      CR_check_every = CR_check_every,
      CR_use_loglik = CR_use_loglik,
      CR_store_path = CR_store_path,
      CR_verb = CR_verb,
      CR_term = CR_term
    )
  ), class = "mc_shock_baseline")
}

mc_fit_shock_from_baseline <- function(baseline_state,
                                       shock_curves, shock_factors, shock_magnitude,
                                       shock_type = c("additive", "multiplicative"),
                                       shock_window = NULL,
                                       shock_window_position = c("first", "last"),
                                       reuse_BS0 = TRUE, reuse_ECM = TRUE,
                                       anchor_shock_curves = FALSE,
                                       ECM_estim = NULL, ECM_type = NULL, ECM_alpha = NULL,
                                       ECM_fallback = c("error", "baseline"),
                                       X_normalize = NULL, X_trunc = NULL,
                                       CR_algo = NULL, CR_init = NULL,
                                       Xt_smooth = NULL, smoother = NULL,
                                       knot_count = NULL, k_count = NULL, k_pass = NULL,
                                       hp_lambda = NULL, henderson_k = NULL,
                                       CR_maxiter = NULL, CR_tol = NULL,
                                       CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = NULL,
                                       CR_check_every = NULL, CR_use_loglik = NULL, CR_store_path = NULL,
                                       CR_verb = NULL, CR_term = NULL) {
  shock_type <- match.arg(shock_type)
  shock_window_position <- match.arg(shock_window_position)
  ECM_fallback <- match.arg(ECM_fallback)
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
  CR_check_every <- st_resolve_config(CR_check_every, config$CR_check_every)
  CR_use_loglik <- st_resolve_config(CR_use_loglik, config$CR_use_loglik)
  CR_store_path <- st_resolve_config(CR_store_path, config$CR_store_path)
  CR_verb <- st_resolve_config(CR_verb, config$CR_verb)
  CR_term <- st_resolve_config(CR_term, config$CR_term)

  shock_fit <- st_apply_factor_shock(
    curves = baseline_state$curves,
    curve_dns_factor = single_curve_state$ns_factor,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window,
    shock_window_position = shock_window_position
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
    baseline_ecm = baseline_state$xt_state$cc_ci_features$cc_ecm,
    ECM_estim = ECM_estim,
    ECM_type = ECM_type,
    ECM_alpha = ECM_alpha,
    ECM_fallback = ECM_fallback,
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
      check_every = CR_check_every,
      use_loglik = CR_use_loglik,
      store_path = CR_store_path,
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

  if (isTRUE(anchor_shock_curves)) {
    anchor_rows <- st_resolve_shock_rows(
      row_count = nrow(single_curve_state$yields[[1]]),
      shock_window = shock_window,
      shock_window_position = shock_window_position
    )
    shock_delta_yields <- st_curve_shock_delta(
      baseline_yields = single_curve_state$yields,
      shocked_yields = shocked_state$yields,
      shock_curves = shock_curves
    )
    shocked_fit_state$curve <- st_anchor_curve_output(
      curve_panel = shocked_fit_state$curve,
      baseline_curve_panel = baseline_state$baseline_fit$curve,
      shock_delta_yields = shock_delta_yields,
      shock_curves = shock_curves,
      curves = baseline_state$curves,
      mats = baseline_state$mats,
      anchor_rows = anchor_rows
    )
    shocked_fit_state$tenor <- st_tenor_panel_from_curve_panel(
      curve_panel = shocked_fit_state$curve,
      curves = baseline_state$curves,
      mats = baseline_state$mats
    )
  }

  return(list(
    tenor = baseline_state$baseline_fit$tenor,
    curve = baseline_state$baseline_fit$curve,
    sigma_JT = shocked_fit_state$sigma_JT,
    Xt = shocked_fit_state$Xt,
    W = shocked_fit_state$W,
    ns_factor = shocked_state$ns_factor,
    shock_matrix = shock_fit$shock_matrix,
    shocked_factors = shock_fit$shocked_factors,
    shocked_yields = shocked_state$yields,
    shocked_tenor = shocked_fit_state$tenor,
    shocked_curve = shocked_fit_state$curve,
    shocked_single_curve_state = shocked_state,
    baseline_state = baseline_state
  ))
}

mc_fit_shock_refit_from_baseline <- function(baseline_state,
                                             shock_curves, shock_factors, shock_magnitude,
                                             shock_type = c("additive", "multiplicative"),
                                             shock_window = NULL,
                                             shock_window_position = c("first", "last"),
                                             ECM_estim = NULL, ECM_type = NULL, ECM_alpha = NULL,
                                             ECM_fallback = c("error", "baseline"),
                                             X_normalize = NULL, X_trunc = NULL,
                                             CR_algo = NULL, CR_init = NULL,
                                             Xt_smooth = NULL, smoother = NULL,
                                             knot_count = NULL, k_count = NULL, k_pass = NULL,
                                             hp_lambda = NULL, henderson_k = NULL,
                                             CR_maxiter = NULL, CR_tol = NULL,
                                             CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = NULL,
                                             CR_check_every = NULL, CR_use_loglik = NULL, CR_store_path = NULL,
                                             CR_verb = NULL, CR_term = NULL) {
  shock_type <- match.arg(shock_type)
  shock_window_position <- match.arg(shock_window_position)
  ECM_fallback <- match.arg(ECM_fallback)

  mc_fit_shock_from_baseline(
    baseline_state = baseline_state,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window,
    shock_window_position = shock_window_position,
    reuse_BS0 = FALSE,
    reuse_ECM = FALSE,
    anchor_shock_curves = FALSE,
    ECM_estim = ECM_estim,
    ECM_type = ECM_type,
    ECM_alpha = ECM_alpha,
    ECM_fallback = ECM_fallback,
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
    CR_check_every = CR_check_every,
    CR_use_loglik = CR_use_loglik,
    CR_store_path = CR_store_path,
    CR_verb = CR_verb,
    CR_term = CR_term
  )
}

mc_fit_shock <- function(yields = NULL, lambdas = NULL, cutoffs = NULL, reference = NULL,
                         shock_curves, shock_factors, shock_magnitude,
                         shock_type = c("additive", "multiplicative"),
                         shock_window = NULL,
                         shock_window_position = c("first", "last"),
                         reuse_BS0 = TRUE, reuse_ECM = TRUE,
                         anchor_shock_curves = FALSE,
                         ECM_estim = "ML", ECM_type = "eigen", ECM_alpha = 0.1,
                         ECM_fallback = c("error", "baseline"),
                         X_normalize = TRUE, X_trunc = FALSE,
                         CR_algo = zero_mean_covreg_em, CR_init = "adaptive",
                         Xt_smooth = FALSE, smoother = c("ns", "bs", "rm", "hp", "henderson"),
                         knot_count = 5, k_count = 9, k_pass = 3,
                         hp_lambda = 1600, henderson_k = 13,
                         CR_maxiter = 1000, CR_tol = 1e-8,
                         CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = 0,
                         CR_check_every = 1, CR_use_loglik = TRUE, CR_store_path = TRUE,
                         CR_verb = FALSE, CR_term = TRUE,
                         curves = NULL, mats = NULL, baseline_state = NULL) {
  ECM_fallback <- match.arg(ECM_fallback)

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
      CR_check_every = CR_check_every,
      CR_use_loglik = CR_use_loglik,
      CR_store_path = CR_store_path,
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
    shock_window_position = shock_window_position,
    reuse_BS0 = reuse_BS0,
    reuse_ECM = reuse_ECM,
    anchor_shock_curves = anchor_shock_curves,
    ECM_estim = ECM_estim,
    ECM_type = ECM_type,
    ECM_alpha = ECM_alpha,
    ECM_fallback = ECM_fallback,
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
    CR_check_every = CR_check_every,
    CR_use_loglik = CR_use_loglik,
    CR_store_path = CR_store_path,
    CR_verb = CR_verb,
    CR_term = CR_term
  )
}
