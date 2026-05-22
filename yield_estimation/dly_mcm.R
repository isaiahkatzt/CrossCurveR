#############################################
#     Diebold-Li-Yue Multi-Curve Model      #
#############################################

dly_phi <- function(lambda_optim, mats) {
  phi <- NS_loadings_ls(lambda_optim, mats)
  cross_phi <- solve(tcrossprod(phi)) %*% phi
  return(list(phi = phi, cross_phi = cross_phi))
}

dly_single_curve <- function(lambda_select, yields, mats) {
  curve_phi <- dly_phi(lambda_select, mats)
  curve_betas <- bl_NSbetas(yields, lambda_select, curve_phi$cross_phi, mats)
  curve_nsfit <- bl_NSfit(
    curve_betas$m_yield,
    curve_phi$phi,
    curve_betas$NSbetas,
    lambda_select,
    mats
  )

  return(list(
    lambda = lambda_select,
    time = yields[[1]],
    yield_names = colnames(yields)[-1],
    phi = curve_phi,
    betas = curve_betas,
    nsfit = curve_nsfit
  ))
}

dly_all_sc <- function(lambdas, yields, mats, curves = names(yields)) {
  if (is.null(curves)) {
    curves <- paste0("curve", seq_along(yields))
  }

  if (length(curves) != length(yields)) {
    stop("curves must have the same length as yields.")
  }

  curve_lambdas <- vapply(yields, function(x) {
    lambda_grid_search(lambdas, x, mats, loss = l2_loss_ls)
  }, numeric(1))

  curve_NSparam <- setNames(lapply(seq_along(yields), function(i) {
    dly_single_curve(curve_lambdas[[i]], yields[[i]], mats)
  }), curves)

  return(curve_NSparam)
}

dly_pull_beta <- function(dly_sclist, curves = names(dly_sclist)) {
  if (is.null(curves)) {
    curves <- names(dly_sclist)
  }

  setNames(lapply(curves, function(curve) {
    dly_sclist[[curve]]$betas$NSbetas
  }), curves)
}

dly_group_factor <- function(sc_full, curves = names(sc_full), latent = c("L", "S")) {
  latent <- match.arg(latent, c("L", "S"))
  latent_idx <- ifelse(latent == "L", 1, 2)
  curve_beta <- dly_pull_beta(sc_full, curves)

  factor_draw <- do.call(cbind, lapply(curve_beta, function(beta_mat) {
    beta_mat[, latent_idx]
  }))

  colnames(factor_draw) <- curves
  return(factor_draw)
}

#############################################
#     Simplified Structure: PCA Loadings    #
#############################################

dly_pca_global <- function(factor_draw, reference = colnames(factor_draw)[1],
                           center = TRUE, scale. = FALSE) {
  factor_draw <- as.matrix(factor_draw)

  if (ncol(factor_draw) < 2) {
    stop("factor_draw must contain at least two curves for PCA.")
  }

  pca_fit <- stats::prcomp(factor_draw, center = center, scale. = scale.)
  orientation <- 1

  if (!is.null(reference) && reference %in% rownames(pca_fit$rotation) &&
      pca_fit$rotation[reference, 1] < 0) {
    orientation <- -1
    pca_fit$x[, 1] <- orientation * pca_fit$x[, 1]
    pca_fit$rotation[, 1] <- orientation * pca_fit$rotation[, 1]
  }

  global_factor <- pca_fit$x[, 1]
  global_loading <- pca_fit$rotation[, 1]

  return(list(
    factor = global_factor,
    loading = global_loading,
    pca = pca_fit,
    explained_variance = (pca_fit$sdev[1]^2) / sum(pca_fit$sdev^2)
  ))
}

dly_country_factor_model <- function(country_factor, global_factor) {
  country_factor <- as.matrix(country_factor)
  global_factor <- as.numeric(global_factor)

  if (nrow(country_factor) != length(global_factor)) {
    stop("country_factor and global_factor must have the same number of observations.")
  }

  X <- cbind(Intercept = 1, Global = global_factor)
  coef_mat <- apply(country_factor, 2, function(y) {
    stats::lm.fit(X, y)$coefficients
  })

  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, ncol = 1)
    colnames(coef_mat) <- colnames(country_factor)
  }

  fitted_factor <- X %*% coef_mat
  residual_factor <- country_factor - fitted_factor

  colnames(fitted_factor) <- colnames(country_factor)
  colnames(residual_factor) <- colnames(country_factor)
  rownames(fitted_factor) <- rownames(country_factor)
  rownames(residual_factor) <- rownames(country_factor)

  total_ss <- colSums((country_factor - matrix(colMeans(country_factor),
                                               nrow = nrow(country_factor),
                                               ncol = ncol(country_factor),
                                               byrow = TRUE))^2)
  resid_ss <- colSums(residual_factor^2)
  r_squared <- 1 - resid_ss / total_ss
  r_squared[!is.finite(r_squared)] <- NA_real_

  return(list(
    intercept = coef_mat["Intercept", ],
    loading = coef_mat["Global", ],
    fitted = fitted_factor,
    residual = residual_factor,
    r_squared = r_squared
  ))
}

dly_bind_betas <- function(level_factor, slope_factor, curves = colnames(level_factor)) {
  setNames(lapply(curves, function(curve) {
    out <- cbind(
      L = level_factor[, curve],
      S = slope_factor[, curve]
    )
    rownames(out) <- rownames(level_factor)
    out
  }), curves)
}

dly_curve_yields <- function(sc_full, level_factor, slope_factor, curves = names(sc_full)) {
  setNames(lapply(curves, function(curve) {
    factor_mat <- cbind(
      level_factor[, curve],
      slope_factor[, curve]
    )
    yield_hat <- factor_mat %*% sc_full[[curve]]$phi$phi
    yield_df <- data.frame(sc_full[[curve]]$time, as.data.frame(yield_hat))
    colnames(yield_df) <- c("time", sc_full[[curve]]$yield_names)
    yield_df
  }), curves)
}

dly_curve_residuals <- function(actual_yields, fitted_yields, curves = names(actual_yields)) {
  setNames(lapply(curves, function(curve) {
    as.matrix(actual_yields[[curve]][, -1, drop = FALSE]) -
      as.matrix(fitted_yields[[curve]][, -1, drop = FALSE])
  }), curves)
}

dly_static_fit_from_single_curve <- function(dly_sc, actual_yields = NULL,
                                             curves = names(dly_sc),
                                             reference = curves[[1]],
                                             center = TRUE,
                                             scale. = FALSE,
                                             mats = NULL) {
  if (is.null(curves)) {
    curves <- names(dly_sc)
  }

  if (!reference %in% curves) {
    stop("reference must be one of the curve names supplied in curves.")
  }

  level_factor <- dly_group_factor(dly_sc, curves = curves, latent = "L")
  slope_factor <- dly_group_factor(dly_sc, curves = curves, latent = "S")

  global_level <- dly_pca_global(level_factor, reference = reference, center = center, scale. = scale.)
  global_slope <- dly_pca_global(slope_factor, reference = reference, center = center, scale. = scale.)

  level_model <- dly_country_factor_model(level_factor, global_level$factor)
  slope_model <- dly_country_factor_model(slope_factor, global_slope$factor)

  raw_betas <- dly_bind_betas(level_factor, slope_factor, curves)
  global_betas <- dly_bind_betas(level_model$fitted, slope_model$fitted, curves)
  idio_betas <- dly_bind_betas(level_model$residual, slope_model$residual, curves)
  global_plus_idio_betas <- dly_bind_betas(
    level_model$fitted + level_model$residual,
    slope_model$fitted + slope_model$residual,
    curves
  )

  raw_yields <- dly_curve_yields(dly_sc, level_factor, slope_factor, curves)
  global_yields <- dly_curve_yields(dly_sc, level_model$fitted, slope_model$fitted, curves)
  idio_yields <- dly_curve_yields(dly_sc, level_model$residual, slope_model$residual, curves)
  global_plus_idio_yields <- dly_curve_yields(
    dly_sc,
    level_model$fitted + level_model$residual,
    slope_model$fitted + slope_model$residual,
    curves
  )

  if (is.null(actual_yields)) {
    actual_yields <- raw_yields
  }

  if (is.null(mats)) {
    mats <- unique(unlist(lapply(raw_yields[curves], function(curve_df) {
      numeric_mats <- suppressWarnings(vapply(colnames(curve_df)[-1], function(x) {
        12 * matname_to_numeric(sub("^X", "", x))
      }, numeric(1)))
      numeric_mats[!is.na(numeric_mats)]
    })))
    mats <- sort(as.numeric(mats))
  }

  return(list(
    curves = curves,
    mats = mats,
    reference = reference,
    center = center,
    scale. = scale.,
    lambdas = vapply(dly_sc, `[[`, numeric(1), "lambda"),
    single_curve = dly_sc,
    country_factors = list(
      level = level_factor,
      slope = slope_factor
    ),
    global_factors = list(
      level = global_level,
      slope = global_slope
    ),
    factor_models = list(
      level = level_model,
      slope = slope_model
    ),
    country_betas = list(
      raw = raw_betas,
      global = global_betas,
      idiosyncratic = idio_betas,
      global_plus_idiosyncratic = global_plus_idio_betas
    ),
    yields = list(
      raw = raw_yields,
      global = global_yields,
      idiosyncratic = idio_yields,
      global_plus_idiosyncratic = global_plus_idio_yields
    ),
    residuals = list(
      raw = dly_curve_residuals(actual_yields, raw_yields, curves),
      global = dly_curve_residuals(actual_yields, global_yields, curves),
      global_plus_idiosyncratic = dly_curve_residuals(
        actual_yields,
        global_plus_idio_yields,
        curves
      )
    )
  ))
}

dly_state_matrix <- function(level_series, slope_series, row_template = NULL) {
  state_mat <- cbind(L = as.numeric(level_series), S = as.numeric(slope_series))
  if (!is.null(row_template)) {
    rownames(state_mat) <- row_template
  }
  state_mat
}

dly_fit_var1_block <- function(state_matrix) {
  state_matrix <- as.matrix(state_matrix)

  if (nrow(state_matrix) < 2) {
    stop("state_matrix must contain at least two rows to estimate dynamics.")
  }

  row_template <- rownames(state_matrix)
  lag_state <- state_matrix[-nrow(state_matrix), , drop = FALSE]
  response <- state_matrix[-1, , drop = FALSE]

  X <- cbind(Intercept = 1, lag_state)
  coef_mat <- apply(response, 2, function(y) {
    stats::lm.fit(X, y)$coefficients
  })

  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, ncol = 1)
    colnames(coef_mat) <- colnames(response)
  }

  coef_mat[!is.finite(coef_mat)] <- 0

  intercept <- as.numeric(coef_mat["Intercept", ])
  names(intercept) <- colnames(response)
  phi <- t(coef_mat[colnames(lag_state), , drop = FALSE])
  colnames(phi) <- colnames(lag_state)
  rownames(phi) <- colnames(response)

  fitted_core <- X %*% coef_mat
  fitted_state <- state_matrix
  fitted_state[-1, ] <- fitted_core
  residual_state <- state_matrix - fitted_state
  rownames(fitted_state) <- row_template
  rownames(residual_state) <- row_template

  list(
    intercept = intercept,
    phi = phi,
    fitted = fitted_state,
    residual = residual_state
  )
}

dly_fit_country_dynamic_block <- function(residual_matrix, global_matrix) {
  residual_matrix <- as.matrix(residual_matrix)
  global_matrix <- as.matrix(global_matrix)

  if (nrow(residual_matrix) != nrow(global_matrix)) {
    stop("residual_matrix and global_matrix must align by time.")
  }

  if (nrow(residual_matrix) < 2) {
    stop("residual_matrix must contain at least two rows to estimate dynamics.")
  }

  row_template <- rownames(residual_matrix)
  residual_lag <- residual_matrix[-nrow(residual_matrix), , drop = FALSE]
  global_lag <- global_matrix[-nrow(global_matrix), , drop = FALSE]
  response <- residual_matrix[-1, , drop = FALSE]

  colnames(residual_lag) <- paste0(colnames(residual_lag), "_lag")
  colnames(global_lag) <- paste0("g", colnames(global_lag), "_lag")
  X_no_intercept <- cbind(residual_lag, global_lag)
  X <- cbind(Intercept = 1, X_no_intercept)

  coef_mat <- apply(response, 2, function(y) {
    stats::lm.fit(X, y)$coefficients
  })

  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, ncol = 1)
    colnames(coef_mat) <- colnames(response)
  }

  coef_mat[!is.finite(coef_mat)] <- 0

  intercept <- as.numeric(coef_mat["Intercept", ])
  names(intercept) <- colnames(response)
  phi <- t(coef_mat[paste0(colnames(residual_matrix), "_lag"), , drop = FALSE])
  kappa <- t(coef_mat[paste0("g", colnames(global_matrix), "_lag"), , drop = FALSE])
  colnames(phi) <- colnames(residual_matrix)
  rownames(phi) <- colnames(response)
  colnames(kappa) <- colnames(global_matrix)
  rownames(kappa) <- colnames(response)

  fitted_core <- X %*% coef_mat
  fitted_state <- residual_matrix
  fitted_state[-1, ] <- fitted_core
  residual_state <- residual_matrix - fitted_state
  rownames(fitted_state) <- row_template
  rownames(residual_state) <- row_template

  list(
    intercept = intercept,
    phi = phi,
    kappa = kappa,
    fitted = fitted_state,
    residual = residual_state
  )
}

dly_transition_var1 <- function(prev_state, intercept, phi) {
  as.numeric(intercept + phi %*% as.numeric(prev_state))
}

dly_transition_country_residual <- function(prev_residual, prev_global, intercept, phi, kappa) {
  as.numeric(intercept + phi %*% as.numeric(prev_residual) + kappa %*% as.numeric(prev_global))
}

dly_country_global_factors <- function(level_global_factor, slope_global_factor, factor_models,
                                       curves) {
  level_intercept <- setNames(as.numeric(factor_models$level$intercept[curves]), curves)
  level_loading <- setNames(as.numeric(factor_models$level$loading[curves]), curves)
  slope_intercept <- setNames(as.numeric(factor_models$slope$intercept[curves]), curves)
  slope_loading <- setNames(as.numeric(factor_models$slope$loading[curves]), curves)

  level_factor <- outer(as.numeric(level_global_factor), level_loading)
  slope_factor <- outer(as.numeric(slope_global_factor), slope_loading)
  level_factor <- sweep(level_factor, 2, level_intercept, FUN = "+")
  slope_factor <- sweep(slope_factor, 2, slope_intercept, FUN = "+")
  colnames(level_factor) <- curves
  colnames(slope_factor) <- curves

  row_template <- NULL
  if (!is.null(rownames(factor_models$level$fitted))) {
    row_template <- rownames(factor_models$level$fitted)
  } else if (!is.null(rownames(factor_models$slope$fitted))) {
    row_template <- rownames(factor_models$slope$fitted)
  }

  if (!is.null(row_template)) {
    rownames(level_factor) <- row_template
    rownames(slope_factor) <- row_template
  }

  list(level = level_factor, slope = slope_factor)
}

dly_propagate_residual_states <- function(residual_models, actual_residuals, global_path,
                                          curves) {
  setNames(lapply(curves, function(curve) {
    model <- residual_models[[curve]]
    actual_state <- as.matrix(actual_residuals[[curve]])
    propagated_state <- actual_state

    for (t in seq_len(nrow(actual_state))[-1]) {
      propagated_state[t, ] <- dly_transition_country_residual(
        prev_residual = propagated_state[t - 1, ],
        prev_global = global_path[t - 1, ],
        intercept = model$intercept,
        phi = model$phi,
        kappa = model$kappa
      )
    }

    rownames(propagated_state) <- rownames(actual_state)
    colnames(propagated_state) <- colnames(actual_state)
    propagated_state
  }), curves)
}

dly_unpack_country_state <- function(country_state_list, curves) {
  level_factor <- do.call(cbind, lapply(curves, function(curve) country_state_list[[curve]][, "L"]))
  slope_factor <- do.call(cbind, lapply(curves, function(curve) country_state_list[[curve]][, "S"]))
  colnames(level_factor) <- curves
  colnames(slope_factor) <- curves

  row_template <- rownames(country_state_list[[curves[[1]]]])
  rownames(level_factor) <- row_template
  rownames(slope_factor) <- row_template

  list(level = level_factor, slope = slope_factor)
}

dly_dynamic_fit_from_static <- function(static_fit, actual_yields = NULL) {
  curves <- static_fit$curves
  row_template <- rownames(static_fit$country_factors$level)

  actual_global <- dly_state_matrix(
    static_fit$global_factors$level$factor,
    static_fit$global_factors$slope$factor,
    row_template = row_template
  )
  actual_residual <- setNames(lapply(curves, function(curve) {
    dly_state_matrix(
      static_fit$factor_models$level$residual[, curve],
      static_fit$factor_models$slope$residual[, curve],
      row_template = row_template
    )
  }), curves)

  global_transition <- dly_fit_var1_block(actual_global)
  residual_transition <- setNames(lapply(curves, function(curve) {
    dly_fit_country_dynamic_block(actual_residual[[curve]], actual_global)
  }), curves)

  propagated_global <- actual_global
  for (t in seq_len(nrow(actual_global))[-1]) {
    propagated_global[t, ] <- dly_transition_var1(
      prev_state = propagated_global[t - 1, ],
      intercept = global_transition$intercept,
      phi = global_transition$phi
    )
  }
  rownames(propagated_global) <- row_template
  colnames(propagated_global) <- colnames(actual_global)

  propagated_residual <- dly_propagate_residual_states(
    residual_models = residual_transition,
    actual_residuals = actual_residual,
    global_path = propagated_global,
    curves = curves
  )

  propagated_global_factors <- dly_country_global_factors(
    level_global_factor = propagated_global[, "L"],
    slope_global_factor = propagated_global[, "S"],
    factor_models = static_fit$factor_models,
    curves = curves
  )
  propagated_residual_factors <- dly_unpack_country_state(propagated_residual, curves)
  propagated_betas <- dly_bind_betas(
    propagated_global_factors$level + propagated_residual_factors$level,
    propagated_global_factors$slope + propagated_residual_factors$slope,
    curves
  )
  propagated_yields <- dly_curve_yields(
    static_fit$single_curve,
    propagated_global_factors$level + propagated_residual_factors$level,
    propagated_global_factors$slope + propagated_residual_factors$slope,
    curves
  )

  if (is.null(actual_yields)) {
    actual_yields <- static_fit$yields$raw
  }

  dynamic_fit <- static_fit
  dynamic_fit$country_betas$dislocation_propagated <- propagated_betas
  dynamic_fit$yields$dislocation_propagated <- propagated_yields
  dynamic_fit$residuals$dislocation_propagated <- dly_curve_residuals(
    actual_yields,
    propagated_yields,
    curves
  )
  dynamic_fit$dynamics <- list(
    global = c(
      global_transition,
      list(actual = actual_global, propagated = propagated_global)
    ),
    residual = setNames(lapply(curves, function(curve) {
      c(
        residual_transition[[curve]],
        list(actual = actual_residual[[curve]], propagated = propagated_residual[[curve]])
      )
    }), curves)
  )

  dynamic_fit
}

dly_fit_from_single_curve <- function(dly_sc, actual_yields = NULL,
                                      curves = names(dly_sc),
                                      reference = curves[[1]],
                                      center = TRUE,
                                      scale. = FALSE,
                                      mats = NULL) {
  static_fit <- dly_static_fit_from_single_curve(
    dly_sc = dly_sc,
    actual_yields = actual_yields,
    curves = curves,
    reference = reference,
    center = center,
    scale. = scale.,
    mats = mats
  )

  dly_dynamic_fit_from_static(static_fit, actual_yields = actual_yields)
}

dly_fit <- function(yields, mats,
                    lambdas = seq(from = 0.001, to = 1, length.out = 100),
                    curves = names(yields),
                    reference = curves[[1]],
                    center = TRUE,
                    scale. = FALSE) {
  if (is.null(curves)) {
    curves <- paste0("curve", seq_along(yields))
  }

  dly_sc <- dly_all_sc(lambdas = lambdas, yields = yields, mats = mats, curves = curves)

  dly_fit_obj <- dly_fit_from_single_curve(
    dly_sc = dly_sc,
    actual_yields = yields,
    curves = curves,
    reference = reference,
    center = center,
    scale. = scale.,
    mats = mats
  )

  dly_fit_obj$lambda_grid <- lambdas
  dly_fit_obj
}
