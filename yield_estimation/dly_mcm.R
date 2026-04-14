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
  global_factor <- pca_fit$x[, 1]
  global_loading <- pca_fit$rotation[, 1]

  if (!is.null(reference) && reference %in% names(global_loading) &&
      global_loading[[reference]] < 0) {
    global_factor <- -global_factor
    global_loading <- -global_loading
  }

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

dly_fit <- function(yields, mats,
                    lambdas = seq(from = 0.001, to = 1, length.out = 100),
                    curves = names(yields),
                    reference = curves[[1]],
                    center = TRUE,
                    scale. = FALSE) {
  if (is.null(curves)) {
    curves <- paste0("curve", seq_along(yields))
  }

  if (!reference %in% curves) {
    stop("reference must be one of the curve names supplied in curves.")
  }

  dly_sc <- dly_all_sc(lambdas = lambdas, yields = yields, mats = mats, curves = curves)

  level_factor <- dly_group_factor(dly_sc, curves = curves, latent = "L")
  slope_factor <- dly_group_factor(dly_sc, curves = curves, latent = "S")

  global_level <- dly_pca_global(level_factor, reference = reference, center = center, scale. = scale.)
  global_slope <- dly_pca_global(slope_factor, reference = reference, center = center, scale. = scale.)

  level_model <- dly_country_factor_model(level_factor, global_level$factor)
  slope_model <- dly_country_factor_model(slope_factor, global_slope$factor)

  raw_betas <- dly_bind_betas(level_factor, slope_factor, curves)
  global_betas <- dly_bind_betas(level_model$fitted, slope_model$fitted, curves)
  idio_betas <- dly_bind_betas(level_model$residual, slope_model$residual, curves)

  raw_yields <- dly_curve_yields(dly_sc, level_factor, slope_factor, curves)
  global_yields <- dly_curve_yields(dly_sc, level_model$fitted, slope_model$fitted, curves)
  idio_yields <- dly_curve_yields(dly_sc, level_model$residual, slope_model$residual, curves)

  return(list(
    curves = curves,
    mats = mats,
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
      idiosyncratic = idio_betas
    ),
    yields = list(
      raw = raw_yields,
      global = global_yields,
      idiosyncratic = idio_yields
    ),
    residuals = list(
      raw = dly_curve_residuals(yields, raw_yields, curves),
      global = dly_curve_residuals(yields, global_yields, curves)
    )
  ))
}
