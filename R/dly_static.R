#############################################
#        Global Factor Decomposition        #
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

dly_fit_loadings <- function(country_factor, global_factor) {
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

  total_ss <- Matrix::colSums((country_factor - matrix(colMeans(country_factor),
                                               nrow = nrow(country_factor),
                                               ncol = ncol(country_factor),
                                               byrow = TRUE))^2)
  resid_ss <- Matrix::colSums(residual_factor^2)
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

#############################################
#            Yield Reconstruction           #
#############################################

dly_bind_betas <- function(level_factor, slope_factor, curves = colnames(level_factor)) {
  stats::setNames(lapply(curves, function(curve) {
    out <- cbind(
      L = level_factor[, curve],
      S = slope_factor[, curve]
    )
    rownames(out) <- rownames(level_factor)
    out
  }), curves)
}

dly_curve_yields <- function(sc_full, level_factor, slope_factor, curves = names(sc_full)) {
  stats::setNames(lapply(curves, function(curve) {
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
  stats::setNames(lapply(curves, function(curve) {
    as.matrix(actual_yields[[curve]][, -1, drop = FALSE]) -
      as.matrix(fitted_yields[[curve]][, -1, drop = FALSE])
  }), curves)
}

#############################################
#           Country Factor Mapping          #
#############################################

dly_map_global <- function(level_global_factor, slope_global_factor, factor_models,
                           curves) {
  level_intercept <- stats::setNames(as.numeric(factor_models$level$intercept[curves]), curves)
  level_loading <- stats::setNames(as.numeric(factor_models$level$loading[curves]), curves)
  slope_intercept <- stats::setNames(as.numeric(factor_models$slope$intercept[curves]), curves)
  slope_loading <- stats::setNames(as.numeric(factor_models$slope$loading[curves]), curves)

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

#############################################
#           Compatibility Aliases           #
#############################################
dly_country_factor_model <- dly_fit_loadings
dly_country_global_factors <- dly_map_global
