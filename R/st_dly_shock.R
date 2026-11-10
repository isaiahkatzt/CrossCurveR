#############################################
#          DLY Shock Reconstruction         #
#############################################

st_dly_sc <- function(dly_single_curve, shocked_factors,
                      curves = names(dly_single_curve)) {
  factor_names <- st_factor_names(ncol(as.matrix(shocked_factors[[1]])))

  stats::setNames(lapply(curves, function(curve_name) {
    baseline_curve <- dly_single_curve[[curve_name]]
    shocked_beta <- .st_factor_matrix(shocked_factors[[curve_name]], factor_names)
    colnames(shocked_beta) <- factor_names
    shocked_nsfit <- shocked_beta %*% baseline_curve$phi$phi

    baseline_curve$betas$NSbetas <- shocked_beta
    baseline_curve$betas$m_yield <- shocked_nsfit
    baseline_curve$nsfit$NSfit <- shocked_nsfit
    baseline_curve
  }), curves)
}

st_dly_project <- function(factor_draw, global_factor_obj) {
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

st_dly_predict <- function(factor_model, global_factor, curves) {
  intercept <- stats::setNames(as.numeric(factor_model$intercept[curves]), curves)
  loading <- stats::setNames(as.numeric(factor_model$loading[curves]), curves)

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

.st_dly_parts <- function(fit, single_curve, global_level, global_slope,
                          level_residual, slope_residual) {
  curves <- fit$curves
  level <- st_dly_predict(fit$factor_models$level, global_level, curves)
  slope <- st_dly_predict(fit$factor_models$slope, global_slope, curves)
  combined_level <- level + level_residual
  combined_slope <- slope + slope_residual

  list(
    level = level,
    slope = slope,
    betas = list(
      global = dly_bind_betas(level, slope, curves),
      idiosyncratic = dly_bind_betas(level_residual, slope_residual, curves),
      global_plus_idiosyncratic = dly_bind_betas(combined_level, combined_slope, curves)
    ),
    yields = list(
      global = dly_curve_yields(single_curve, level, slope, curves),
      global_plus_idiosyncratic = dly_curve_yields(single_curve, combined_level, combined_slope, curves)
    )
  )
}

st_dly_global_shock <- function(dly_fit_obj, shock_factors,
                                shock_profile,
                                shock_type = c("additive", "multiplicative")) {
  shock_type <- match.arg(shock_type)
  if (!identical(shock_type, "additive")) {
    stop("DLY global shocks currently support additive shocks only.")
  }
  if (!all(shock_factors %in% c("L", "S"))) {
    stop("DLY global shocks support `L` and/or `S` shock factors.")
  }

  shock_magnitudes <- if (is.data.frame(shock_profile)) {
    shock_profile$shock_magnitude
  } else {
    shock_profile
  }
  shock_magnitudes <- as.numeric(shock_magnitudes)

  if (any(!is.finite(shock_magnitudes))) {
    stop("`shock_profile` magnitudes must contain only finite values.", call. = FALSE)
  }

  baseline_global_level <- as.numeric(dly_fit_obj$global_factors$level$factor)
  baseline_global_slope <- as.numeric(dly_fit_obj$global_factors$slope$factor)
  if (length(shock_magnitudes) != length(baseline_global_level)) {
    stop("Shock profile length must match the DLY global factor length.")
  }

  shocked_global_level <- baseline_global_level
  shocked_global_slope <- baseline_global_slope

  if ("L" %in% shock_factors) {
    loading_mean <- mean(as.numeric(dly_fit_obj$factor_models$level$loading[dly_fit_obj$curves]))
    if (!is.finite(loading_mean) || abs(loading_mean) <= sqrt(.Machine$double.eps)) {
      stop("Average DLY level loading is too close to zero for global shock scaling.")
    }
    shocked_global_level <- baseline_global_level + shock_magnitudes / loading_mean
  }

  if ("S" %in% shock_factors) {
    loading_mean <- mean(as.numeric(dly_fit_obj$factor_models$slope$loading[dly_fit_obj$curves]))
    if (!is.finite(loading_mean) || abs(loading_mean) <= sqrt(.Machine$double.eps)) {
      stop("Average DLY slope loading is too close to zero for global shock scaling.")
    }
    shocked_global_slope <- baseline_global_slope + shock_magnitudes / loading_mean
  }


  baseline_level_residual <- as.matrix(
    dly_fit_obj$factor_models$level$residual[, dly_fit_obj$curves, drop = FALSE]
  )
  baseline_slope_residual <- as.matrix(
    dly_fit_obj$factor_models$slope$residual[, dly_fit_obj$curves, drop = FALSE]
  )

  parts <- .st_dly_parts(
    dly_fit_obj, dly_fit_obj$single_curve, shocked_global_level, shocked_global_slope,
    baseline_level_residual, baseline_slope_residual
  )
  raw_betas <- dly_fit_obj$country_betas$raw
  raw_yields <- dly_fit_obj$yields$raw
  idiosyncratic_yields <- dly_fit_obj$yields$idiosyncratic

  shocked_fit <- dly_fit_obj
  shocked_fit$global_factors$level$factor <- shocked_global_level
  shocked_fit$global_factors$slope$factor <- shocked_global_slope
  shocked_fit$factor_models$level$fitted <- parts$level
  shocked_fit$factor_models$slope$fitted <- parts$slope
  shocked_fit$factor_models$level$residual <- baseline_level_residual
  shocked_fit$factor_models$slope$residual <- baseline_slope_residual
  shocked_fit$country_betas$raw <- raw_betas
  shocked_fit$country_betas$global <- parts$betas$global
  shocked_fit$country_betas$idiosyncratic <- parts$betas$idiosyncratic
  shocked_fit$country_betas$global_plus_idiosyncratic <- parts$betas$global_plus_idiosyncratic
  shocked_fit$yields$raw <- raw_yields
  shocked_fit$yields$global <- parts$yields$global
  shocked_fit$yields$idiosyncratic <- idiosyncratic_yields
  shocked_fit$yields$global_plus_idiosyncratic <- parts$yields$global_plus_idiosyncratic
  shocked_fit$residuals$global <- dly_curve_residuals(raw_yields, parts$yields$global, dly_fit_obj$curves)
  shocked_fit$residuals$global_plus_idiosyncratic <- dly_curve_residuals(
    raw_yields,
    parts$yields$global_plus_idiosyncratic,
    dly_fit_obj$curves
  )

  list(
    fit = shocked_fit,
    shock_profile = shock_profile,
    baseline_fit = dly_fit_obj
  )
}

st_dly_shock <- function(dly_fit_obj, shock_curves, shock_factors,
                         shock_profile,
                         shock_type = c("additive", "multiplicative")) {
  shock_type <- match.arg(shock_type)
  shock_magnitudes <- if (is.data.frame(shock_profile)) {
    shock_profile$shock_magnitude
  } else {
    shock_profile
  }

  shock_fit <- st_factor_profile(
    curves = dly_fit_obj$curves,
    curve_dns_factor = dly_pull_beta(dly_fit_obj$single_curve, curves = dly_fit_obj$curves),
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitudes = shock_magnitudes,
    shock_type = shock_type
  )

  shocked_single_curve <- st_dly_sc(
    dly_single_curve = dly_fit_obj$single_curve,
    shocked_factors = shock_fit$shocked_factors,
    curves = dly_fit_obj$curves
  )

  raw_level_factor <- dly_group_factor(shocked_single_curve, curves = dly_fit_obj$curves, latent = "L")
  raw_slope_factor <- dly_group_factor(shocked_single_curve, curves = dly_fit_obj$curves, latent = "S")

  shocked_global_level <- st_dly_project(
    factor_draw = raw_level_factor,
    global_factor_obj = dly_fit_obj$global_factors$level
  )
  shocked_global_slope <- st_dly_project(
    factor_draw = raw_slope_factor,
    global_factor_obj = dly_fit_obj$global_factors$slope
  )


  baseline_level_residual <- as.matrix(
    dly_fit_obj$factor_models$level$residual[, dly_fit_obj$curves, drop = FALSE]
  )
  baseline_slope_residual <- as.matrix(
    dly_fit_obj$factor_models$slope$residual[, dly_fit_obj$curves, drop = FALSE]
  )

  shocked_level_residual <- baseline_level_residual
  shocked_slope_residual <- baseline_slope_residual
  directly_shocked_curves <- intersect(shock_curves, dly_fit_obj$curves)
  if (length(directly_shocked_curves) > 0) {
    baseline_level_factor <- as.matrix(
      dly_fit_obj$country_factors$level[, dly_fit_obj$curves, drop = FALSE]
    )
    baseline_slope_factor <- as.matrix(
      dly_fit_obj$country_factors$slope[, dly_fit_obj$curves, drop = FALSE]
    )
    shocked_level_residual[, directly_shocked_curves] <-
      baseline_level_residual[, directly_shocked_curves, drop = FALSE] +
      raw_level_factor[, directly_shocked_curves, drop = FALSE] -
      baseline_level_factor[, directly_shocked_curves, drop = FALSE]
    shocked_slope_residual[, directly_shocked_curves] <-
      baseline_slope_residual[, directly_shocked_curves, drop = FALSE] +
      raw_slope_factor[, directly_shocked_curves, drop = FALSE] -
      baseline_slope_factor[, directly_shocked_curves, drop = FALSE]
  }

  parts <- .st_dly_parts(
    dly_fit_obj, shocked_single_curve, shocked_global_level, shocked_global_slope,
    shocked_level_residual, shocked_slope_residual
  )
  raw_betas <- dly_bind_betas(raw_level_factor, raw_slope_factor, curves = dly_fit_obj$curves)
  raw_yields <- dly_curve_yields(
    shocked_single_curve, raw_level_factor, raw_slope_factor, curves = dly_fit_obj$curves
  )
  idiosyncratic_yields <- dly_curve_yields(
    shocked_single_curve, shocked_level_residual, shocked_slope_residual,
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
        fitted = parts$level,
        residual = shocked_level_residual,
        r_squared = dly_fit_obj$factor_models$level$r_squared
      ),
      slope = list(
        intercept = dly_fit_obj$factor_models$slope$intercept,
        loading = dly_fit_obj$factor_models$slope$loading,
        fitted = parts$slope,
        residual = shocked_slope_residual,
        r_squared = dly_fit_obj$factor_models$slope$r_squared
      )
    ),
    country_betas = list(
      raw = raw_betas,
      global = parts$betas$global,
      idiosyncratic = parts$betas$idiosyncratic,
      global_plus_idiosyncratic = parts$betas$global_plus_idiosyncratic
    ),
    yields = list(
      raw = raw_yields,
      global = parts$yields$global,
      idiosyncratic = idiosyncratic_yields,
      global_plus_idiosyncratic = parts$yields$global_plus_idiosyncratic
    ),
    residuals = list(
      raw = dly_curve_residuals(raw_yields, raw_yields, dly_fit_obj$curves),
      global = dly_curve_residuals(raw_yields, parts$yields$global, dly_fit_obj$curves),
      global_plus_idiosyncratic = dly_curve_residuals(
        raw_yields,
        parts$yields$global_plus_idiosyncratic,
        dly_fit_obj$curves
      )
    )
  )

  list(
    fit = shocked_fit,
    shock_profile = shock_profile,
    shocked_factors = shock_fit$shocked_factors,
    baseline_fit = dly_fit_obj
  )
}

#############################################
#           Compatibility Aliases           #
#############################################

st_build_shocked_dly_single_curve <- st_dly_sc
st_dly_project_global_factor <- st_dly_project
st_dly_predict_country_factor <- st_dly_predict
st_dly_global_shock_from_fit_profile <- st_dly_global_shock
st_dly_shock_from_fit_profile <- st_dly_shock
