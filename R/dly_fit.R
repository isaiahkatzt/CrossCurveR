#############################################
#         Static Multi-Curve Fitting        #
#############################################

dly_fit_static_sc <- function(dly_sc, actual_yields = NULL,
                              curves = names(dly_sc),
                              reference = curves[[1]],
                              center = TRUE,
                              scale. = FALSE,
                              mats = NULL) {
  if (is.null(curves)) {
    curves <- names(dly_sc)
  }
  validated <- .cc_dly_sc(dly_sc, mats, curves)
  dly_sc <- validated$sc
  mats <- validated$mats

  .cc_reference(reference, curves)

  level_factor <- dly_factor_panel(dly_sc, curves = curves, latent = "L")
  slope_factor <- dly_factor_panel(dly_sc, curves = curves, latent = "S")

  global_level <- dly_pca_global(level_factor, reference = reference, center = center, scale. = scale.)
  global_slope <- dly_pca_global(slope_factor, reference = reference, center = center, scale. = scale.)

  level_model <- dly_fit_loadings(level_factor, global_level$factor)
  slope_model <- dly_fit_loadings(slope_factor, global_slope$factor)

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
  actual_yields <- .cc_yields(actual_yields, mats, curves)
  if (!identical(.cc_index(actual_yields[[1]][[1]]), .cc_index(dly_sc[[1]]$time))) {
    stop("Observed yields must match the country-fit observation index.", call. = FALSE)
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
    input_spec = .cc_input_spec(actual_yields, mats, curves,
      decay = list(selected = vapply(dly_sc, `[[`, numeric(1), "lambda"),
                   parameterization = "inverse months")),
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

#############################################
#        Dynamic Multi-Curve Fitting        #
#############################################

dly_fit_dynamic <- function(static_fit, actual_yields = NULL) {
  curves <- static_fit$curves
  row_template <- rownames(static_fit$country_factors$level)

  actual_global <- dly_state_matrix(
    static_fit$global_factors$level$factor,
    static_fit$global_factors$slope$factor,
    row_template = row_template
  )
  actual_residual <- stats::setNames(lapply(curves, function(curve) {
    dly_state_matrix(
      static_fit$factor_models$level$residual[, curve],
      static_fit$factor_models$slope$residual[, curve],
      row_template = row_template
    )
  }), curves)

  global_transition <- dly_fit_var1(actual_global)
  residual_transition <- stats::setNames(lapply(curves, function(curve) {
    dly_fit_rs(actual_residual[[curve]], actual_global)
  }), curves)

  propagated_global <- actual_global
  for (t in seq_len(nrow(actual_global))[-1]) {
    propagated_global[t, ] <- dly_step_var1(
      prev_state = propagated_global[t - 1, ],
      intercept = global_transition$intercept,
      phi = global_transition$phi
    )
  }
  rownames(propagated_global) <- row_template
  colnames(propagated_global) <- colnames(actual_global)

  propagated_residual <- dly_propagate_rs(
    residual_models = residual_transition,
    actual_residuals = actual_residual,
    global_path = propagated_global,
    curves = curves
  )

  propagated_global_factors <- dly_map_global(
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
  actual_yields <- .cc_yields(actual_yields, static_fit$mats, curves)
  if (!identical(.cc_index(actual_yields[[1]][[1]]),
                 .cc_index(static_fit$single_curve[[1]]$time))) {
    stop("Observed yields must match the static-fit observation index.", call. = FALSE)
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
    residual = stats::setNames(lapply(curves, function(curve) {
      c(
        residual_transition[[curve]],
        list(actual = actual_residual[[curve]], propagated = propagated_residual[[curve]])
      )
    }), curves)
  )

  dynamic_fit
}

#############################################
#      Multi-Curve Fitting Entrypoints      #
#############################################

dly_fit_from_sc <- function(dly_sc, actual_yields = NULL,
                            curves = names(dly_sc),
                            reference = curves[[1]],
                            center = TRUE,
                            scale. = FALSE,
                            mats = NULL) {
  static_fit <- dly_fit_static_sc(
    dly_sc = dly_sc,
    actual_yields = actual_yields,
    curves = curves,
    reference = reference,
    center = center,
    scale. = scale.,
    mats = mats
  )

  dly_fit_dynamic(static_fit, actual_yields = actual_yields)
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
  .cc_reference(reference, curves)
  yields <- .cc_yields(yields, mats, curves, reference)

  dly_sc <- dly_sc_fit_all(lambdas = lambdas, yields = yields, mats = mats, curves = curves)

  dly_fit_obj <- dly_fit_from_sc(
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

#############################################
#           Compatibility Aliases           #
#############################################
dly_static_fit_from_single_curve <- dly_fit_static_sc
dly_dynamic_fit_from_static <- dly_fit_dynamic
dly_fit_from_single_curve <- dly_fit_from_sc
