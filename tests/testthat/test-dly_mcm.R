make_ls_curve_df <- function(level_factor, slope_factor, lambda = 0.35,
                             mats = fixture_mats, time = fixture_time) {
  yield_mat <- cbind(level_factor, slope_factor) %*% NS_loadings_ls(lambda, mats)
  curve_df <- data.frame(time = time, as.data.frame(yield_mat))
  colnames(curve_df) <- c("time", sapply(mats, numeric_to_matname))
  curve_df
}

make_dly_fixture <- function() {
  idx <- seq_along(fixture_time) - 1
  global_level <- numeric(length(idx))
  global_slope <- numeric(length(idx))
  global_level[[1]] <- 2
  global_slope[[1]] <- -1

  for (t in 2:length(idx)) {
    global_level[[t]] <- 0.15 + 0.92 * global_level[[t - 1]]
    global_slope[[t]] <- -0.05 + 0.65 * global_slope[[t - 1]]
  }

  usa_level_resid <- numeric(length(idx))
  usa_slope_resid <- numeric(length(idx))
  gbr_level_resid <- numeric(length(idx))
  gbr_slope_resid <- numeric(length(idx))
  jpn_level_resid <- numeric(length(idx))
  jpn_slope_resid <- numeric(length(idx))

  for (t in 2:length(idx)) {
    usa_level_resid[[t]] <- 0.70 * usa_level_resid[[t - 1]] + 0.08 * global_level[[t - 1]]
    usa_slope_resid[[t]] <- 0.50 * usa_slope_resid[[t - 1]] + 0.05 * global_slope[[t - 1]]
    gbr_level_resid[[t]] <- 0.55 * gbr_level_resid[[t - 1]] - 0.04 * global_level[[t - 1]]
    gbr_slope_resid[[t]] <- 0.45 * gbr_slope_resid[[t - 1]] + 0.03 * global_slope[[t - 1]]
    jpn_level_resid[[t]] <- 0.35 * jpn_level_resid[[t - 1]] + 0.02 * global_level[[t - 1]]
    jpn_slope_resid[[t]] <- 0.40 * jpn_slope_resid[[t - 1]] - 0.02 * global_slope[[t - 1]]
  }

  list(
    global_level = global_level,
    global_slope = global_slope,
    yields = list(
      usa = make_ls_curve_df(
        level_factor = 0.4 + 1.00 * global_level + usa_level_resid,
        slope_factor = -0.2 + 0.85 * global_slope + usa_slope_resid
      ),
      gbr = make_ls_curve_df(
        level_factor = 0.2 + 1.15 * global_level + gbr_level_resid,
        slope_factor = 0.1 + 0.70 * global_slope + gbr_slope_resid
      ),
      jpn = make_ls_curve_df(
        level_factor = -0.3 + 0.90 * global_level + jpn_level_resid,
        slope_factor = -0.1 + 0.45 * global_slope + jpn_slope_resid
      )
    )
  )
}

test_that("dly_all_sc selects the exact lambda for two-factor NS data", {
  fixture <- make_dly_fixture()
  lambda_grid <- c(0.2, 0.35, 0.5)

  dly_sc <- dly_all_sc(lambda_grid, fixture$yields, fixture_mats, curves = names(fixture$yields))

  expect_equal(unname(vapply(dly_sc, `[[`, numeric(1), "lambda")), rep(0.35, 3))
  expect_equal(names(dly_sc), names(fixture$yields))
})

test_that("dly_pca_global returns a reference-aligned first principal component", {
  fixture <- make_dly_fixture()
  dly_sc <- dly_all_sc(c(0.35), fixture$yields, fixture_mats, curves = names(fixture$yields))
  level_factor <- dly_group_factor(dly_sc, curves = names(fixture$yields), latent = "L")

  global_level <- dly_pca_global(level_factor, reference = "usa")

  expect_gt(global_level$loading[["usa"]], 0)
  expect_equal(
    unname(global_level$loading),
    unname(global_level$pca$rotation[, 1]),
    tolerance = 1e-8
  )
  expect_gt(global_level$explained_variance, 0.95)
  expect_equal(length(global_level$factor), nrow(level_factor))
})

test_that("dly_fit preserves raw fits, static decomposition identities, and dynamic state blocks", {
  fixture <- make_dly_fixture()
  dly_model <- dly_fit(
    yields = fixture$yields,
    mats = fixture_mats,
    lambdas = c(0.2, 0.35, 0.5),
    curves = names(fixture$yields),
    reference = "usa"
  )

  expect_true("dislocation_propagated" %in% names(dly_model$yields))
  expect_true("dislocation_propagated" %in% names(dly_model$country_betas))
  expect_true("dynamics" %in% names(dly_model))
  expect_equal(colnames(dly_model$dynamics$global$actual), c("L", "S"))
  expect_equal(colnames(dly_model$dynamics$global$propagated), c("L", "S"))

  for (curve in names(fixture$yields)) {
    actual_mat <- as.matrix(fixture$yields[[curve]][, -1, drop = FALSE])
    raw_fit_mat <- as.matrix(dly_model$yields$raw[[curve]][, -1, drop = FALSE])
    global_fit_mat <- as.matrix(dly_model$yields$global[[curve]][, -1, drop = FALSE])
    idio_fit_mat <- as.matrix(dly_model$yields$idiosyncratic[[curve]][, -1, drop = FALSE])
    propagated_fit_mat <- as.matrix(dly_model$yields$dislocation_propagated[[curve]][, -1, drop = FALSE])

    expect_equal(raw_fit_mat, actual_mat, tolerance = 1e-6)
    expect_equal(raw_fit_mat, global_fit_mat + idio_fit_mat, tolerance = 1e-6)
    expect_equal(dim(propagated_fit_mat), dim(actual_mat))

    raw_beta <- dly_model$country_betas$raw[[curve]]
    global_beta <- dly_model$country_betas$global[[curve]]
    idio_beta <- dly_model$country_betas$idiosyncratic[[curve]]
    propagated_beta <- dly_model$country_betas$dislocation_propagated[[curve]]

    expect_equal(raw_beta, global_beta + idio_beta, tolerance = 1e-6)
    expect_equal(dim(propagated_beta), dim(raw_beta))
    expect_equal(colnames(dly_model$dynamics$residual[[curve]]$propagated), c("L", "S"))
  }

  raw_resid_norm <- max(abs(unlist(dly_model$residuals$raw)))
  expect_lt(raw_resid_norm, 1e-6)
  expect_gt(
    max(abs(
      as.matrix(dly_model$yields$dislocation_propagated$usa[, -1, drop = FALSE]) -
        as.matrix(dly_model$yields$raw$usa[, -1, drop = FALSE])
    )),
    0
  )
})

test_that("dly_fit_from_single_curve rebuilds the same dynamic DLY outputs from baseline single-curve fits", {
  fixture <- make_dly_fixture()
  dly_model <- dly_fit(
    yields = fixture$yields,
    mats = fixture_mats,
    lambdas = c(0.2, 0.35, 0.5),
    curves = names(fixture$yields),
    reference = "usa"
  )

  rebuilt_model <- dly_fit_from_single_curve(
    dly_sc = dly_model$single_curve,
    actual_yields = fixture$yields,
    curves = dly_model$curves,
    reference = dly_model$reference,
    center = dly_model$center,
    scale. = dly_model$scale.,
    mats = dly_model$mats
  )

  expect_equal(rebuilt_model$curves, dly_model$curves)
  expect_equal(rebuilt_model$mats, dly_model$mats)
  expect_equal(rebuilt_model$dynamics$global$intercept, dly_model$dynamics$global$intercept, tolerance = 1e-8)
  expect_equal(rebuilt_model$dynamics$global$phi, dly_model$dynamics$global$phi, tolerance = 1e-8)

  for (curve in names(fixture$yields)) {
    expect_equal(rebuilt_model$country_betas$raw[[curve]], dly_model$country_betas$raw[[curve]], tolerance = 1e-8)
    expect_equal(
      as.matrix(rebuilt_model$yields$raw[[curve]][, -1, drop = FALSE]),
      as.matrix(dly_model$yields$raw[[curve]][, -1, drop = FALSE]),
      tolerance = 1e-8
    )
    expect_equal(
      as.matrix(rebuilt_model$yields$global[[curve]][, -1, drop = FALSE]),
      as.matrix(dly_model$yields$global[[curve]][, -1, drop = FALSE]),
      tolerance = 1e-8
    )
    expect_equal(
      as.matrix(rebuilt_model$yields$dislocation_propagated[[curve]][, -1, drop = FALSE]),
      as.matrix(dly_model$yields$dislocation_propagated[[curve]][, -1, drop = FALSE]),
      tolerance = 1e-8
    )
  }
})
