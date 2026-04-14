make_ls_curve_df <- function(level_factor, slope_factor, lambda = 0.35,
                             mats = fixture_mats, time = fixture_time) {
  yield_mat <- cbind(level_factor, slope_factor) %*% NS_loadings_ls(lambda, mats)
  curve_df <- data.frame(time = time, as.data.frame(yield_mat))
  colnames(curve_df) <- c("time", sapply(mats, numeric_to_matname))
  curve_df
}

make_dly_fixture <- function() {
  idx <- seq_along(fixture_time) - 1
  global_level <- 2 + 0.08 * idx
  global_slope <- -1 + sin(idx / 3)

  list(
    global_level = global_level,
    global_slope = global_slope,
    yields = list(
      usa = make_ls_curve_df(
        level_factor = 0.4 + 1.00 * global_level + 0.05 * cos(idx / 2),
        slope_factor = -0.2 + 0.85 * global_slope + 0.04 * sin(idx / 4)
      ),
      gbr = make_ls_curve_df(
        level_factor = 0.2 + 1.15 * global_level - 0.04 * sin(idx / 5),
        slope_factor = 0.1 + 0.70 * global_slope + 0.03 * cos(idx / 4)
      ),
      jpn = make_ls_curve_df(
        level_factor = -0.3 + 0.90 * global_level + 0.03 * cos(idx / 6),
        slope_factor = -0.1 + 0.45 * global_slope - 0.02 * sin(idx / 5)
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
  expect_gt(global_level$explained_variance, 0.95)
  expect_equal(length(global_level$factor), nrow(level_factor))
})

test_that("dly_fit preserves exact two-factor fits and additive decomposition identities", {
  fixture <- make_dly_fixture()
  dly_model <- dly_fit(
    yields = fixture$yields,
    mats = fixture_mats,
    lambdas = c(0.2, 0.35, 0.5),
    curves = names(fixture$yields),
    reference = "usa"
  )

  for (curve in names(fixture$yields)) {
    actual_mat <- as.matrix(fixture$yields[[curve]][, -1, drop = FALSE])
    raw_fit_mat <- as.matrix(dly_model$yields$raw[[curve]][, -1, drop = FALSE])
    global_fit_mat <- as.matrix(dly_model$yields$global[[curve]][, -1, drop = FALSE])
    idio_fit_mat <- as.matrix(dly_model$yields$idiosyncratic[[curve]][, -1, drop = FALSE])

    expect_equal(raw_fit_mat, actual_mat, tolerance = 1e-6)
    expect_equal(raw_fit_mat, global_fit_mat + idio_fit_mat, tolerance = 1e-6)

    raw_beta <- dly_model$country_betas$raw[[curve]]
    global_beta <- dly_model$country_betas$global[[curve]]
    idio_beta <- dly_model$country_betas$idiosyncratic[[curve]]

    expect_equal(raw_beta, global_beta + idio_beta, tolerance = 1e-6)
  }

  raw_resid_norm <- max(abs(unlist(dly_model$residuals$raw)))
  expect_lt(raw_resid_norm, 1e-6)
})
