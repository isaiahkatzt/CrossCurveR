test_that("reused coefficients retain scenario-sample standardization", {
  baseline <- make_feature_baseline()
  result <- mc_fit_shock_cached(baseline, "gbr", "L", .25)
  # A constant spread translation cancels under the pre-existing normalization.
  expect_equal(result$Xt, baseline$xt_state$Xt, tolerance = 1e-12)
  expect_equal(result$sigma_JT, baseline$baseline_fit$sigma_JT, tolerance = 1e-12)
  expect_identical(result$BS0, baseline$BS0)
  expect_equal(result$xt_state$feature_spec$relationships$gbr$center,
    baseline$xt_state$feature_spec$relationships$gbr$center + c(.25, 0),
    tolerance = 1e-12)
  expect_equal(colMeans(result$Xt), c(gbr.level = 0, gbr.slope = 0),
    tolerance = 1e-12)

  raw <- make_feature_baseline(normalize = FALSE)
  shifted <- mc_fit_shock_cached(raw, "gbr", "L", .25)
  expected <- raw$xt_state$Xt
  expected[, "gbr.level"] <- expected[, "gbr.level"] + .25
  expect_equal(shifted$Xt, expected, tolerance = 1e-12)
})

test_that("zero shock preserves features and conditional covariance", {
  for (normalize in c(FALSE, TRUE)) {
    baseline <- make_feature_baseline(normalize = normalize)
    result <- mc_fit_shock_cached(baseline, "gbr", "L", 0)
    expect_equal(result$Xt, baseline$xt_state$Xt, tolerance = 1e-12)
    expect_equal(result$sigma_JT, baseline$baseline_fit$sigma_JT, tolerance = 1e-12)
  }
})

test_that("scenario feature controls remain independent of coefficient reuse", {
  baseline <- make_feature_baseline()
  raw <- mc_fit_shock_cached(baseline, "gbr", "L", .25,
    X_normalize = FALSE)
  expect_identical(raw$BS0, baseline$BS0)
  expect_false(raw$xt_state$feature_spec$normalized)
  expect_equal(colMeans(raw$Xt),
    baseline$xt_state$feature_spec$relationships$gbr$center + c(.25, 0),
    ignore_attr = TRUE, tolerance = 1e-12)

  smooth <- mc_fit_shock_cached(baseline, "gbr", "L", .25,
    Xt_smooth = TRUE, smoother = "hp", hp_lambda = 20)
  expect_identical(smooth$BS0, baseline$BS0)
  expect_identical(smooth$xt_state$feature_spec$smoothing$method, "hp")
  expect_equal(smooth$xt_state$feature_spec$smoothing$controls$hp_lambda, 20)
  expect_gt(max(abs(smooth$Xt - baseline$xt_state$Xt)), 1e-6)
})

test_that("the ECM can be re-estimated with covariance coefficients retained or refitted", {
  baseline <- make_feature_baseline()
  ecm <- baseline$xt_state$cc_ci_features$cc_ecm$gbr
  ecm$Gmatrix <- ecm$Gmatrix %*% diag(c(-1, 1))
  colnames(ecm$Gmatrix) <- c("level", "slope")
  with_temp_bindings(list(blmc_ECM = function(...) ecm), {
    retained <- mc_fit_shock_cached(baseline, "gbr", "L", .25,
      reuse_ECM = FALSE, reuse_BS0 = TRUE)
    expect_identical(retained$BS0, baseline$BS0)
    expect_equal(retained$xt_state$feature_spec$relationships$gbr$basis, ecm$Gmatrix)
    expected <- baseline$xt_state$Xt
    expected[, "gbr.level"] <- -expected[, "gbr.level"]
    expect_equal(retained$Xt, expected, tolerance = 1e-12)
    expect_gt(max(abs(as.matrix(retained$sigma_JT$jSigma[[1]][[1]]) -
      as.matrix(baseline$baseline_fit$sigma_JT$jSigma[[1]][[1]]))), 1e-6)

    set.seed(42)
    refitted <- quiet_covreg(mc_fit_shock_cached(baseline, "gbr", "L", .25,
      reuse_ECM = FALSE, reuse_BS0 = FALSE, CR_maxiter = 2))
    expect_equal(refitted$Xt, retained$Xt, tolerance = 1e-12)
    expect_false(identical(refitted$BS0, baseline$BS0))
  })
})

test_that("cached coefficient layouts are checked without requiring frozen transformations", {
  baseline <- make_feature_baseline()
  damaged <- baseline
  damaged$xt_state$feature_spec$relationships$gbr$scale <- c(1, 1)
  expect_error(mc_fit_shock_cached(damaged, "gbr", "L", .25), "stored input specification")
  damaged <- baseline
  damaged$BS0[[1]]$B <- damaged$BS0[[1]]$B[, 2:1]
  expect_error(mc_fit_shock_cached(damaged, "gbr", "L", .25), "coefficient identities/order")
  damaged <- baseline
  damaged$BS0[[1]]$B <- damaged$BS0[[1]]$B[, 1, drop = FALSE]
  expect_error(mc_fit_shock_cached(damaged, "gbr", "L", .25), "coefficient dimensions")
  damaged <- baseline
  damaged$BS0 <- rev(damaged$BS0)
  expect_error(mc_fit_shock_cached(damaged, "gbr", "L", .25), "maturity order")

  legacy <- baseline
  legacy$input_spec <- NULL
  legacy$xt_state$feature_spec <- NULL
  result <- mc_fit_shock_cached(legacy, "gbr", "L", .25)
  expect_equal(result$Xt, baseline$xt_state$Xt, tolerance = 1e-12)
  expect_identical(result$BS0, baseline$BS0)
})

test_that("cached convenience calls inherit settings and honor explicit feature overrides", {
  baseline <- make_feature_baseline(normalize = FALSE, trunc = TRUE)
  direct <- mc_fit_shock_cached(baseline, "gbr", "S", .25)
  wrapped <- mc_fit_shock(baseline_state = baseline, shock_curves = "gbr",
    shock_factors = "S", shock_magnitude = .25)
  expect_identical(wrapped, direct)
  overridden <- mc_fit_shock(baseline_state = baseline, shock_curves = "gbr",
    shock_factors = "S", shock_magnitude = .25, X_normalize = TRUE)
  expect_identical(overridden$BS0, baseline$BS0)
  expect_equal(mean(overridden$Xt), 0, tolerance = 1e-12)
  expect_identical(overridden, mc_fit_shock_cached(baseline, "gbr", "S", .25,
    X_normalize = TRUE))

  captured <- NULL
  with_temp_bindings(list(mc_fit_shock_cached = function(..., ECM_alpha = NULL,
      CR_maxiter = NULL, X_normalize = NULL, X_trunc = NULL) {
    captured <<- list(ECM_alpha, CR_maxiter, X_normalize, X_trunc)
  }), mc_fit_shock(baseline_state = baseline, shock_curves = "gbr", shock_factors = "S",
    shock_magnitude = .25, ECM_alpha = .05, CR_maxiter = 7))
  expect_identical(captured, list(.05, 7, NULL, NULL))
})

test_that("cached wrappers reject conflicting baseline-only inputs", {
  baseline <- make_feature_baseline()
  call <- list(baseline_state = baseline, shock_curves = "gbr", shock_factors = "L",
               shock_magnitude = .25)
  expect_error(do.call(mc_fit_shock, c(call, list(lambdas = .5))), "lambdas")
  expect_error(do.call(mc_fit_shock, c(call, list(cutoffs = c(1, 3)))), "cutoffs")
  yields <- baseline$yields
  yields$gbr[1, 2] <- yields$gbr[1, 2] + .1
  expect_error(do.call(mc_fit_shock, c(call, list(yields = yields))), "Supplied yields")
  expect_equal(do.call(mc_fit_shock, c(call, list(yields = baseline$yields,
    lambdas = baseline$config$lambdas, cutoffs = baseline$config$cutoffs)))$Xt,
    mc_fit_shock_cached(baseline, "gbr", "L", .25)$Xt)
})

test_that("fresh convenience calls preserve ordinary baseline defaults", {
  fixture <- make_feature_baseline()
  ecm <- fixture$xt_state$cc_ci_features$cc_ecm$gbr
  with_temp_bindings(list(blmc_ECM = function(...) ecm), {
    set.seed(31)
    baseline <- quiet_h(quiet_covreg(mc_fit_shock_base(fixture$yields, .35, c(2, 4), "usa",
      curves = fixture$curves, mats = fixture$mats, CR_maxiter = 2)))
    set.seed(31)
    wrapped <- quiet_h(quiet_covreg(mc_fit_shock(fixture$yields, .35, c(2, 4), "usa",
      shock_curves = "gbr", shock_factors = "L", shock_magnitude = .25,
      curves = fixture$curves, mats = fixture$mats, CR_maxiter = 2)))
    direct <- mc_fit_shock_cached(baseline, "gbr", "L", .25)
    expect_identical(wrapped, direct)
    expect_equal(wrapped$baseline_state$config$ECM_alpha, .1)
    expect_identical(wrapped$baseline_state$config$X_trunc, FALSE)
    expect_identical(wrapped$baseline_state$config$CR_store_path, TRUE)
  })
})

test_that("baseline fallback covers errors during ECM estimation", {
  baseline <- make_feature_baseline()
  with_temp_bindings(list(blcc_cointegration = function(...) stop("VECM estimation failed")), {
    expect_warning(result <- st_xt_state(baseline$single_curve_state$bl_yields,
      baseline$curves, baseline$reference, reuse_ECM = FALSE,
      baseline_ecm = baseline$xt_state$cc_ci_features$cc_ecm,
      ECM_fallback = "baseline", ytime = baseline$single_curve_state$time), "falling back")
    expect_equal(result$Xt, baseline$xt_state$Xt, tolerance = 1e-12)
    expect_match(result$cc_ci_features$ECM_fallback_reason, "VECM estimation failed")
  })
})
