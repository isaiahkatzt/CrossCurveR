test_that("application profile scenarios retain the core scenario-sample convention", {
  baseline <- make_feature_baseline()
  spec <- list(shock_curves = "gbr", shock_factors = "L")
  n <- nrow(baseline$xt_state$Xt)
  for (window in c(n, n / 2)) {
    profile <- data.frame(shock_magnitude = c(rep(0, n - window), rep(.25, window)))
    application <- tontine_mce_shock_from_baseline_profile(baseline, spec, profile)
    core <- mc_fit_shock_from_baseline(baseline, "gbr", "L", .25,
      shock_window = window, shock_window_position = "last")
    expect_identical(application$Xt, core$Xt)
    expect_identical(application$BS0, core$BS0)
    expect_identical(application$shocked_curve, core$shocked_curve)
    expect_identical(application$xt_state$feature_spec, core$xt_state$feature_spec)
  }
  raw <- tontine_mce_shock_from_baseline_profile(baseline, spec, profile,
    controls = list(X_normalize = FALSE))
  expect_identical(raw$BS0, baseline$BS0)
  expect_false(raw$xt_state$feature_spec$normalized)
})
