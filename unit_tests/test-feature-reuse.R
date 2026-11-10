test_that("application profile scenarios retain the core scenario-sample convention", {
  baseline <- make_feature_baseline()
  spec <- list(shock_curves = "gbr", shock_factors = "L")
  n <- nrow(baseline$xt_state$Xt)
  for (window in c(n, n / 2)) {
    profile <- data.frame(shock_magnitude = c(rep(0, n - window), rep(.25, window)))
    application <- tontine_mce_shock_from_baseline_profile(baseline, spec, profile)
    core <- mc_fit_shock_cached(baseline, "gbr", "L", .25,
      shock_window = window, shock_window_position = "last")
    expect_identical(application$Xt, core$Xt)
    expect_identical(application$BS0, core$BS0)
    expect_identical(application$shocked_curve, core$shocked_curve)
    expect_identical(application$xt_state$feature_spec, core$xt_state$feature_spec)
    core$shock_matrix <- NULL
    core$input_spec <- NULL
    expected <- append(core, list(shock_profile = profile), after = match("ns_factor", names(core)))
    expect_identical(application, expected)
  }
  raw <- tontine_mce_shock_from_baseline_profile(baseline, spec, profile,
    controls = list(X_normalize = FALSE))
  expect_identical(raw$BS0, baseline$BS0)
  expect_false(raw$xt_state$feature_spec$normalized)
})

test_that("application profiles delegate fitting while preserving defaults and profile metadata", {
  baseline <- make_feature_baseline()
  profile <- data.frame(time = baseline$single_curve_state$time,
                        shock_magnitude = seq(0, .1, length.out = nrow(baseline$xt_state$Xt)))
  spec <- list(shock_curves = "gbr", shock_factors = "L", shock_basis = "average_level_fraction")
  controls <- list(X_normalize = FALSE)
  original <- CrossCurveR::mc_fit_shocked
  calls <- list()
  with_temp_bindings(list(mc_fit_shocked = function(baseline_state, shocked_factors,
      reuse_BS0, reuse_ECM, ECM_fallback, controls) {
    calls[[length(calls) + 1L]] <<- list(reuse_BS0, reuse_ECM, ECM_fallback, controls)
    original(baseline_state, shocked_factors, reuse_BS0, reuse_ECM, ECM_fallback, controls)
  }), result <- tontine_mce_shock_from_baseline_profile(baseline, spec, profile, controls = controls))
  expect_identical(calls, list(list(TRUE, TRUE, "baseline", controls)))
  expect_identical(result$shock_profile$time, profile$time)
  expect_identical(result$shock_profile$shock_fraction, profile$shock_magnitude)
  expect_false("input_spec" %in% names(result))
  expect_false(result$xt_state$feature_spec$normalized)
})
