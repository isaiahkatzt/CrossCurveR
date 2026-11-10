test_that("the DNS workflow retains baseline calibration and scalar shock results", {
  single <- make_feature_baseline()$single_curve_state
  fitted <- suppressWarnings(sc_fit_shock(
    single$yields, .35, c(2, 4), single$curves, single$mats
  ))
  expect_identical(fitted$baseline_state, single)
  expect_null(fitted$shock_profile)
  expect_identical(fitted$ns_factor, single$ns_factor)

  # A saved baseline must bypass calibration entirely.
  local_mocked_bindings(st_sc_state = function(...) stop("Unexpected refit"),
                       .package = "CrossCurveR")
  n <- length(single$time)
  for (factor in c("L", "S", "C")) {
    for (kind in c("additive", "multiplicative")) {
      expected <- st_dns_shock(single, "gbr", factor, .1, kind)
      actual <- sc_fit_shock(baseline_state = single, shock_curves = "gbr",
                            shock_factors = factor, shock_profile = rep(.1, n),
                            shock_type = kind)
      expected$shock_matrix <- NULL
      expected <- append(expected, list(shock_profile = rep(.1, n)),
                         after = match("ns_factor", names(expected)))
      expect_identical(actual, expected)
    }
  }
})

test_that("DNS profiles preserve metadata and reject invalid magnitudes", {
  single <- make_feature_baseline()$single_curve_state
  n <- length(single$time)
  profile <- data.frame(time = single$time, shock_magnitude = seq(0, .1, length.out = n))
  result <- sc_fit_shock(baseline_state = single, shock_curves = "gbr",
                         shock_profile = profile)
  expect_identical(result$shock_profile, profile)
  expect_identical(result$shocked_factors$usa, single$ns_factor$usa)
  expect_equal(result$shocked_factors$gbr[, "L"] - single$ns_factor$gbr[, "L"],
               profile$shock_magnitude)
  for (bad in c(NA_real_, Inf, -Inf, NaN)) {
    expect_error(sc_fit_shock(baseline_state = single, shock_curves = "gbr",
                             shock_profile = rep(bad, n)), "finite")
  }
  expect_error(sc_fit_shock(baseline_state = single, shock_curves = "gbr",
                           shock_profile = .1), "length")
})
