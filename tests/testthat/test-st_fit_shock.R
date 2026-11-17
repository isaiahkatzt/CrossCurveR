test_that("exported shock names preserve their compatibility aliases", {
  aliases <- c(
    st_apply_factor_shock_profile = "st_factor_profile",
    st_average_level_shock_profile = "st_mean_profile",
    st_dly_global_shock_from_fit_profile = "st_dly_global_shock",
    st_dly_shock_from_fit_profile = "st_dly_shock",
    mc_fit_shock_baseline = "mc_fit_shock_base",
    mc_fit_shock_from_baseline = "mc_fit_shock_cached",
    mc_fit_shock_refit_from_baseline = "mc_fit_shock_refit"
  )
  for (old in names(aliases)) {
    expect_identical(getExportedValue("CrossCurveR", old),
                     getExportedValue("CrossCurveR", aliases[[old]]))
  }
})

test_that("cached baselines reject inconsistent state and conflicting inputs", {
  baseline <- make_feature_baseline()
  expect_error(st_check_baseline(baseline, curves = rev(baseline$curves)), "Supplied .*curves")
  expect_error(st_check_baseline(baseline, mats = baseline$mats[-1]), "Supplied .*mats")
  expect_error(st_check_baseline(baseline, reference = "gbr"), "Supplied .*reference")

  wrong <- baseline
  wrong$BS0 <- NULL
  expect_error(st_check_baseline(wrong), "must be created")
  wrong <- baseline
  wrong$xt_state$Xt <- wrong$xt_state$Xt[-1, , drop = FALSE]
  expect_error(st_check_baseline(wrong), "one row")
  wrong <- baseline
  wrong$input_spec$features$columns[1] <- "different"
  expect_error(st_check_baseline(wrong), "stored input specification")
  wrong <- baseline
  wrong$single_curve_state$yields$gbr$time <- wrong$single_curve_state$time + 1
  expect_error(st_check_baseline(wrong), "observation index")
})

test_that("cached MCE shocks keep the calibrated curve state and residuals", {
  baseline <- make_feature_baseline()
  with_temp_bindings(list(st_sc_state = function(...) stop("baseline should not be rebuilt")), {
    result <- mc_fit_shock(baseline_state = baseline, shock_curves = "gbr",
                            shock_factors = "L", shock_magnitude = .2)
  })
  expect_identical(result$W, baseline$single_curve_state$W)
  expect_equal(result$ns_factor$gbr[, "L"], baseline$single_curve_state$ns_factor$gbr[, "L"] + .2)
})

test_that("an ECM refit with no relationships falls back to the baseline features", {
  baseline <- make_feature_baseline()
  refit <- list(cc_cspread = list(gbr = NULL), cc_ecm = list(gbr = list(rank = 0)))
  with_temp_bindings(list(blcc_cointegration = function(...) refit), {
    expect_warning(result <- st_xt_state(
      baseline$single_curve_state$bl_yields, baseline$curves, baseline$reference,
      reuse_ECM = FALSE, baseline_ecm = baseline$xt_state$cc_ci_features$cc_ecm,
      ECM_fallback = "baseline", ytime = baseline$single_curve_state$time
    ), "falling back to baseline ECM")
    expect_equal(result$Xt, baseline$xt_state$Xt)
    expect_equal(result$cc_ci_features$refit_cc_ci_features, refit)
    expect_match(result$cc_ci_features$ECM_fallback_reason, "no usable relationships")
  })
})

test_that("windowed and profile shocks match factor names in reordered panels", {
  for (factors in list(c("L", "S"), c("L", "S", "C"))) {
    x <- matrix(seq_len(5 * length(factors)), 5, dimnames = list(NULL, factors))
    input <- list(usa = x, gbr = x + 20)
    shuffled <- lapply(input, function(m) m[, rev(factors), drop = FALSE])
    target <- tail(factors, 1)
    expected <- input
    expected$usa[, target] <- x[, target] + .2
    expect_equal(st_factor_profile(names(input), shuffled, "usa", target,
                                   rep(.2, 5))$shocked_factors, expected)

    expected <- input
    expected$usa[4:5, target] <- x[4:5, target] * 1.2
    result <- st_factor_shock(names(input), shuffled, "usa", target, .2,
                              "multiplicative", 2, "last")
    expect_equal(result$shocked_factors, expected)
    expect_equal(st_mean_profile(shuffled, rep(.1, 5)), rowMeans(cbind(x[, "L"], x[, "L"] + 20)) * .1)
  }
})

test_that("shock factor labels reject ambiguity and retain positional compatibility", {
  x <- cbind(L = 1:3, S = 4:6, C = 7:9)
  for (labels in list(c("L", "L", "C"), c("L", "l1", "C"), c("L", "S", "unknown"),
                      c("L", "", "C"), c("L", NA, "C"))) {
    bad <- x
    colnames(bad) <- labels
    input <- list(usa = x, gbr = bad) # Validate even an unshocked curve.
    expect_error(st_factor_profile(names(input), input, "usa", "L", rep(.1, 3)),
                 "Factor columns")
  }
  input <- list(usa = x, gbr = x[, c(1, 1, 3)])
  expect_error(st_factor_shock(names(input), input, "usa", "L", .1), "Factor columns")
  expect_error(st_mean_profile(input, rep(.1, 3)), "Factor columns")
  alias <- x[, c("C", "L", "S")]
  colnames(alias) <- c("l3", "l1", "l2")
  expect_identical(.st_factor_matrix(alias), x)
  unnamed <- x
  colnames(unnamed) <- NULL
  expect_identical(.st_factor_matrix(unnamed), unnamed)
  expect_error(.st_factor_matrix(x[, c("L", "C")]), "Factor columns")
  expect_error(st_factor_profile(c("usa", "gbr"),
    list(usa = x, gbr = x[, 1:2]), "usa", "L", rep(.1, 3)), "same number of columns")
})

test_that("additive and multiplicative shocks preserve factors outside the window", {
  factors <- list(usa = cbind(L = 1:5, S = 6:10), gbr = cbind(L = 11:15, S = 16:20))
  additive <- st_factor_shock(fixture_curves, factors, "usa", "L", .5,
                               shock_window = 2, shock_window_position = "last")
  expected <- factors
  expected$usa[4:5, "L"] <- factors$usa[4:5, "L"] + .5
  expect_equal(additive$shocked_factors, expected)

  multiplied <- st_multiplicative(fixture_curves, factors, "usa", "L", .1,
                                  shock_window = 2, shock_window_position = "first")
  expected <- factors
  expected$usa[1:2, "L"] <- factors$usa[1:2, "L"] * 1.1
  expect_equal(multiplied$shocked_factors, expected)
  expect_equal(multiplied$shock_matrix,
               rbind(L = c(usa = 1.1, gbr = 1), S = c(usa = 1, gbr = 1)))
})

test_that("DLY country shocks propagate through the frozen global factor", {
  baseline <- dly_fit(make_dly_pair_yields(), fixture_mats, lambdas = .35, reference = "usa")
  profile <- c(rep(.25, 5), rep(0, length(fixture_time) - 5))
  result <- st_dly_shock(baseline, "usa", "L", profile)
  fit <- result$fit
  level <- baseline$global_factors$level
  score_shift <- profile * level$pca$rotation["usa", 1]
  expect_equal(fit$global_factors$level$factor, level$factor + score_shift)
  expect_equal(fit$global_factors$slope$factor, baseline$global_factors$slope$factor)
  expect_identical(fit$factor_models$level$loading, baseline$factor_models$level$loading)

  expected_raw <- baseline$country_betas$raw
  expected_raw$usa[, "L"] <- expected_raw$usa[, "L"] + profile
  expect_equal(fit$country_betas$raw, expected_raw)
  for (curve in fixture_curves) {
    direct <- if (curve == "usa") profile else rep(0, length(profile))
    level_shift <- score_shift * baseline$factor_models$level$loading[[curve]] + direct
    before <- baseline$yields$global_plus_idiosyncratic[[curve]]
    after <- fit$yields$global_plus_idiosyncratic[[curve]]
    expect_equal(as.matrix(after[-1] - before[-1]),
                 outer(level_shift, rep(1, length(fixture_mats))),
                 ignore_attr = TRUE, tolerance = 1e-8)
  }
})

test_that("st_dly_global_shock_from_fit_profile shifts average fitted country level by target profile", {
  dly_fit_obj <- dly_fit(
    yields = make_dly_pair_yields(),
    mats = fixture_mats,
    lambdas = c(0.35),
    curves = fixture_curves,
    reference = "usa",
    center = TRUE,
    scale. = FALSE
  )

  shock_fraction <- rep(-0.3, length(fixture_time))
  target_profile <- st_mean_profile(
    curve_factor = dly_fit_obj$country_betas$raw,
    shock_fraction_profile = shock_fraction,
    shock_factor = "L",
    curves = fixture_curves
  )

  result <- st_dly_global_shock(
    dly_fit_obj = dly_fit_obj,
    shock_factors = "L",
    shock_profile = target_profile
  )

  baseline_avg_level <- rowMeans(do.call(
    cbind,
    lapply(dly_fit_obj$country_betas$global_plus_idiosyncratic, function(x) x[, "L"])
  ))
  shocked_avg_level <- rowMeans(do.call(
    cbind,
    lapply(result$fit$country_betas$global_plus_idiosyncratic, function(x) x[, "L"])
  ))

  expect_equal(
    as.numeric(shocked_avg_level - baseline_avg_level),
    as.numeric(target_profile),
    tolerance = 1e-8
  )
})

test_that("the MCE refit entry point refits both covariance and ECM", {
  baseline <- make_feature_baseline()
  settings <- with_temp_bindings(
    list(mc_fit_shock_cached = function(..., reuse_BS0, reuse_ECM) c(reuse_BS0, reuse_ECM)),
    mc_fit_shock_refit(baseline, "usa", "L", .2)
  )
  expect_identical(settings, c(FALSE, FALSE))
})

test_that("shock matrices, profiles, and global shifts reject nonfinite values", {
  factors <- list(usa = cbind(L = seq_along(fixture_time), S = 1))
  fit <- dly_fit(make_dly_pair_yields(), mats = fixture_mats, lambdas = .35)
  profile <- rep(0, length(fixture_time))
  profile[2] <- Inf
  expect_error(st_build_shock_matrix(fixture_curves, c("L", "S"), "usa", "L", NA_real_), "finite")
  expect_error(st_factor_profile("usa", factors, "usa", "L", profile), "finite")
  expect_error(st_mean_profile(factors, profile), "finite")
  expect_error(st_dly_global_shock(fit, "L", profile), "finite")
})

test_that("st_build_shock_matrix validates inputs and places named shocks", {
  shock_matrix <- st_build_shock_matrix(
    curves = fixture_curves,
    model_factors = c("L", "S", "C"),
    shock_curves = "gbr",
    shock_factors = c("L", "C"),
    magnitude = c(C = -0.2, L = 0.5)
  )

  expect_equal(shock_matrix, cbind(usa = c(L = NA_real_, S = NA_real_, C = NA_real_),
                                   gbr = c(L = .5, S = NA_real_, C = -.2)))

  expect_error(
    st_build_shock_matrix(fixture_curves, c("L", "S", "C"), "usa", c("L", "S"), 0.1),
    "Scalar `magnitude`"
  )
  expect_error(
    st_build_shock_matrix(fixture_curves, c("L", "S", "C"), "cad", "L", 0.1),
    "`shock_curves` must be a subset"
  )
  expect_error(
    st_build_shock_matrix(c("usa", "usa"), c("L", "S", "C"), "usa", "L", 0.1),
    "`curves` must contain unique"
  )
})

test_that("st_build_shocked_curve_state rebuilds curve panels from shocked factors", {
  single_curve_state <- st_sc_state(
    yields = make_pair_yields(),
    lambdas = c(0.35),
    cutoffs = c(2, 4),
    curves = fixture_curves,
    mats = fixture_mats
  )
  shocked_factors <- single_curve_state$ns_factor
  shocked_factors$usa[, "L"] <- shocked_factors$usa[, "L"] + 0.25

  shocked_state <- st_shock_state(
    single_curve_state = single_curve_state,
    shocked_factors = shocked_factors
  )

  expected_usa_fit <- shocked_factors$usa %*% single_curve_state$bl_yields$phi$usa$phi

  expect_equal(names(shocked_state$yields), fixture_curves)
  expect_equal(colnames(shocked_state$yields$usa), colnames(single_curve_state$yields$usa))
  expect_equal(
    unname(as.matrix(shocked_state$yields$usa[, -1])),
    unname(expected_usa_fit),
    tolerance = 1e-8
  )
  expect_equal(shocked_state$ns_factor$usa, shocked_factors$usa)
  expect_equal(names(shocked_state$W), paste0("X", sapply(fixture_mats, numeric_to_matname)))
})

test_that("st_align_Xt_to_baseline reorders compatible columns and rejects incompatible designs", {
  Xt_new <- cbind(b = 1:3, a = 4:6)
  Xt_reference <- cbind(a = 10:12, b = 13:15)

  aligned <- st_align_Xt(Xt_new, Xt_reference)

  expect_identical(aligned, Xt_new[, c("a", "b")])
  expect_error(
    st_align_Xt(cbind(a = 1:3), Xt_reference),
    "incompatible"
  )
  expect_error(
    st_align_Xt(cbind(a = 1:3, c = 4:6), Xt_reference),
    "columns do not match"
  )
})
