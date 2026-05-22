make_fake_shock_baseline <- function() {
  structure(list(
    curves = fixture_curves,
    mats = fixture_mats,
    reference = "usa",
    yields = make_pair_yields(),
    single_curve_state = list(
      curves = fixture_curves,
      mats = fixture_mats,
      time = fixture_time,
      yields = make_pair_yields(),
      bl_yields = list(),
      phi_hat = matrix(1, nrow = length(fixture_curves) * length(fixture_mats), ncol = 3),
      W = list(X01M = matrix(0, nrow = length(fixture_time), ncol = length(fixture_curves))),
      ns_factor = list(
        usa = cbind(L = seq_along(fixture_time), S = seq_along(fixture_time) / 10),
        gbr = cbind(L = seq_along(fixture_time) + 1, S = seq_along(fixture_time) / 10 + 0.1)
      )
    ),
    xt_state = list(
      Xt = matrix(1, nrow = length(fixture_time), ncol = 1, dimnames = list(NULL, "gbr.r1")),
      cc_ci_features = list(
        cc_ecm = list(gbr = list(rank = 1, Gmatrix = matrix(1, nrow = 1, ncol = 1)))
      )
    ),
    BS0 = list(X01M = list(S0 = diag(2), B = matrix(1, nrow = 2, ncol = 1))),
    baseline_fit = list(
      tenor = data.frame(time = fixture_time, value = seq_along(fixture_time)),
      curve = data.frame(time = fixture_time, value = seq_along(fixture_time))
    ),
    config = list(
      lambdas = c(0.1, 0.2),
      cutoffs = c(1, 3),
      ECM_estim = "ML",
      ECM_type = "eigen",
      ECM_alpha = 0.1,
      X_normalize = TRUE,
      X_trunc = TRUE,
      CR_algo = zero_mean_covreg_em,
      CR_init = "adaptive",
      Xt_smooth = FALSE,
      smoother = NULL,
      knot_count = 5,
      k_count = 9,
      k_pass = 3,
      hp_lambda = 1600,
      henderson_k = 13,
      CR_maxiter = 1000,
      CR_tol = 1e-8,
      CR_Binit = NULL,
      CR_S0init = NULL,
      CR_S0_shrink_diag = 0.4,
      CR_verb = FALSE,
      CR_term = TRUE
    )
  ), class = "mc_shock_baseline")
}

test_that("mc_fit_shock_baseline builds and returns a reusable baseline state", {
  fake_single_curve_state <- make_fake_shock_baseline()$single_curve_state
  fake_xt_state <- make_fake_shock_baseline()$xt_state
  fake_BS0 <- make_fake_shock_baseline()$BS0
  fake_fit <- list(
    tenor = data.frame(time = fixture_time, value = seq_along(fixture_time)),
    curve = data.frame(time = fixture_time, value = seq_along(fixture_time))
  )

  baseline <- with_temp_bindings(
    list(
      st_prepare_single_curve_state = function(...) fake_single_curve_state,
      st_build_xt_state = function(...) fake_xt_state,
      blcc_fe = function(...) fake_BS0,
      st_finalize_shocked_fit = function(...) fake_fit
    ),
    mc_fit_shock_baseline(
      yields = make_pair_yields(),
      lambdas = c(0.1, 0.2),
      cutoffs = c(1, 3),
      reference = "usa",
      curves = fixture_curves,
      mats = fixture_mats,
      X_trunc = TRUE,
      CR_S0_shrink_diag = 0.4
    )
  )

  expect_s3_class(baseline, "mc_shock_baseline")
  expect_equal(baseline$reference, "usa")
  expect_equal(baseline$xt_state$Xt, fake_xt_state$Xt)
  expect_equal(baseline$baseline_fit$curve, fake_fit$curve)
  expect_equal(baseline$config$CR_S0_shrink_diag, 0.4)
})

test_that("mc_fit_shock_from_baseline uses cached DNS state instead of rebuilding baseline curves", {
  baseline_state <- make_fake_shock_baseline()
  captured_dns <- NULL

  shock_result <- with_temp_bindings(
    list(
      lambda_grid_search = function(...) stop("baseline lambda search should not run"),
      bl_single_curve = function(...) stop("baseline single-curve fit should not run"),
      st_additive = function(curves, curve_dns_factor, ...) {
        captured_dns <<- curve_dns_factor
        list(
          shock_matrix = matrix(c(0.2, 0, NA, NA), nrow = 2, dimnames = list(c("L", "S"), curves)),
          shocked_factors = lapply(curve_dns_factor, function(x) {
            out <- x
            out[, "L"] <- out[, "L"] + 0.2
            out
          })
        )
      },
      st_build_shocked_curve_state = function(single_curve_state, shocked_factors) {
        list(
          yields = single_curve_state$yields,
          bl_yields = list(),
          phi_hat = single_curve_state$phi_hat,
          W = single_curve_state$W,
          ns_factor = shocked_factors
        )
      },
      st_build_xt_state = function(...) list(
        Xt = matrix(2, nrow = length(fixture_time), ncol = 1, dimnames = list(NULL, "gbr.r1")),
        cc_ci_features = list(cc_ecm = list(gbr = list(rank = 1, Gmatrix = matrix(1, nrow = 1, ncol = 1))))
      ),
      st_finalize_shocked_fit = function(yields, time, phi_hat, Xt, BS0, curves, mats, W) {
        list(
          tenor = data.frame(time = time, shocked = seq_along(time)),
          curve = data.frame(time = time, shocked = seq_along(time)),
          sigma_JT = "shock_sigma",
          Xt = Xt,
          W = W
        )
      }
    ),
    mc_fit_shock_from_baseline(
      baseline_state = baseline_state,
      shock_curves = "usa",
      shock_factors = "L",
      shock_magnitude = 0.2
    )
  )

  expect_equal(captured_dns, baseline_state$single_curve_state$ns_factor)
  expect_equal(shock_result$curve, baseline_state$baseline_fit$curve)
  expect_equal(shock_result$shocked_curve$shocked, seq_along(fixture_time))
  expect_equal(shock_result$baseline_state$reference, "usa")
})

test_that("st_build_xt_state can fall back to baseline ECM when a refit yields no usable Xt", {
  baseline_ecm <- list(gbr = list(rank = 1, Gmatrix = matrix(1, nrow = 1, ncol = 1)))
  captured_ecm <- NULL

  result <- with_temp_bindings(
    list(
      blcc_cointegration = function(...) list(
        cc_cspread = "refit_cspread",
        cc_ecm = "refit_ecm",
        cc_beta = "refit_beta"
      ),
      blcc_build_Xt = function(...) {
        stop("No usable cointegration relationships were identified for Xt construction.")
      },
      st_build_xt_from_ecm = function(bl_yields, curves, reference, curve_ecm,
                                      X_normalize, X_trunc) {
        captured_ecm <<- curve_ecm
        list(
          Xt = matrix(9, nrow = length(fixture_time), ncol = 1, dimnames = list(NULL, "gbr.r1")),
          cc_ci_features = list(cc_ecm = curve_ecm)
        )
      }
    ),
    expect_warning(
      st_build_xt_state(
        bl_yields = list(),
        curves = fixture_curves,
        reference = "usa",
        reuse_ECM = FALSE,
        baseline_ecm = baseline_ecm,
        ECM_fallback = "baseline",
        ytime = fixture_time
      ),
      "falling back to baseline ECM"
    )
  )

  expect_equal(captured_ecm, baseline_ecm)
  expect_equal(result$Xt, matrix(9, nrow = length(fixture_time), ncol = 1, dimnames = list(NULL, "gbr.r1")))
  expect_equal(result$cc_ci_features$cc_ecm, baseline_ecm)
  expect_equal(result$cc_ci_features$refit_cc_ci_features$cc_ecm, "refit_ecm")
  expect_equal(result$cc_ci_features$ECM_fallback, "baseline")
  expect_match(result$cc_ci_features$ECM_fallback_reason, "No usable cointegration")
})

test_that("mc_fit_shock reuses a supplied baseline_state", {
  baseline_state <- make_fake_shock_baseline()

  result <- with_temp_bindings(
    list(
      mc_fit_shock_baseline = function(...) stop("baseline should not be rebuilt"),
      mc_fit_shock_from_baseline = function(baseline_state, shock_curves, shock_factors, shock_magnitude, ...) {
        list(
          baseline_reference = baseline_state$reference,
          shock_curves = shock_curves,
          shock_factors = shock_factors,
          shock_magnitude = shock_magnitude
        )
      }
    ),
    mc_fit_shock(
      baseline_state = baseline_state,
      shock_curves = "usa",
      shock_factors = "L",
      shock_magnitude = 0.15
    )
  )

  expect_equal(result$baseline_reference, "usa")
  expect_equal(result$shock_curves, "usa")
  expect_equal(result$shock_factors, "L")
  expect_equal(result$shock_magnitude, 0.15)
})

test_that("st_additive supports placing a finite shock window at the end of the sample", {
  curve_dns_factor <- list(
    usa = cbind(L = 1:5, S = seq(0.1, 0.5, by = 0.1)),
    gbr = cbind(L = 6:10, S = seq(0.6, 1, by = 0.1))
  )

  result <- st_additive(
    curves = fixture_curves,
    curve_dns_factor = curve_dns_factor,
    shock_curves = "usa",
    shock_factors = "L",
    shock_magnitude = 0.5,
    shock_window = 2,
    shock_window_position = "last"
  )

  expect_equal(result$shocked_factors$usa[1:3, "L"], curve_dns_factor$usa[1:3, "L"])
  expect_equal(result$shocked_factors$usa[4:5, "L"], curve_dns_factor$usa[4:5, "L"] + 0.5)
  expect_equal(result$shocked_factors$gbr, curve_dns_factor$gbr)
})

test_that("shock row resolution selects full, leading, or trailing windows", {
  expect_equal(
    st_resolve_shock_rows(
      row_count = 10,
      shock_window = 3,
      shock_window_position = "last"
    ),
    8:10
  )
  expect_equal(
    st_resolve_shock_rows(
      row_count = 10,
      shock_window = NULL,
      shock_window_position = "first"
    ),
    1:10
  )
})

test_that("st_dns_shock_from_single_curve_state applies the requested factor shock only to the target curve", {
  single_curve_state <- make_fake_shock_baseline()$single_curve_state

  result <- with_temp_bindings(
    list(
      st_build_shocked_curve_state = function(single_curve_state, shocked_factors) {
        list(
          yields = single_curve_state$yields,
          ns_factor = shocked_factors
        )
      }
    ),
    st_dns_shock_from_single_curve_state(
      single_curve_state = single_curve_state,
      shock_curves = "usa",
      shock_factors = "L",
      shock_magnitude = 0.5
    )
  )

  expect_equal(result$shock_matrix["L", "usa"], 0.5)
  expect_equal(result$shock_matrix["L", "gbr"], 0)
  expect_equal(result$shocked_factors$usa[, "L"], single_curve_state$ns_factor$usa[, "L"] + 0.5)
  expect_equal(result$shocked_factors$gbr, single_curve_state$ns_factor$gbr)
})

test_that("st_dly_shock_from_fit_profile transmits shocks through the DLY global factor", {
  dly_fit_obj <- dly_fit(
    yields = make_pair_yields(),
    mats = fixture_mats,
    lambdas = c(0.35),
    curves = fixture_curves,
    reference = "usa",
    center = TRUE,
    scale. = FALSE
  )

  shock_profile <- rep(0, length(fixture_time))
  shock_profile[1:5] <- 0.25
  result <- st_dly_shock_from_fit_profile(
    dly_fit_obj = dly_fit_obj,
    shock_curves = "usa",
    shock_factors = "L",
    shock_profile = shock_profile
  )

  expected_shocked_single_curve <- st_build_shocked_dly_single_curve(
    dly_single_curve = dly_fit_obj$single_curve,
    shocked_factors = result$shocked_factors,
    curves = fixture_curves
  )
  expected_level_draw <- dly_group_factor(expected_shocked_single_curve, curves = fixture_curves, latent = "L")
  expected_slope_draw <- dly_group_factor(expected_shocked_single_curve, curves = fixture_curves, latent = "S")

  expected_global_level <- as.numeric(
    predict(dly_fit_obj$global_factors$level$pca, newdata = expected_level_draw)[, 1]
  )
  expected_global_slope <- as.numeric(
    predict(dly_fit_obj$global_factors$slope$pca, newdata = expected_slope_draw)[, 1]
  )
  shock_rows <- 1:5
  expected_shocked_level <- st_dly_predict_country_factor(
    factor_model = dly_fit_obj$factor_models$level,
    global_factor = expected_global_level,
    curves = fixture_curves
  )
  expected_shocked_slope <- st_dly_predict_country_factor(
    factor_model = dly_fit_obj$factor_models$slope,
    global_factor = expected_global_slope,
    curves = fixture_curves
  )
  expected_level_residual <- dly_fit_obj$factor_models$level$residual[, fixture_curves]
  expected_slope_residual <- dly_fit_obj$factor_models$slope$residual[, fixture_curves]
  expected_level_residual[, "usa"] <-
    dly_fit_obj$factor_models$level$residual[, "usa"] +
    expected_level_draw[, "usa"] - dly_fit_obj$country_factors$level[, "usa"]
  expected_slope_residual[, "usa"] <-
    dly_fit_obj$factor_models$slope$residual[, "usa"] +
    expected_slope_draw[, "usa"] - dly_fit_obj$country_factors$slope[, "usa"]
  expected_global_plus_idio_level <- expected_shocked_level + expected_level_residual
  expected_global_plus_idio_slope <- expected_shocked_slope + expected_slope_residual
  expected_global_plus_idio_betas <- dly_bind_betas(
    expected_global_plus_idio_level,
    expected_global_plus_idio_slope,
    curves = fixture_curves
  )

  expect_equal(as.matrix(result$shocked_factors$usa)[shock_rows, 1], dly_fit_obj$country_betas$raw$usa[shock_rows, "L"] + 0.25)
  expect_equal(as.matrix(result$shocked_factors$usa)[-shock_rows, 1], dly_fit_obj$country_betas$raw$usa[-shock_rows, "L"])
  expect_equal(unname(as.matrix(result$shocked_factors$gbr)), unname(dly_fit_obj$country_betas$raw$gbr))

  expect_equal(result$fit$global_factors$level$factor, expected_global_level)
  expect_equal(result$fit$global_factors$slope$factor, expected_global_slope)
  expect_equal(result$fit$factor_models$level$loading, dly_fit_obj$factor_models$level$loading)
  expect_equal(result$fit$factor_models$slope$loading, dly_fit_obj$factor_models$slope$loading)
  expect_equal(result$fit$country_betas$idiosyncratic$usa, cbind(L = expected_level_residual[, "usa"], S = expected_slope_residual[, "usa"]))
  expect_equal(result$fit$country_betas$idiosyncratic$gbr, dly_fit_obj$country_betas$idiosyncratic$gbr)
  expect_gt(
    max(abs(
      result$fit$country_betas$global_plus_idiosyncratic$usa -
        result$fit$country_betas$raw$usa
    )),
    0
  )
  expect_equal(
    result$fit$country_betas$global_plus_idiosyncratic$usa,
    expected_global_plus_idio_betas$usa,
    tolerance = 1e-8
  )
  expect_equal(
    result$fit$country_betas$global_plus_idiosyncratic$gbr,
    expected_global_plus_idio_betas$gbr,
    tolerance = 1e-8
  )
  expect_gt(
    max(abs(
      result$fit$country_betas$global_plus_idiosyncratic$gbr -
        dly_fit_obj$country_betas$global_plus_idiosyncratic$gbr
    )),
    0
  )
  expect_equal(
    as.matrix(result$fit$yields$global_plus_idiosyncratic$usa[, -1, drop = FALSE]),
    as.matrix(dly_curve_yields(
      expected_shocked_single_curve,
      expected_global_plus_idio_level,
      expected_global_plus_idio_slope,
      curves = fixture_curves
    )$usa[, -1, drop = FALSE]),
    tolerance = 1e-8
  )
})

test_that("st_dly_global_shock_from_fit_profile shifts average fitted country level by target profile", {
  dly_fit_obj <- dly_fit(
    yields = make_pair_yields(),
    mats = fixture_mats,
    lambdas = c(0.35),
    curves = fixture_curves,
    reference = "usa",
    center = TRUE,
    scale. = FALSE
  )

  shock_fraction <- rep(-0.3, length(fixture_time))
  target_profile <- st_average_level_shock_profile(
    curve_factor = dly_fit_obj$country_betas$raw,
    shock_fraction_profile = shock_fraction,
    shock_factor = "L",
    curves = fixture_curves
  )

  result <- st_dly_global_shock_from_fit_profile(
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

test_that("mc_fit_shock_from_baseline uses shocked factors without refitting single-curve state", {
  baseline_state <- make_fake_shock_baseline()
  shocked_yields <- make_pair_yields()
  direct_shocked_factors <- baseline_state$single_curve_state$ns_factor
  captured_bl_yields <- NULL
  captured_W <- NULL

  result <- with_temp_bindings(
    list(
      st_additive = function(curves, curve_dns_factor, ...) {
        list(
          shock_matrix = matrix(c(0.2, 0, NA, NA), nrow = 2, dimnames = list(c("L", "S"), curves)),
          shocked_factors = direct_shocked_factors
        )
      },
      st_build_shocked_curve_state = function(single_curve_state, shocked_factors) {
        list(
          yields = shocked_yields,
          bl_yields = "shock_bl_yields",
          phi_hat = "shock_phi_hat",
          W = "shock_W",
          ns_factor = "shock_ns_factor"
        )
      },
      st_prepare_single_curve_state = function(yields, lambdas, cutoffs, curves, mats) {
        stop("single-curve state should not be refit during MCE shocks")
      },
      st_build_xt_state = function(bl_yields, curves, reference, reuse_ECM, baseline_ecm, ..., ytime) {
        captured_bl_yields <<- bl_yields
        expect_true(reuse_ECM)
        expect_equal(reference, "usa")
        expect_equal(baseline_ecm, baseline_state$xt_state$cc_ci_features$cc_ecm)
        expect_equal(ytime, fixture_time)
        list(
          Xt = matrix(2, nrow = length(fixture_time), ncol = 1, dimnames = list(NULL, "gbr.r1")),
          cc_ci_features = list(cc_ecm = list(gbr = list(rank = 1, Gmatrix = matrix(1, nrow = 1, ncol = 1))))
        )
      },
      blcc_fe = function(Xt, Wt, mats, ...) {
        captured_W <<- Wt
        "refit_BS0"
      },
      st_finalize_shocked_fit = function(yields, time, phi_hat, Xt, BS0, curves, mats, W) {
        expect_equal(yields, shocked_yields)
        expect_equal(time, fixture_time)
        expect_equal(phi_hat, "shock_phi_hat")
        expect_equal(BS0, "refit_BS0")
        expect_equal(W, "shock_W")
        list(
          tenor = data.frame(time = time, shocked = seq_along(time)),
          curve = data.frame(time = time, shocked = seq_along(time)),
          sigma_JT = "shock_sigma",
          Xt = Xt,
          W = W
        )
      }
    ),
    mc_fit_shock_from_baseline(
      baseline_state = baseline_state,
      shock_curves = "usa",
      shock_factors = "L",
      shock_magnitude = 0.2,
      reuse_ECM = TRUE,
      reuse_BS0 = FALSE
    )
  )

  expect_equal(captured_bl_yields, "shock_bl_yields")
  expect_equal(captured_W, "shock_W")
  expect_equal(result$curve, baseline_state$baseline_fit$curve)
  expect_equal(result$shocked_curve$shocked, seq_along(fixture_time))
  expect_equal(result$shock_matrix["L", "usa"], 0.2)
  expect_equal(result$ns_factor, "shock_ns_factor")
  expect_equal(result$shocked_factors, direct_shocked_factors)
  expect_equal(result$shocked_yields, shocked_yields)
  expect_equal(result$shocked_single_curve_state$ns_factor, "shock_ns_factor")
})

test_that("mc_fit_shock_refit_from_baseline delegates to MCE refits without single-curve refit", {
  baseline_state <- make_fake_shock_baseline()
  captured <- NULL

  result <- with_temp_bindings(
    list(
      mc_fit_shock_from_baseline = function(baseline_state, shock_curves, shock_factors,
                                            shock_magnitude, shock_type, shock_window,
                                            shock_window_position,
                                            reuse_BS0, reuse_ECM,
                                            ...) {
        captured <<- list(
          shock_curves = shock_curves,
          shock_factors = shock_factors,
          shock_magnitude = shock_magnitude,
          shock_type = shock_type,
          shock_window = shock_window,
          shock_window_position = shock_window_position,
          reuse_BS0 = reuse_BS0,
          reuse_ECM = reuse_ECM
        )
        list(tag = "delegated", baseline_reference = baseline_state$reference)
      }
    ),
    mc_fit_shock_refit_from_baseline(
      baseline_state = baseline_state,
      shock_curves = "usa",
      shock_factors = "L",
      shock_magnitude = 0.2
    )
  )

  expect_equal(result$tag, "delegated")
  expect_equal(result$baseline_reference, "usa")
  expect_equal(captured$shock_curves, "usa")
  expect_equal(captured$shock_factors, "L")
  expect_equal(captured$shock_magnitude, 0.2)
  expect_equal(captured$shock_type, "additive")
  expect_null(captured$shock_window)
  expect_equal(captured$shock_window_position, "first")
  expect_false(captured$reuse_BS0)
  expect_false(captured$reuse_ECM)
})
