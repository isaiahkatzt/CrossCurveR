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

test_that("st_dly_shock_from_fit rebuilds a shocked DLY fit from shocked single-curve factors", {
  baseline_betas <- list(
    usa = cbind(L = c(1, 2), S = c(0.2, 0.3)),
    gbr = cbind(L = c(3, 4), S = c(0.4, 0.5))
  )
  captured_shocked_factors <- NULL

  dly_fit_obj <- list(
    curves = fixture_curves,
    reference = "usa",
    center = TRUE,
    scale. = FALSE,
    mats = fixture_mats,
    single_curve = "baseline_single_curve"
  )

  result <- with_temp_bindings(
    list(
      dly_pull_beta = function(...) baseline_betas,
      st_build_shocked_dly_single_curve = function(dly_single_curve, shocked_factors, curves) {
        captured_shocked_factors <<- shocked_factors
        list(tag = "shocked_single_curve")
      },
      dly_fit_from_single_curve = function(dly_sc, curves, reference, center, scale., mats) {
        list(
          tag = "rebuilt_dly",
          dly_sc = dly_sc,
          curves = curves,
          reference = reference,
          center = center,
          scale. = scale.,
          mats = mats
        )
      }
    ),
    st_dly_shock_from_fit(
      dly_fit_obj = dly_fit_obj,
      shock_curves = "usa",
      shock_factors = "L",
      shock_magnitude = 0.25
    )
  )

  expect_equal(result$shock_matrix["L", "usa"], 0.25)
  expect_equal(result$shock_matrix["L", "gbr"], 0)
  expect_equal(captured_shocked_factors$usa[, "L"], baseline_betas$usa[, "L"] + 0.25)
  expect_equal(captured_shocked_factors$gbr, baseline_betas$gbr)
  expect_equal(result$fit$tag, "rebuilt_dly")
  expect_equal(result$fit$reference, "usa")
  expect_equal(result$fit$mats, fixture_mats)
})
