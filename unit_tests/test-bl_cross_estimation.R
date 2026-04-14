test_that("blce_rescale_yield returns correctly labeled tenor and curve views", {
  normalized <- list(
    usa = as.matrix(make_ns_curve_df()[, -1]),
    gbr = as.matrix(make_ns_curve_df(curve_shift = 0.2)[, -1])
  )
  P <- bln_build_P(fixture_curves, fixture_mats)
  sqrt_sigma_t <- rep(list(diag(nrow(P))), nrow(normalized$usa))

  tenor_view <- blce_rescale_yield(normalized, sqrt_sigma_t, P, "tenor", fixture_time, fixture_curves, fixture_mats)
  curve_view <- blce_rescale_yield(normalized, sqrt_sigma_t, P, "curve", fixture_time, fixture_curves, fixture_mats)

  expect_equal(nrow(tenor_view), length(fixture_time))
  expect_equal(nrow(curve_view), length(fixture_time))
  expect_equal(colnames(curve_view)[1], "time")
  expect_true(any(grepl("\\.", colnames(tenor_view)[-1])))
})

test_that("blce_smooth_Xt dispatches to spline, RM, HP, and Henderson smoothers", {
  Xt <- matrix(seq_len(20), ncol = 2)
  ytime <- fixture_time
  
  ns_result <- with_temp_bindings(
    list(
      fe_kp_quantile = function(...) c(1, 2, 3),
      fe_xt_spline = function(col, time, spline_fit, knot_points, ...) rep(if (spline_fit == "ns") 11 else 22, length(col))
    ),
    blce_smooth_Xt(Xt, ytime, smoother = "ns", knot_count = 3)
  )
  
  rm_result <- with_temp_bindings(
    list(
      fe_xt_rm = function(X, k, passes) matrix(33, nrow = nrow(X), ncol = ncol(X))
    ),
    blce_smooth_Xt(Xt, ytime, smoother = "rm", k_count = 7, k_pass = 2)
  )
  
  hp_result <- with_temp_bindings(
    list(
      fe_xt_hp = function(X, lambda) matrix(lambda, nrow = nrow(X), ncol = ncol(X))
    ),
    blce_smooth_Xt(Xt, ytime, smoother = "hp", hp_lambda = 1600)
  )
  
  henderson_result <- with_temp_bindings(
    list(
      fe_xt_henderson = function(X, k) matrix(k, nrow = nrow(X), ncol = ncol(X))
    ),
    blce_smooth_Xt(Xt, ytime, smoother = "henderson", henderson_k = 13)
  )
  
  expect_true(all(ns_result == 11))
  expect_true(all(rm_result == 33))
  expect_true(all(hp_result == 1600))
  expect_true(all(henderson_result == 13))
})

test_that("mc_fit_end assembles the expected top-level output structure", {
  yields <- make_pair_yields()
  fake_time <- yields$usa$time
  fake_tY <- data.frame(
    usa.01M = seq_along(fake_time),
    usa.03M = seq_along(fake_time) + 1,
    gbr.01M = seq_along(fake_time) + 2,
    gbr.03M = seq_along(fake_time) + 3
  )

  fake_bl_yields <- list(
    phi = list(usa = list(phi = matrix(1, nrow = 3, ncol = 2)), gbr = list(phi = matrix(1, nrow = 3, ncol = 2))),
    betas = list(usa = list(NSbetas = matrix(1, nrow = length(fake_time), ncol = 3)), gbr = list(NSbetas = matrix(1, nrow = length(fake_time), ncol = 3))),
    nsfit = list(usa = list(NSfit = matrix(1, nrow = length(fake_time), ncol = 2), W = matrix(0, nrow = length(fake_time), ncol = 2)),
                 gbr = list(NSfit = matrix(1, nrow = length(fake_time), ncol = 2), W = matrix(0, nrow = length(fake_time), ncol = 2))),
    H = list(usa = diag(2), gbr = diag(2))
  )

  result <- with_temp_bindings(
    list(
      vecY = function(...) list(tY = fake_tY, time = fake_time),
      bln_build_P = function(...) diag(4),
      PY_full = function(...) list(py = fake_tY),
      lambda_grid_search = function(...) 0.1,
      bl_single_curve = function(...) list(),
      transpose = function(...) fake_bl_yields,
      bln_phi_hat = function(...) matrix(1, nrow = 4, ncol = 3),
      blcc_build_Wj = function(...) list(X01M = matrix(1, nrow = length(fake_time), ncol = 2), X03M = matrix(1, nrow = length(fake_time), ncol = 2)),
      bdiag = function(x) diag(4),
      blcc_cointegration = function(...) list(cc_cspread = list(gbr = matrix(1, nrow = length(fake_time), ncol = 1)),
                                             cc_ecm = list(gbr = list(Fmatrix = matrix(1, nrow = 1, ncol = 1))),
                                             cc_beta = list()),
      blcc_build_Xt = function(...) matrix(1, nrow = length(fake_time), ncol = 1),
      blcc_fe = function(...) list(X01M = list(S0 = diag(2), B = matrix(1, nrow = 2, ncol = 1))),
      tSigma_optim = function(...) list(tSigma = rep(list(diag(4)), length(fake_time)), jSigma = list()),
      sqrt_inv_build = function(...) list(sqrt = rep(list(diag(4)), length(fake_time)), inverse = rep(list(diag(4)), length(fake_time))),
      bln_YB = function(...) fake_tY,
      bln_phi_BT = function(...) rep(list(matrix(1, nrow = 4, ncol = 3)), length(fake_time)),
      bln_Hb = function(...) diag(4),
      bln_phi_CT = function(...) rep(list(matrix(1, nrow = 4, ncol = 3)), length(fake_time)),
      cf_curve_cut = function(...) list(usa = matrix(1, nrow = length(fake_time), ncol = 2), gbr = matrix(1, nrow = length(fake_time), ncol = 2)),
      cf_phi_cut = function(...) list(usa = rep(list(matrix(1, nrow = 2, ncol = 3)), length(fake_time)),
                                      gbr = rep(list(matrix(1, nrow = 2, ncol = 3)), length(fake_time))),
      bln_normalized_fit = function(...) list(n_yield = list(usa = matrix(1, nrow = length(fake_time), ncol = 2),
                                                             gbr = matrix(1, nrow = length(fake_time), ncol = 2)),
                                              n_beta = list()),
      blce_rescale_yield = function(...) data.frame(time = fake_time, value = seq_along(fake_time))
    ),
    mc_fit_end(yields, lambdas = c(0.1, 0.2), cutoffs = c(1, 1), reference = "usa", curves = fixture_curves, mats = c(1, 3))
  )

  expect_equal(names(result), c("tenor", "curve", "sigma_JT", "curveECM", "BS0", "Xt"))
  expect_s3_class(result$tenor, "data.frame")
  expect_equal(nrow(result$Xt), length(fake_time))
})

test_that("mc_fit_exo assembles the expected top-level output structure", {
  yields <- make_pair_yields()
  fake_time <- yields$usa$time
  fake_tY <- data.frame(
    usa.01M = seq_along(fake_time),
    usa.03M = seq_along(fake_time) + 1,
    gbr.01M = seq_along(fake_time) + 2,
    gbr.03M = seq_along(fake_time) + 3
  )
  Xt <- matrix(1, nrow = length(fake_time), ncol = 2)
  fake_bl_yields <- list(
    phi = list(usa = list(phi = matrix(1, nrow = 3, ncol = 2)), gbr = list(phi = matrix(1, nrow = 3, ncol = 2))),
    betas = list(usa = list(NSbetas = matrix(1, nrow = length(fake_time), ncol = 3)), gbr = list(NSbetas = matrix(1, nrow = length(fake_time), ncol = 3))),
    nsfit = list(usa = list(NSfit = matrix(1, nrow = length(fake_time), ncol = 2), W = matrix(0, nrow = length(fake_time), ncol = 2)),
                 gbr = list(NSfit = matrix(1, nrow = length(fake_time), ncol = 2), W = matrix(0, nrow = length(fake_time), ncol = 2))),
    H = list(usa = diag(2), gbr = diag(2))
  )

  result <- with_temp_bindings(
    list(
      vecY = function(...) list(tY = fake_tY, time = fake_time),
      bln_build_P = function(...) diag(4),
      PY_full = function(...) list(py = fake_tY),
      lambda_grid_search = function(...) 0.1,
      bl_single_curve = function(...) list(),
      transpose = function(...) fake_bl_yields,
      bln_phi_hat = function(...) matrix(1, nrow = 4, ncol = 3),
      blcc_build_Wj = function(...) list(X01M = matrix(1, nrow = length(fake_time), ncol = 2), X03M = matrix(1, nrow = length(fake_time), ncol = 2)),
      bdiag = function(x) diag(4),
      blcc_fe = function(...) list(X01M = list(S0 = diag(2), B = matrix(1, nrow = 2, ncol = 2))),
      tSigma_optim = function(...) list(tSigma = rep(list(diag(4)), length(fake_time)), jSigma = list()),
      sqrt_inv_build = function(...) list(sqrt = rep(list(diag(4)), length(fake_time)), inverse = rep(list(diag(4)), length(fake_time))),
      bln_YB = function(...) fake_tY,
      bln_phi_BT = function(...) rep(list(matrix(1, nrow = 4, ncol = 3)), length(fake_time)),
      bln_Hb = function(...) diag(4),
      bln_phi_CT = function(...) rep(list(matrix(1, nrow = 4, ncol = 3)), length(fake_time)),
      cf_curve_cut = function(...) list(usa = matrix(1, nrow = length(fake_time), ncol = 2), gbr = matrix(1, nrow = length(fake_time), ncol = 2)),
      cf_phi_cut = function(...) list(usa = rep(list(matrix(1, nrow = 2, ncol = 3)), length(fake_time)),
                                      gbr = rep(list(matrix(1, nrow = 2, ncol = 3)), length(fake_time))),
      bln_normalized_fit = function(...) list(n_yield = list(usa = matrix(1, nrow = length(fake_time), ncol = 2),
                                                             gbr = matrix(1, nrow = length(fake_time), ncol = 2)),
                                              n_beta = list()),
      blce_rescale_yield = function(...) data.frame(time = fake_time, value = seq_along(fake_time))
    ),
    mc_fit_exo(yields, lambdas = c(0.1, 0.2), cutoffs = c(1, 1), reference = "usa", Xt = Xt, curves = fixture_curves, mats = c(1, 3))
  )

  expect_equal(names(result), c("tenor", "curve", "sigma_JT", "BS0", "Xt"))
  expect_s3_class(result$curve, "data.frame")
  expect_equal(ncol(result$Xt), 2)
})

test_that("mc_fit_end forwards S0 diagonal shrinkage to covariance regression", {
  yields <- make_pair_yields()
  fake_time <- yields$usa$time
  fake_tY <- data.frame(
    usa.01M = seq_along(fake_time),
    usa.03M = seq_along(fake_time) + 1,
    gbr.01M = seq_along(fake_time) + 2,
    gbr.03M = seq_along(fake_time) + 3
  )
  captured_shrink <- NULL
  fake_bl_yields <- list(
    phi = list(usa = list(phi = matrix(1, nrow = 3, ncol = 2)), gbr = list(phi = matrix(1, nrow = 3, ncol = 2))),
    betas = list(usa = list(NSbetas = matrix(1, nrow = length(fake_time), ncol = 3)), gbr = list(NSbetas = matrix(1, nrow = length(fake_time), ncol = 3))),
    nsfit = list(usa = list(NSfit = matrix(1, nrow = length(fake_time), ncol = 2), W = matrix(0, nrow = length(fake_time), ncol = 2)),
                 gbr = list(NSfit = matrix(1, nrow = length(fake_time), ncol = 2), W = matrix(0, nrow = length(fake_time), ncol = 2))),
    H = list(usa = diag(2), gbr = diag(2))
  )

  with_temp_bindings(
    list(
      vecY = function(...) list(tY = fake_tY, time = fake_time),
      bln_build_P = function(...) diag(4),
      PY_full = function(...) list(py = fake_tY),
      lambda_grid_search = function(...) 0.1,
      bl_single_curve = function(...) list(),
      transpose = function(...) fake_bl_yields,
      bln_phi_hat = function(...) matrix(1, nrow = 4, ncol = 3),
      blcc_build_Wj = function(...) list(X01M = matrix(1, nrow = length(fake_time), ncol = 2)),
      bdiag = function(x) diag(4),
      blcc_cointegration = function(...) list(cc_cspread = list(gbr = matrix(1, nrow = length(fake_time), ncol = 1)),
                                             cc_ecm = list(gbr = list(Fmatrix = matrix(1, nrow = 1, ncol = 1))),
                                             cc_beta = list()),
      blcc_build_Xt = function(...) matrix(1, nrow = length(fake_time), ncol = 1),
      blcc_fe = function(..., S0_shrink_diag) {
        captured_shrink <<- S0_shrink_diag
        list(X01M = list(S0 = diag(2), B = matrix(1, nrow = 2, ncol = 1)))
      },
      tSigma_optim = function(...) list(tSigma = rep(list(diag(2)), length(fake_time)), jSigma = list()),
      sqrt_inv_build = function(...) list(sqrt = rep(list(diag(4)), length(fake_time)), inverse = rep(list(diag(4)), length(fake_time))),
      bln_YB = function(...) fake_tY,
      bln_phi_BT = function(...) rep(list(matrix(1, nrow = 4, ncol = 3)), length(fake_time)),
      bln_Hb = function(...) diag(4),
      bln_phi_CT = function(...) rep(list(matrix(1, nrow = 4, ncol = 3)), length(fake_time)),
      cf_curve_cut = function(...) list(usa = matrix(1, nrow = length(fake_time), ncol = 2), gbr = matrix(1, nrow = length(fake_time), ncol = 2)),
      cf_phi_cut = function(...) list(usa = rep(list(matrix(1, nrow = 2, ncol = 3)), length(fake_time)),
                                      gbr = rep(list(matrix(1, nrow = 2, ncol = 3)), length(fake_time))),
      bln_normalized_fit = function(...) list(n_yield = list(usa = matrix(1, nrow = length(fake_time), ncol = 2),
                                                             gbr = matrix(1, nrow = length(fake_time), ncol = 2)),
                                              n_beta = list()),
      blce_rescale_yield = function(...) data.frame(time = fake_time, value = seq_along(fake_time))
    ),
    mc_fit_end(yields, lambdas = c(0.1, 0.2), cutoffs = c(1, 1), reference = "usa",
               curves = fixture_curves, mats = c(1, 3), CR_S0_shrink_diag = 0.6)
  )

  expect_equal(captured_shrink, 0.6)
})
