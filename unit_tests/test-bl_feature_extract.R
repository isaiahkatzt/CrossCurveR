test_that("blmc_phi stacks curve phi matrices from its input list", {
  blsc_list <- make_mock_blsc_list()
  phi_hat <- blmc_phi(blsc_list)
  expected_phi <- do.call(rbind, lapply(blsc_list$phi, function(curve) t(curve$phi)))

  expect_equal(dim(phi_hat), c(length(fixture_curves) * length(fixture_mats), 3))
  expect_equal(unname(phi_hat), unname(expected_phi), tolerance = tol)
})

test_that("blmc_VAR stacks reference and marginal beta matrices with labeled columns", {
  blsc_list <- make_mock_blsc_list()
  cross_beta <- blmc_VAR(blsc_list, reference = "usa", curves = fixture_curves)

  expect_equal(names(cross_beta), "gbr")
  expect_equal(ncol(cross_beta$gbr), 6)
  expect_equal(
    colnames(cross_beta$gbr),
    c("usa.L", "usa.S", "usa.C", "gbr.L", "gbr.S", "gbr.C")
  )
})

test_that("blmc_ECM returns restricted matrices when rank is positive", {
  cc_betas <- matrix(seq_len(60), ncol = 6)

  result <- with_temp_bindings(
    list(
      VARselect = function(...) list(selection = c("SC(n)" = 2)),
      VECM = function(...) structure(list(), class = "mock_vecm"),
      rank.test = function(...) list(r = 1),
      coefA = function(...) matrix(c(1, 0, 0, 1), nrow = 2),
      coefB = function(...) matrix(c(1, -1, 0, 0, 0, 0), nrow = 6)
    ),
    blmc_ECM(cc_betas)
  )

  expect_equal(result$rank, 1)
  expect_true(is.matrix(result$Fmatrix))
  expect_true(is.matrix(result$Gmatrix))
})

test_that("blmc_ECM handles the zero-rank case", {
  cc_betas <- matrix(seq_len(60), ncol = 6)

  result <- with_temp_bindings(
    list(
      VARselect = function(...) list(selection = c("SC(n)" = 2)),
      VECM = function(...) structure(list(), class = "mock_vecm"),
      rank.test = function(...) list(r = 0)
    ),
    blmc_ECM(cc_betas)
  )

  expect_equal(result, list(Gmatrix = 0, Fmatrix = 0, rank = 0))
})

test_that("spread and reversion helpers normalize and rank feature blocks", {
  cc_betas <- matrix(1:12, ncol = 3)
  G <- diag(3)
  spread <- blmc_cspread(cc_betas, G, normalize = TRUE)

  expect_equal(unname(round(colMeans(spread), 8)), c(0, 0, 0))
  expect_equal(unname(round(apply(spread, 2, sd), 8)), c(1, 1, 1))
  expect_equal(blmc_Fselect(matrix(c(3, 4, 0, 1), nrow = 2)), 1)
})

test_that("blmc_cspread prefixes spread columns with the non-reference curve name", {
  cc_betas <- matrix(seq_len(18), ncol = 6)
  colnames(cc_betas) <- c("usa.L", "usa.S", "usa.C", "gbr.L", "gbr.S", "gbr.C")
  G <- diag(6)[, 1:2, drop = FALSE]
  colnames(G) <- c("r1", "r2")

  spread <- blmc_cspread(cc_betas, G, curve_name = "gbr", normalize = FALSE)

  expect_equal(colnames(spread), c("gbr.r1", "gbr.r2"))
})

test_that("blcc_cointegration, Xt builders, and Wj builders assemble expected structures", {
  blsc_list <- make_mock_blsc_list()

  cointegration <- with_temp_bindings(
    list(
      blmc_VAR = function(...) list(gbr = matrix(1:12, ncol = 3)),
      blmc_ECM = function(...) list(Fmatrix = matrix(c(3, 0, 0, 1), nrow = 2), Gmatrix = diag(3), rank = 1),
      blmc_cspread = function(..., curve_name) {
        out <- matrix(1:9, ncol = 3)
        colnames(out) <- paste0(curve_name, ".", paste0("r", 1:3))
        out
      }
    ),
    blcc_cointegration(blsc_list, curves = fixture_curves, reference = "usa")
  )

  expect_equal(names(cointegration), c("cc_cspread", "cc_ecm", "cc_beta"))
  expect_equal(names(cointegration$cc_cspread), "gbr")
  expect_equal(colnames(cointegration$cc_cspread$gbr), c("gbr.r1", "gbr.r2", "gbr.r3"))

  full_xt <- blcc_build_Xt(
    cc_cspread = list(gbr = matrix(1:9, ncol = 3), jpn = matrix(10:18, ncol = 3)),
    cc_ecm = list(
      gbr = list(Fmatrix = matrix(c(3, 4, 0, 0), nrow = 2), rank = 1),
      jpn = list(Fmatrix = matrix(c(0, 1, 0, 5), nrow = 2), rank = 1)
    ),
    trunc = TRUE
  )
  expect_equal(ncol(full_xt), 2)

  Wj <- blcc_build_Wj(blsc_list$nsfit, fixture_mats)
  expect_equal(length(Wj), length(fixture_mats))
  expect_true(all(grepl("^X", names(Wj))))
})

test_that("feature-extraction wrappers and sigma optimizers return structured outputs", {
  Xt <- matrix(seq_len(20), ncol = 2)
  Wt <- list(X01M = matrix(seq_len(20), ncol = 2), X03M = matrix(seq_len(20), ncol = 2))

  BS0 <- blcc_fe(
    Xt,
    Wt,
    fixture_mats,
    covreg = function(W, X, ...) list(S0 = diag(ncol(W)), B = matrix(1, nrow = ncol(W), ncol = ncol(X)))
  )

  expect_equal(names(BS0), names(Wt))

  sigma_j <- sigma_jt_optim(list(S0 = diag(2), B = matrix(c(1, 0, 0, 1), nrow = 2)), c(1, 2))
  sigma_path <- jSigma_optim(list(S0 = diag(2), B = matrix(c(1, 0, 0, 1), nrow = 2)), Xt)
  sigma_time <- tSigma_optim(
    list(
      a = list(S0 = diag(2), B = matrix(c(1, 0, 0, 1), nrow = 2)),
      b = list(S0 = diag(2), B = matrix(c(1, 0, 0, 1), nrow = 2))
    ),
    Xt
  )

  expect_equal(dim(sigma_j), c(2, 2))
  expect_equal(length(sigma_path), nrow(Xt))
  expect_equal(length(sigma_time$tSigma), nrow(Xt))
})

test_that("time smoothers build knot points, spline bases, and rolling averages", {
  knot_points <- fe_kp_quantile(fixture_time, knot_count = 3)
  ns_basis <- fe_spline_basis("ns")
  bs_basis <- fe_spline_basis("bs")
  xt <- seq_len(length(fixture_time))

  ns_smooth <- fe_xt_spline(fixture_time, xt, knot_points, spline_fit = "ns")
  bs_smooth <- fe_xt_spline(fixture_time, xt, knot_points, spline_fit = "bs")
  rm_smooth <- fe_xt_rm(cbind(a = xt, b = xt + 1), k = 5, passes = 2)

  expect_equal(length(knot_points), 3)
  expect_true(is.function(ns_basis))
  expect_true(is.function(bs_basis))
  expect_equal(length(ns_smooth), length(fixture_time))
  expect_equal(length(bs_smooth), length(fixture_time))
  expect_equal(dim(rm_smooth), c(length(fixture_time), 2))
})

test_that("HP smoother returns same-shaped smoothed features", {
  xt <- sin(seq(0, 4 * pi, length.out = length(fixture_time))) + rep(c(-0.2, 0.2), length.out = length(fixture_time))
  X <- cbind(a = xt, b = xt + seq_along(xt) / 10)
  
  hp_smooth <- fe_xt_hp(X, lambda = 1600)
  hp_vector <- fe_xt_hp(xt, lambda = 1600)
  
  expect_equal(dim(hp_smooth), dim(X))
  expect_equal(colnames(hp_smooth), colnames(X))
  expect_equal(length(hp_vector), length(xt))
  expect_lt(sum(diff(hp_smooth[, "a"], differences = 2)^2), sum(diff(X[, "a"], differences = 2)^2))
  expect_error(fe_xt_hp(X, lambda = -1), "lambda must be a non-negative numeric scalar")
})

test_that("Henderson smoother returns same-shaped smoothed features", {
  xt <- sin(seq(0, 4 * pi, length.out = length(fixture_time))) + rep(c(-0.2, 0.2), length.out = length(fixture_time))
  X <- cbind(a = xt, b = xt + seq_along(xt) / 10)
  
  henderson_smooth <- fe_xt_henderson(X, k = 13)
  henderson_vector <- fe_xt_henderson(xt, k = 13)
  
  expect_equal(dim(henderson_smooth), dim(X))
  expect_equal(colnames(henderson_smooth), colnames(X))
  expect_equal(length(henderson_vector), length(xt))
  expect_lt(sum(diff(henderson_smooth[, "a"], differences = 2)^2), sum(diff(X[, "a"], differences = 2)^2))
  expect_error(fe_xt_henderson(X, k = 12), "k must be a positive odd integer")
  expect_error(fe_xt_henderson(X, k = 0), "k must be a positive odd integer")
  expect_error(fe_xt_henderson(X, k = c(9, 13)), "k must be a positive odd integer")
})
