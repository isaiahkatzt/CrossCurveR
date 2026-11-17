test_that("blmc_VAR orders the reference factors before the marginal factors", {
  fits <- list(betas = list(
    gbr = list(NSbetas = matrix(1:9, 3)),
    usa = list(NSbetas = matrix(10:18, 3))
  ))
  expected <- matrix(c(10:18, 1:9), 3, dimnames = list(NULL,
    c("usa.L", "usa.S", "usa.C", "gbr.L", "gbr.S", "gbr.C")))
  expect_equal(blmc_VAR(fits, "usa", fixture_curves), list(gbr = expected))
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

test_that("ECM errors preserve dependency details and available sample and lag context", {
  betas <- matrix(seq_len(24), ncol = 6)
  for (dependency in c("VARselect", "VECM", "rank.test")) {
    bindings <- list(VARselect = function(...) list(selection = c("SC(n)" = 3)),
                     VECM = function(...) list())
    bindings[[dependency]] <- function(...) stop("original dependency error")
    expected <- paste0("ECM estimation: failed.*observations=4, variables=6",
                       if (dependency != "VARselect") ", VECM lag=2" else "",
                       ".*original dependency error")
    with_temp_bindings(bindings, expect_error(blmc_ECM(betas), expected))
  }
})

test_that("spread normalization identifies zero and undefined standard deviations", {
  betas <- cbind(varying = 1:4, constant = rep(2, 4))
  expect_error(blmc_cspread(betas, diag(2), curve_name = "gbr"),
               "Spread normalization: zero or nonfinite standard deviation for gbr.r2")
  expect_error(blmc_cspread(betas[1, , drop = FALSE], diag(2), curve_name = "gbr"),
               "Spread normalization: zero or nonfinite standard deviation for gbr.r1, gbr.r2")
  raw <- blmc_cspread(betas, diag(2), curve_name = "gbr", normalize = FALSE)
  expect_identical(unname(raw), unname(betas))
  expect_error(blcc_build_Xt(list(gbr = NULL), list(gbr = list(rank = 0))),
               "Cointegration feature construction: no usable relationships")
})

test_that("blmc_cspread standardizes and names each spread", {
  cc_betas <- matrix(seq_len(18), ncol = 6)
  colnames(cc_betas) <- c("usa.L", "usa.S", "usa.C", "gbr.L", "gbr.S", "gbr.C")
  G <- Matrix::diag(6)[, 1:2, drop = FALSE]
  colnames(G) <- c("r1", "r2")

  spread <- blmc_cspread(cc_betas, G, curve_name = "gbr")
  standardized <- as.matrix(scale(cc_betas[, 1:2]))
  expect_equal(spread, matrix(standardized, ncol = 2,
                             dimnames = list(NULL, c("gbr.r1", "gbr.r2"))))
})

test_that("feature truncation selects the strongest adjustment for each curve", {
  spreads <- list(gbr = matrix(1:9, ncol = 3), jpn = matrix(10:18, ncol = 3))
  ecms <- list(
    gbr = list(Fmatrix = matrix(c(3, 4, 0, 0, 0, 1), nrow = 2), rank = 3),
    jpn = list(Fmatrix = matrix(c(0, 1, 0, 5, 2, 1), nrow = 2), rank = 3)
  )
  expect_equal(blcc_build_Xt(spreads, ecms, trunc = TRUE),
               cbind(spreads$gbr[, 1], spreads$jpn[, 2]))
})

test_that("covariance builders preserve values, ordering, and singleton outputs", {
  object <- list(B = diag(2), S0 = diag(c(2, 3)))
  X <- rbind(c(1, 2), c(3, 4))
  expected <- array(c(3, 2, 2, 7, 11, 12, 12, 19), c(2, 2, 2))
  expect_equal(cr_build(object, X), expected)
  expect_equal(jSigma_optim(object, X), list(expected[, , 1], expected[, , 2]))
  for (as_array in c(TRUE, FALSE)) for (drop in c(TRUE, FALSE)) {
    singleton <- cr_build(object, X[1, , drop = FALSE], return_array = as_array, drop_single = drop)
    if (!as_array && !drop) singleton <- singleton[[1]]
    expect_equal(singleton, expected[, , 1, drop = FALSE])
  }
  simulated <- cr_synthetic(N = 2, B = object$B, S0 = object$S0, X = X, seed = 1)
  expect_equal(simulated$Sigma_true, expected)

  assembled <- tSigma_optim(list(X02Y = object,
    X01Y = list(B = matrix(0, 2, 2), S0 = diag(c(5, 6)))), X)
  expected_block <- diag(c(3, 7, 5, 6))
  expected_block[1, 2] <- expected_block[2, 1] <- 2
  expect_named(assembled$jSigma, c("X02Y", "X01Y"))
  expect_s4_class(assembled$tSigma[[1]], "sparseMatrix")
  expect_equal(as.matrix(assembled$tSigma[[1]]), expected_block)
})

test_that("running means use shorter windows at the boundaries", {
  input <- cbind(a = c(0, 0, 9, 0, 0), b = rep(2, 5))
  expected <- cbind(a = c(0, 3, 3, 3, 0), b = rep(2, 5))
  expect_equal(fe_xt_rm(input, k = 3, passes = 1), expected)
})

test_that("spline smoothing uses numeric observation coordinates without a date origin", {
  for (time in list(-10L:9L, -5 + cumsum(seq_len(20)^1.3) / 10)) {
    xt <- sin(seq_along(time))
    knots <- stats::quantile(time, probs = seq(0, 1, length.out = 5))[2:4]
    expect_identical(fe_kp_quantile(time, knot_count = 3), knots)
    for (kind in c("ns", "bs")) {
      basis <- if (kind == "ns") splines::ns(time, knots = knots) else
        splines::bs(time, knots = knots, degree = 3)
      expected <- stats::predict(stats::lm(xt ~ basis))
      expect_equal(fe_xt_spline(time, xt, knots, spline_fit = kind), expected)
      expect_equal(unname(blce_smooth_Xt(cbind(x = xt), time, kind, knot_count = 3)),
                   matrix(unname(expected), ncol = 1))
    }
  }
})

test_that("HP smoother returns same-shaped smoothed features", {
  xt <- sin(seq(0, 4 * pi, length.out = length(fixture_time))) + rep(c(-0.2, 0.2), length.out = length(fixture_time))
  X <- cbind(a = xt, b = xt + seq_along(xt) / 10)

  hp_smooth <- fe_xt_hp(X, lambda = 1600)
  expect_identical(dimnames(hp_smooth), dimnames(X))
  expect_equal(fe_xt_hp(xt, lambda = 1600), as.numeric(hp_smooth[, "a"]))
  expect_lt(sum(diff(hp_smooth[, "a"], differences = 2)^2), sum(diff(X[, "a"], differences = 2)^2))
  expect_error(fe_xt_hp(X, lambda = -1), "lambda must be a non-negative numeric scalar")
})

test_that("Henderson smoother returns same-shaped smoothed features", {
  xt <- sin(seq(0, 4 * pi, length.out = length(fixture_time))) + rep(c(-0.2, 0.2), length.out = length(fixture_time))
  X <- cbind(a = xt, b = xt + seq_along(xt) / 10)

  henderson_smooth <- fe_xt_henderson(X, k = 13)
  expect_identical(dimnames(henderson_smooth), dimnames(X))
  expect_equal(fe_xt_henderson(xt, k = 13), as.numeric(henderson_smooth[, "a"]))
  expect_lt(sum(diff(henderson_smooth[, "a"], differences = 2)^2), sum(diff(X[, "a"], differences = 2)^2))
  expect_error(fe_xt_henderson(X, k = 12), "k must be a positive odd integer")
  expect_error(fe_xt_henderson(X, k = 0), "k must be a positive odd integer")
  expect_error(fe_xt_henderson(X, k = c(9, 13)), "k must be a positive odd integer")
})
