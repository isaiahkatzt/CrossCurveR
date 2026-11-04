test_that("zero_mean_covreg returns correctly shaped covariance-regression outputs", {
  set.seed(123)
  W <- cbind(seq(1, 10), seq(1, 10)^2)
  X <- cbind(1, seq(-1, 1, length.out = 10))

  fit <- zero_mean_covreg(W, X, init = "static", max_iter = 5)

  expect_true(is.matrix(fit$S0))
  expect_true(is.matrix(fit$B))
  expect_equal(dim(fit$S0), c(ncol(W), ncol(W)))
  expect_equal(dim(fit$B), c(ncol(W), ncol(X)))
})

test_that("zero_mean_covreg defaults to adaptive initialization", {
  W <- cbind(seq(1, 10), seq(1, 10)^2)
  X <- cbind(1, seq(-1, 1, length.out = 10))

  set.seed(123)
  default_fit <- zero_mean_covreg(W, X, max_iter = 5)
  set.seed(123)
  adaptive_fit <- zero_mean_covreg(W, X, init = "adaptive", max_iter = 5)

  expect_identical(default_fit, adaptive_fit)
})

test_that("zero_mean_covreg supports adaptive initialization and early termination mode", {
  set.seed(123)
  W <- cbind(seq(1, 12), seq(1, 12)^2)
  X <- cbind(1, seq(-1, 1, length.out = 12))

  fit <- zero_mean_covreg(W, X, init = "adaptive", max_iter = 20, tol = 1e8, term = TRUE)

  expect_true(is.list(fit$b_est))
  expect_true(is.list(fit$s_est))
  expect_lte(fit$iter, 20)
})

test_that("covariance regression shrinkage can diagonalize S0", {
  set.seed(123)
  W <- cbind(seq(1, 10), seq(1, 10) + c(rep(0, 5), rep(3, 5)))
  X <- cbind(1, seq(-1, 1, length.out = 10))

  fit_dense <- zero_mean_covreg_em(W, X, init = "static", max_iter = 3, S0_shrink_diag = 0)
  fit_diag <- zero_mean_covreg_em(W, X, init = "static", max_iter = 3, S0_shrink_diag = 1)

  expect_gt(abs(fit_dense$S0[1, 2]), 0)
  expect_equal(fit_diag$S0[1, 2], 0, tolerance = 1e-10)
  expect_equal(fit_diag$S0_shrink_diag, 1)
})

test_that("regularized estimators use the stored covariance in E-steps and likelihoods", {
  W <- matrix(c(-1, .2, .5, -.4, 1.2, .8, -.3, 1.1, .7, -.6), ncol = 2, byrow = TRUE)
  X <- cbind(1, seq(-1, 1, length.out = nrow(W)))
  B0 <- matrix(c(.3, .2, -.1, .25), nrow = 2)
  S0 <- diag(c(.4, .6))
  ridge <- .2
  initial_covariance <- S0 + diag(ridge, 2)
  u <- X %*% t(B0)
  precision <- solve(initial_covariance)
  expected_vz <- 1 / (1 + rowSums((u %*% precision) * u))
  expected_mz <- expected_vz * rowSums((u %*% precision) * W)

  for (estimator in list(zero_mean_covreg_em, covreg_comp_estimator)) {
    fit <- estimator(W, X, B = B0, S0 = S0, ridge_S0 = ridge, max_iter = 1)
    expect_equal(fit$vz, expected_vz, tolerance = 1e-12)
    expect_equal(fit$mz, expected_mz, tolerance = 1e-12)

    direct_loglik <- sum(vapply(seq_len(nrow(W)), function(i) {
      covariance <- fit$S0 + tcrossprod(fit$B %*% X[i, ])
      -.5 * (ncol(W) * log(2 * pi) +
        as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
        as.numeric(crossprod(W[i, ], solve(covariance, W[i, ]))))
    }, numeric(1)))
    expect_equal(fit$loglik, direct_loglik, tolerance = 1e-12)
  }
})
