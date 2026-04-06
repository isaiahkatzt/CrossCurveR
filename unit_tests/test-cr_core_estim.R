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

test_that("zero_mean_covreg supports adaptive initialization and early termination mode", {
  set.seed(123)
  W <- cbind(seq(1, 12), seq(1, 12)^2)
  X <- cbind(1, seq(-1, 1, length.out = 12))

  fit <- zero_mean_covreg(W, X, init = "adaptive", max_iter = 20, tol = 1e8, term = TRUE)

  expect_true(is.list(fit$b_est))
  expect_true(is.list(fit$s_est))
  expect_lte(fit$iter, 20)
})

test_that("dynamic_mean_covreg is documented as pending implementation", {
  skip("dynamic_mean_covreg is not implemented yet")
})
