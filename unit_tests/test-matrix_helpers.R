test_that("symmPSD_sqrt reconstructs positive semidefinite matrices", {
  Sigma <- matrix(c(2, 1, 1, 2), nrow = 2)
  sqrt_sigma <- symmPSD_sqrt(Sigma)

  expect_equal(sqrt_sigma %*% sqrt_sigma, Sigma, tolerance = 1e-5)
  expect_equal(sqrt_sigma, t(sqrt_sigma), tolerance = tol)
})

test_that("sqrt_inv_build returns square roots and inverses for each matrix", {
  sigma_list <- list(diag(2), matrix(c(2, 0.5, 0.5, 1), nrow = 2))
  built <- sqrt_inv_build(sigma_list)

  expect_equal(length(built$sqrt), 2)
  expect_equal(length(built$inverse), 2)
  expect_equal(built$sqrt[[1]] %*% built$inverse[[1]], diag(2), tolerance = 1e-5)
  expect_equal(built$sqrt[[2]] %*% built$inverse[[2]], diag(2), tolerance = 1e-5)
})
