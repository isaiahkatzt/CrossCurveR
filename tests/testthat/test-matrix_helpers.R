test_that("loading inverse diagnostics retain the solver's conditioning policy", {
  gram <- diag(c(1, 1e-14))
  expect_identical(.bl_loading_inverse(gram, "Test loading regression"), Matrix::solve(gram))
  expect_error(.bl_loading_inverse(matrix(0, 1, 1), "Test loading regression"),
               "Test loading regression.*could not be inverted.*rcond=")
  expect_error(.bl_loading_inverse(Matrix::Diagonal(x = c(1, 1e-320)), "Test loading regression"),
               "Test loading regression.*inverse contains nonfinite values")
})

test_that("symmPSD_sqrt reconstructs positive semidefinite matrices", {
  Sigma <- matrix(c(2, 1, 1, 2), nrow = 2)
  sqrt_sigma <- symmPSD_sqrt(Sigma)

  expect_equal(sqrt_sigma %*% sqrt_sigma, Sigma, tolerance = 1e-5)
  expect_equal(sqrt_sigma, Matrix::t(sqrt_sigma), tolerance = tol)
})

test_that("sqrt_inv_build returns square roots and inverses for each matrix", {
  sigma_list <- list(Matrix::diag(2), matrix(c(2, 0.5, 0.5, 1), nrow = 2))
  built <- sqrt_inv_build(sigma_list)

  expect_equal(length(built$sqrt), 2)
  expect_equal(length(built$inverse), 2)
  expect_equal(built$sqrt[[1]] %*% built$inverse[[1]], Matrix::diag(2), tolerance = 1e-5)
  expect_equal(built$sqrt[[2]] %*% built$inverse[[2]], Matrix::diag(2), tolerance = 1e-5)
})

test_that("square roots handle scalar, singular, and rounding-boundary covariances", {
  expect_equal(symmPSD_sqrt(matrix(4, 1, 1)), matrix(2, 1, 1))
  expect_equal(symmPSD_sqrt(matrix(0, 1, 1)), matrix(0, 1, 1))
  expect_equal(symmPSD_sqrt(diag(c(1, 0))), diag(c(1, 0)))
  for (scale in c(1, 1e-12)) {
    expect_equal(symmPSD_sqrt(scale * diag(c(1, -1e-9))), sqrt(scale) * diag(c(1, 0)))
    expect_error(symmPSD_sqrt(scale * diag(c(1, -1e-5))), "minimum eigenvalue .*; tolerance")
  }
  expect_equal(symmPSD_sqrt(diag(c(1, 1e-12)))[2, 2], 1e-6, tolerance = 1e-12)
})

test_that("square roots reconcile negligible asymmetry and reject material asymmetry", {
  A <- matrix(c(2, .5, .5 + 1e-10, 1), 2)
  root <- symmPSD_sqrt(A)
  expect_equal(root %*% root, (A + t(A)) / 2, tolerance = 1e-12)
  expect_error(symmPSD_sqrt(matrix(c(2, 0, 1, 2), 2)), "materially asymmetric")
})

test_that("whitening preserves invertible inputs and identifies unusable covariances", {
  built <- sqrt_inv_build(list(scalar = matrix(4, 1, 1), small = diag(c(1, 1e-12))))
  expect_named(built$inverse, c("scalar", "small"))
  expect_equal(built$inverse$scalar, matrix(.5, 1, 1))
  expect_equal(built$inverse$small[2, 2], 1e6)
  for (small in c(0, 1e-20)) {
    expect_error(sqrt_inv_build(list(diag(2), diag(c(1, small)))), "element 2 .*singular")
  }
})
