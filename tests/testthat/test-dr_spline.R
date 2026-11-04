test_that("dr_row_spline interpolates missing points on the requested grid", {
  row <- c("2020-01-01", 1, NA, 3)
  interpolated <- dr_row_spline(row, method = "natural", fmat = c(1, 2, 3))

  expect_equal(length(interpolated), 3)
  expect_equal(interpolated[c(1, 3)], c(1, 3), tolerance = tol)
  expect_true(is.finite(interpolated[2]))
})

test_that("dr_row_spline accepts hermite and warns on insufficient inputs", {
  hermite_fit <- dr_row_spline(c("2020-01-01", 1, NA, 4), method = "hermite", fmat = c(1, 2, 3))
  expect_equal(length(hermite_fit), 3)

  expect_warning(
    sparse_fit <- dr_row_spline(c("2020-01-01", 1, NA, NA), fmat = c(1, 2, 3)),
    "Insufficient data"
  )
  expect_equal(unname(sparse_fit), c("1", NA, NA))
})

test_that("dr_row_loglin performs log-linear interpolation and warns when needed", {
  row <- c("2020-01-01", 1, NA, 4)
  interpolated <- dr_row_loglin(row, fmat = c(1, 2, 3))

  expect_equal(length(interpolated), 3)
  expect_equal(interpolated[c(1, 3)], c(1, 4), tolerance = tol)
  expect_equal(interpolated[2], 2, tolerance = 1e-4)

  expect_warning(
    sparse_fit <- dr_row_loglin(c("2020-01-01", 1, NA, NA), fmat = c(1, 2, 3)),
    "Insufficient data"
  )
  expect_equal(unname(sparse_fit), c("1", NA, NA))
})
