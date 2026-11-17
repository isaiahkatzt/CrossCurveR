test_that("dr_row_model returns original values when the row is already complete", {
  complete_row <- c(0, 1, 1.1, 1.2, 1.3)
  fitted <- dr_row_model(complete_row, method = "nelson", fmat = c(1, 3, 6, 12), complete = FALSE)

  expect_equal(as.numeric(fitted), as.numeric(complete_row[-1]), tolerance = tol)
})

test_that("dr_row_model fills missing values while preserving observed points when complete is FALSE", {
  mats <- c(1, 3, 6, 12, 24)
  curve <- make_ns_curve_df(time = as.Date("2020-01-01"), mats = mats)
  yields <- as.numeric(curve[1, -1])
  yields[2] <- NA
  fitted <- as.numeric(dr_row_model(c(0, yields), fmat = mats, complete = FALSE))
  expect_equal(fitted[-2], yields[-2])
  expect_gt(fitted[2], 0)
})

test_that("parametric reconstruction preserves rows with too few observations", {
  expect_warning(
    sparse <- dr_row_model(c(0, 1, NA, NA, NA), fmat = c(1, 3, 6, 12)),
    "Insufficient data"
  )
  expect_equal(unname(sparse), c(1, NA, NA, NA))
})
