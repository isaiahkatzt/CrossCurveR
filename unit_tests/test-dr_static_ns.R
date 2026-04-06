test_that("dr_row_model returns original values when the row is already complete", {
  complete_row <- c(0, 1, 1.1, 1.2, 1.3)
  fitted <- dr_row_model(complete_row, method = "nelson", fmat = c(1, 3, 6, 12), complete = FALSE)

  expect_equal(as.numeric(fitted), as.numeric(complete_row[-1]), tolerance = tol)
})

test_that("dr_row_model fills missing values while preserving observed points when complete is FALSE", {
  base_curve <- make_ns_curve_df(time = as.Date("2020-01-01"), mats = c(1, 3, 6, 12))
  row <- c(0, as.numeric(base_curve[1, -1]))
  row[3] <- NA

  fitted <- as.numeric(dr_row_model(row, method = "nelson", fmat = c(1, 3, 6, 12), complete = FALSE))

  expect_equal(fitted[c(1, 3, 4)], as.numeric(row[c(2, 4, 5)]), tolerance = 1e-4)
  expect_true(is.finite(fitted[2]))
})

test_that("dr_row_model supports complete fits and warns on too few observations", {
  base_curve <- make_ns_curve_df(time = as.Date("2020-01-01"), mats = c(1, 3, 6, 12))
  row <- c(0, as.numeric(base_curve[1, -1]))
  row[4] <- NA

  full_fit <- dr_row_model(row, method = "svensson", fmat = c(1, 3, 6, 12), complete = TRUE)
  expect_equal(length(full_fit), 4)
  expect_true(all(is.finite(full_fit)))

  expect_warning(
    sparse_fit <- dr_row_model(c(0, 1, NA, NA, NA), fmat = c(1, 3, 6, 12)),
    "Insufficient data"
  )
  expect_equal(unname(sparse_fit), c(1, NA, NA, NA))
})
