test_that("dr_normalize_timestamp handles supported timestamp types", {
  posix_input <- as.POSIXct("2020-01-01 12:30:00", tz = "UTC")
  date_input <- as.Date("2020-01-02")
  numeric_input <- 86400

  expect_s3_class(dr_normalize_timestamp(posix_input), "POSIXct")
  expect_equal(as.Date(dr_normalize_timestamp(date_input)), date_input)
  expect_equal(as.Date(dr_normalize_timestamp(numeric_input)), as.Date("1970-01-02"))
  expect_equal(as.Date(dr_normalize_timestamp("2020-01-03")), as.Date("2020-01-03"))
  expect_equal(as.Date(dr_normalize_timestamp("04/01/2020")), as.Date("2020-01-04"))
})

test_that("dr_normalize_timestamp rejects unsupported or invalid timestamps", {
  expect_error(
    dr_normalize_timestamp("not-a-date"),
    "could not be parsed"
  )
  expect_error(
    dr_normalize_timestamp(list("2020-01-01")),
    "Unsupported timestamp type"
  )
})

test_that("validate_input_data checks maturity counts and returns normalized time", {
  input_data <- data.frame(
    date = c("2020-01-01", "2020-01-02"),
    x1 = c(1, 2),
    x2 = c(3, 4)
  )

  expect_s3_class(.validate_input_data(input_data, mats = c(1, 12)), "POSIXct")
  expect_error(
    .validate_input_data(input_data, mats = c(1)),
    "columns must match input maturities"
  )
})

test_that("dr_yield_format renames columns, filters dates, and converts NaN to NA", {
  input_data <- data.frame(
    date = c("2020-01-01", "2020-01-02", "2020-01-03"),
    short = c(1, NaN, 3),
    long = c(2, 4, 6)
  )

  formatted <- dr_yield_format(
    data = input_data,
    mats = c(1, 12),
    start_date = "2020-01-02",
    end_date = "2020-01-03"
  )

  expect_equal(colnames(formatted), c("time", "01M", "01Y"))
  expect_equal(formatted$time, c("2020-01-02", "2020-01-03"))
  expect_true(is.na(formatted$`01M`[1]))
})

test_that("dr_daily interpolates sparse rows and can augment to a reference calendar", {
  curve_data <- data.frame(
    time = as.Date("2020-01-01") + 0:2,
    `01M` = c(1, NA, 3),
    `03M` = c(2, NA, 4),
    `06M` = c(3, 10, 5)
  )

  interpolated <- dr_daily(curve_data, threshold = 2)
  expect_equal(interpolated[[2]][2], 2, tolerance = tol)
  expect_equal(interpolated[[3]][2], 3, tolerance = tol)
  expect_equal(interpolated[[4]][2], 10, tolerance = tol)

  reference <- as.Date("2019-12-31") + 0:4
  augmented <- dr_daily(curve_data, reference = reference, threshold = 2)
  expect_equal(nrow(augmented), length(reference))
  expect_equal(as.Date(augmented$time), reference)
})

test_that("dr_build_curve returns ordered tenor grids for spline, log, and parametric builders", {
  base_curve <- make_ns_curve_df(time = as.Date("2020-01-01") + 0:2, mats = c(1, 6, 12))

  spline_curve <- dr_build_curve(base_curve, nmats = c(1, 3, 6, 12), build = "spline", complete = TRUE)
  log_curve <- dr_build_curve(base_curve, nmats = c(1, 3, 6, 12), build = "log", complete = TRUE)
  param_curve <- dr_build_curve(base_curve, nmats = c(1, 3, 6, 12), build = "parametric", complete = FALSE)

  expected_names <- c("time", "01M", "03M", "06M", "01Y")
  expect_equal(colnames(spline_curve), expected_names)
  expect_equal(colnames(log_curve), expected_names)
  expect_equal(colnames(param_curve), expected_names)
  expect_equal(nrow(spline_curve), 3)
  expect_true(all(is.finite(as.matrix(spline_curve[, -1]))))
  expect_true(all(is.finite(as.matrix(log_curve[, -1]))))
})
