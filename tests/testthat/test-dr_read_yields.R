test_that("CSV input preserves downstream names and reconstructs through each method", {
  panel <- make_ns_curve_df(time = as.Date("2020-01-01") + 0:2, mats = c(1, 6, 12))
  file <- tempfile(fileext = ".csv")
  withr::defer(unlink(file))
  utils::write.csv(panel, file, row.names = FALSE)

  prefixed <- dr_read_yields(file)
  canonical <- utils::read.csv(file, check.names = FALSE)
  expect_identical(prefixed, utils::read.csv(file, stringsAsFactors = FALSE))

  for (method in c("spline", "log", "parametric")) {
    rebuilt <- dr_build_curve(prefixed, c(1, 3, 6, 12), build = method)
    expect_equal(rebuilt, dr_build_curve(canonical, c(1, 3, 6, 12), build = method),
                 tolerance = 1e-12)
    expect_true(all(is.finite(as.matrix(rebuilt[, -1]))))
    expect_equal(as.matrix(rebuilt[, c("01M", "06M", "01Y")]),
                 round(as.matrix(canonical[, -1]), 4))
  }
})

test_that("invalid or duplicate maturity headers fail before reconstruction", {
  panel <- data.frame(time = "2020-01-01", X01M = 1, metadata = 2)
  expect_error(dr_build_curve(panel), 'Unrecognized maturity name "metadata"')

  file <- tempfile(fileext = ".csv")
  withr::defer(unlink(file))
  for (header in c('"time","01Y","01Y"', '"time","12M","01Y"')) {
    writeLines(c(header, '"2020-01-01",1,2'), file)
    expect_error(dr_read_yields(file), "duplicate tenors")
  }
})
