test_that("numeric_to_matname formats monthly and yearly maturities", {
  expect_equal(numeric_to_matname(1), "01M")
  expect_equal(numeric_to_matname(6), "06M")
  expect_equal(numeric_to_matname(12), "01Y")
  expect_equal(numeric_to_matname(60), "05Y")
})

test_that("matname_to_numeric parses tenor strings into years", {
  expect_equal(matname_to_numeric("01M"), 1 / 12)
  expect_equal(matname_to_numeric("06M"), 6 / 12)
  expect_equal(matname_to_numeric("02Y"), 2)
  expect_equal(matname_to_numeric("10Y"), 10)
})

test_that("Nelson-Siegel loading helpers return expected dimensions and values", {
  lambda <- 0.5
  mats_local <- c(1, 3, 6)

  ls_loadings <- NS_loadings_ls(lambda, mats_local)
  full_loadings <- NS_loadings(lambda, mats_local)

  expect_equal(dim(ls_loadings), c(2, 3))
  expect_equal(dim(full_loadings), c(3, 3))
  expect_equal(ls_loadings[1, ], rep(1, 3))
  expect_equal(full_loadings[1, ], rep(1, 3))
  expect_equal(
    unname(full_loadings[3, 2]),
    unname((1 - exp(-mats_local[2] * lambda)) / (lambda * mats_local[2]) - exp(-lambda * mats_local[2])),
    tolerance = tol
  )
})

test_that("Svensson_loadings returns a four-factor matrix with unit level loading", {
  loadings <- Svensson_loadings(1.5, 3, c(1, 3, 6, 12))

  expect_equal(dim(loadings), c(4, 4))
  expect_equal(loadings[1, ], rep(1, 4))
  expect_true(all(is.finite(loadings)))
})
