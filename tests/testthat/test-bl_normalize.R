test_that("the permutation converts curve-major values to tenor order", {
  P <- bln_build_P(fixture_curves, fixture_mats)
  indices <- matrix(seq_len(ncol(P)), nrow = length(fixture_mats))
  expect_equal(as.numeric(P %*% as.vector(indices)), as.vector(t(indices)))
})

test_that("transformed loading errors identify the failing observation", {
  curve <- make_ns_curve_df(time = fixture_time[1:2])
  phi <- make_check_phi(2)
  phi[[2]][, 3] <- 0
  expect_error(bln_rescale_beta_NS(curve[, -1], phi, fixture_mats),
               "Transformed NS.*observation 2.*mats=.*could not be inverted.*rcond=")
  phi[[2]][1, 1] <- Inf
  expect_error(bln_rescale_beta_NS(curve[, -1], phi, fixture_mats),
               "Transformed NS.*observation 2.*matrix contains nonfinite values")
})
