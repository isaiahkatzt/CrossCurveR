test_that("curve and phi cutters retain each curve's block", {
  yields <- make_pair_yields()
  embedded <- vecY(yields, fixture_curves, fixture_mats)$tY
  split <- cf_curve_cut(embedded, fixture_curves, fixture_mats)
  expect_identical(names(split), fixture_curves)
  for (curve in fixture_curves) {
    expect_equal(unname(as.matrix(split[[curve]])), unname(as.matrix(yields[[curve]][-1])))
  }

  n <- length(fixture_mats)
  loadings <- matrix(seq_len(6 * n), nrow = 2 * n)
  split_phi <- cf_phi_cut(list(loadings), fixture_curves, fixture_mats)
  expect_equal(split_phi$gbr[[1]], loadings[n + seq_len(n), ])
})
