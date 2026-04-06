test_that("full_matrix_yield strips time columns and preserves curve names", {
  yields <- make_pair_yields()
  full_matrix <- full_matrix_yield(yields, fixture_curves)

  expect_equal(names(full_matrix), fixture_curves)
  expect_equal(ncol(full_matrix$usa), length(fixture_mats))
  expect_equal(nrow(full_matrix$usa), length(fixture_time))
})

test_that("vecY builds curve-major embeddings with and without timestamps", {
  yields <- make_pair_yields()
  ts_vec <- vecY(yields, fixture_curves, fixture_mats, ts = TRUE)
  no_ts_vec <- vecY(lapply(yields, function(x) as.matrix(x[, -1])), fixture_curves, fixture_mats, ts = FALSE)

  expect_equal(length(ts_vec$time), length(fixture_time))
  expect_equal(ncol(ts_vec$tY), length(fixture_curves) * length(fixture_mats))
  expect_equal(ncol(no_ts_vec$tY), length(fixture_curves) * length(fixture_mats))
  expect_true(all(grepl("\\.", colnames(ts_vec$tY))))
})

test_that("PY_full permutes columns using the supplied permutation matrix", {
  yields <- make_pair_yields()
  vec_y <- vecY(yields, fixture_curves, fixture_mats, ts = TRUE)$tY
  P <- bln_build_P(fixture_curves, fixture_mats)
  permuted <- PY_full(vec_y, P)

  expect_equal(ncol(permuted$py), ncol(vec_y))
  expect_equal(permuted$cn, colnames(permuted$py))
  expect_false(identical(colnames(permuted$py), colnames(vec_y)))
})

test_that("curve and phi cutters split blocks by curve", {
  yields <- make_pair_yields()
  vec_y <- vecY(yields, fixture_curves, fixture_mats, ts = TRUE)$tY
  cut_y <- cf_curve_cut(vec_y, fixture_curves, fixture_mats)

  phi_list <- rep(list(matrix(seq_len(length(fixture_curves) * length(fixture_mats) * 3),
                              nrow = length(fixture_curves) * length(fixture_mats), ncol = 3)),
                  2)
  cut_phi <- cf_phi_cut(phi_list, fixture_curves, fixture_mats)

  expect_equal(names(cut_y), fixture_curves)
  expect_equal(names(cut_phi), fixture_curves)
  expect_equal(ncol(cut_y$usa), length(fixture_mats))
  expect_equal(nrow(cut_phi$gbr[[1]]), length(fixture_mats))
})
