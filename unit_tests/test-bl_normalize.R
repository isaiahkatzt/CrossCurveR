test_that("bln_build_P constructs a valid permutation matrix", {
  P <- bln_build_P(fixture_curves, fixture_mats)

  expect_equal(dim(P), c(length(fixture_curves) * length(fixture_mats), length(fixture_curves) * length(fixture_mats)))
  expect_equal(as.numeric(rowSums(P)), rep(1, nrow(P)))
  expect_equal(as.numeric(colSums(P)), rep(1, ncol(P)))
})

test_that("Hb and YB transforms preserve dimensions and identity behavior", {
  P <- bln_build_P(fixture_curves, fixture_mats)
  HT <- diag(nrow(P))
  sqrt_inv_sigmaT <- diag(nrow(P))
  Hb <- bln_Hb(HT, P, sqrt_inv_sigmaT)

  expect_equal(Hb, P %*% HT %*% t(P), tolerance = tol)

  vec_y <- vecY(make_pair_yields(), fixture_curves, fixture_mats, ts = TRUE)$tY
  yb <- bln_YB(vec_y, rep(list(diag(ncol(vec_y))), nrow(vec_y)))
  expect_equal(as.matrix(yb), as.matrix(vec_y), tolerance = tol)
})

test_that("phi stacking and permutation helpers preserve row counts and labels", {
  blsc_list <- make_mock_blsc_list()
  mat_str <- sapply(fixture_mats, numeric_to_matname)
  mds <- as.vector(t(outer(fixture_curves, mat_str, paste, sep = ".")))
  phi_hat <- bln_phi_hat(blsc_list$phi, mds)
  P <- bln_build_P(fixture_curves, fixture_mats)

  phi_bt <- bln_phi_BT(phi_hat, P, rep(list(diag(nrow(phi_hat))), 2))
  phi_ct <- bln_phi_CT(phi_bt, P)

  expect_equal(rownames(phi_hat), mds)
  expect_equal(length(phi_bt), 2)
  expect_equal(length(phi_ct), 2)
  expect_equal(nrow(phi_ct[[1]]), nrow(phi_hat))
})

test_that("bln_rescale_beta_NS and bln_normalized_fit reconstruct curve-level fits", {
  curve <- make_ns_curve_df()
  check_phi <- make_check_phi(nrow(curve))
  rescaled <- bln_rescale_beta_NS(curve[, -1], check_phi, fixture_mats)

  expect_equal(dim(rescaled$yhat), c(nrow(curve), length(fixture_mats)))
  expect_equal(dim(rescaled$betahat), c(nrow(curve), 3))

  cut_yield <- lapply(make_pair_yields(), function(x) as.matrix(x[, -1]))
  cut_phi <- list(
    usa = make_check_phi(nrow(curve)),
    gbr = make_check_phi(nrow(curve))
  )
  normalized <- bln_normalized_fit(cut_yield, cut_phi, fixture_curves, fixture_mats)

  expect_equal(names(normalized), c("n_yield", "n_beta"))
  expect_equal(names(normalized$n_yield), fixture_curves)
  expect_equal(dim(normalized$n_yield$usa), c(nrow(curve), length(fixture_mats)))
})
