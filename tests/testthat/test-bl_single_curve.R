test_that("lambda_grid_search identifies the grid minimum for exact NS data", {
  curve <- make_ns_curve_df(lambda = 0.35)
  lambda_grid <- c(0.2, 0.35, 0.5)

  expect_equal(lambda_grid_search(lambda_grid, curve, fixture_mats), 0.35)
})

test_that("phi builders produce valid cross-sectional inverses", {
  phi_obj <- bl_phi(0.35, fixture_mats)
  sven_obj <- sven_phi(1.5, 3, fixture_mats)

  expect_equal(dim(phi_obj$phi), c(3, length(fixture_mats)))
  expect_equal(dim(sven_obj$phi), c(4, length(fixture_mats)))
  expect_equal(unname(phi_obj$cross_phi %*% Matrix::t(phi_obj$phi)), Matrix::diag(3), tolerance = 1e-5)
  expect_equal(unname(sven_obj$cross_phi %*% Matrix::t(sven_obj$phi)), Matrix::diag(4), tolerance = 1e-5)
})

test_that("loading regressions report noninvertible designs with model context", {
  expect_error(bl_phi(.35, c(1, 3)),
               "NS loading regression.*lambda=0.35.*mats=.*could not be inverted.*rcond=")
  expect_error(dly_phi(.35, 1),
               "DLY country loading regression.*lambda=0.35.*could not be inverted.*rcond=")
  expect_error(sven_phi(2, 2, fixture_mats),
               "Svensson loading regression.*lambda1=2, lambda2=2.*could not be inverted.*rcond=")
})

test_that("NS and DLY country residuals agree with independent OLS residuals", {
  curve <- make_ns_curve_df(lambda = 0.35)
  rownames(curve) <- paste0("observation", seq_len(nrow(curve)))
  perturbation <- outer(seq_len(nrow(curve)), seq_along(fixture_mats),
                        function(i, j) 0.05 * sin(i * j))
  curve[, -1] <- as.matrix(curve[, -1]) + perturbation

  fits <- list(
    NS = quiet_h(bl_single_curve(0.35, curve, c(2, 4), fixture_mats)),
    DLY = dly_sc_fit(0.35, curve, fixture_mats)
  )
  observed <- as.matrix(curve[, -1])

  for (fit in fits) {
    independent <- stats::lm.fit(Matrix::t(fit$phi$phi), Matrix::t(observed))
    expected <- Matrix::t(independent$residuals)
    expect_equal(unname(fit$nsfit$W), unname(expected), tolerance = 1e-10)
    expect_equal(unname(fit$nsfit$NSfit + fit$nsfit$W), unname(observed), tolerance = 1e-10)
    expect_identical(dimnames(fit$nsfit$W), dimnames(fit$nsfit$NSfit - observed))
  }
})

test_that("bl_bdiag_H validates cutoffs and returns a square covariance block matrix", {
  residuals <- as.matrix(make_ns_curve_df()[, -1]) - mean(as.matrix(make_ns_curve_df()[, -1]))
  H <- quiet_h(bl_bdiag_H(residuals, c(2, 4), fixture_mats))

  expect_equal(dim(H), c(length(fixture_mats), length(fixture_mats)))
  expect_error(bl_bdiag_H(residuals, c(2, 6), fixture_mats), "Rate cutoffs")
  expect_error(bl_bdiag_H(residuals, "bad", fixture_mats), "Invalid rate cutoffs")
})

test_that("sc_fit returns model outputs for daily, nelson, and svensson modes", {
  curve <- make_ns_curve_df(lambda = 0.35)

  daily_fit <- sc_fit(curve, fixture_mats, type = "daily", ts = FALSE)
  nelson_fit <- sc_fit(curve, fixture_mats, type = "nelson", ts = FALSE, lambdas = c(0.2, 0.35, 0.5))
  sven_phi_obj <- sven_phi(1.5, 3, fixture_mats)
  sven_betas <- cbind(rep(2, nrow(curve)), rep(-1, nrow(curve)), rep(0.5, nrow(curve)), rep(0.2, nrow(curve)))
  sven_curve <- data.frame(time = curve$time, as.data.frame(sven_betas %*% sven_phi_obj$phi))
  colnames(sven_curve) <- c("time", sapply(fixture_mats, numeric_to_matname))
  svensson_fit <- sc_fit(sven_curve, fixture_mats, type = "svensson", ts = FALSE, slambdas = c(1.5, 3))

  expect_true(all(is.finite(daily_fit$yields)))
  expect_true(is.na(daily_fit$dynamic))
  expect_equal(unname(nelson_fit$yields), unname(as.matrix(curve[-1])), tolerance = 1e-8)
  expect_equal(unname(svensson_fit$betas), unname(sven_betas), tolerance = 1e-8)
})
