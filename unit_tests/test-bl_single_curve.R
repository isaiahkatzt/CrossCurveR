test_that("loss functions return one loss per lambda candidate", {
  curve <- make_ns_curve_df()
  daily_yield <- as.numeric(curve[1, -1])
  lambda_grid <- c(0.2, 0.35, 0.5)

  ls_loss <- l2_loss_ls(lambda_grid, daily_yield, fixture_mats)
  full_loss <- l2_loss_lambda(lambda_grid, daily_yield, fixture_mats)

  expect_equal(length(ls_loss), length(lambda_grid))
  expect_equal(length(full_loss), length(lambda_grid))
  expect_equal(which.min(full_loss), 2)
})

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
  expect_equal(unname(phi_obj$cross_phi %*% t(phi_obj$phi)), diag(3), tolerance = 1e-5)
  expect_equal(unname(sven_obj$cross_phi %*% t(sven_obj$phi)), diag(4), tolerance = 1e-5)
})

test_that("NS beta extraction and fit reconstruction recover exact NS yields", {
  curve <- make_ns_curve_df(lambda = 0.35)
  phi_obj <- bl_phi(0.35, fixture_mats)
  betas_obj <- bl_NSbetas(curve, 0.35, phi_obj$cross_phi, fixture_mats)
  fit_obj <- bl_NSfit(betas_obj$m_yield, phi_obj$phi, betas_obj$NSbetas, 0.35, fixture_mats)

  expect_equal(unname(fit_obj$NSfit), unname(betas_obj$m_yield), tolerance = 1e-5)
  expect_equal(unname(fit_obj$W), matrix(0, nrow = nrow(fit_obj$W), ncol = ncol(fit_obj$W)), tolerance = 1e-5)
})

test_that("bl_bdiag_H validates cutoffs and returns a square covariance block matrix", {
  residuals <- as.matrix(make_ns_curve_df()[, -1]) - mean(as.matrix(make_ns_curve_df()[, -1]))
  H <- bl_bdiag_H(residuals, c(2, 4), fixture_mats)

  expect_equal(dim(H), c(length(fixture_mats), length(fixture_mats)))
  expect_error(bl_bdiag_H(residuals, c(2, 6), fixture_mats), "Rate cutoffs")
  expect_error(bl_bdiag_H(residuals, "bad", fixture_mats), "Invalid rate cutoffs")
})

test_that("bl_single_curve assembles phi, betas, fit, and covariance outputs", {
  curve <- make_ns_curve_df(lambda = 0.35)
  single_curve <- bl_single_curve(0.35, curve, c(2, 4), fixture_mats)

  expect_equal(names(single_curve), c("phi", "betas", "nsfit", "H"))
  expect_equal(dim(single_curve$nsfit$NSfit), c(nrow(curve), length(fixture_mats)))
})

test_that("single-curve NS and Svensson helpers return fitted yields and betas", {
  curve <- make_ns_curve_df(lambda = 0.35)
  nt_yields <- as.matrix(curve[, -1])

  ns_daily <- sc_NSdaily(nt_yields, fixture_mats)
  ns_window <- sc_NSwindow(0.35, nt_yields, fixture_mats)

  sven_phi_obj <- sven_phi(1.5, 3, fixture_mats)
  sven_betas <- cbind(rep(2, nrow(nt_yields)), rep(-1, nrow(nt_yields)), rep(0.5, nrow(nt_yields)), rep(0.2, nrow(nt_yields)))
  sven_yields <- sven_betas %*% sven_phi_obj$phi
  sven_fit <- sc_Svensson(1.5, 3, sven_yields, fixture_mats)

  expect_equal(dim(ns_daily$NSyields), dim(nt_yields))
  expect_equal(dim(ns_window$NSbetas), c(nrow(nt_yields), 3))
  expect_equal(dim(sven_fit$Svenbetas), c(nrow(nt_yields), 4))
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

  expect_equal(names(daily_fit), c("yields", "betas", "dynamic"))
  expect_equal(dim(nelson_fit$betas), c(nrow(curve), 3))
  expect_equal(dim(svensson_fit$betas), c(nrow(curve), 4))
  expect_true(is.na(daily_fit$dynamic))
})
