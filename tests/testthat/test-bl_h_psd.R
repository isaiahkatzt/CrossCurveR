test_that("indefinite within-curve H warns without changing its construction", {
  W <- cbind(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1))
  expected <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), 3)
  for (scale in c(1e-6, 1, 1e6)) {
    result <- capture_h_warnings(bl_bdiag_H(scale * W, c(1, 2), c(1, 3, 6)))
    expect_length(result$warnings, 1)
    w <- result$warnings[[1]]
    expect_s3_class(w, "crosscurver_h_not_psd")
    expect_equal(result$value, expected * scale^2, tolerance = 1e-14)
    expect_equal(w$minimum_eigenvalue / scale^2, 1 - sqrt(2), tolerance = 1e-12)
    expect_equal(w$psd_tolerance / scale^2, 1e-8 * (1 + sqrt(2)), tolerance = 1e-12)
    expect_true(w$psd_assessed)
    expect_identical(w$cutoffs, c(1, 2))
    expect_identical(w$mats, c(1, 3, 6))
    expect_equal(w$n_observations, 3)
    expect_match(conditionMessage(w), "H is not PSD under this cutoff configuration")
    expect_match(conditionMessage(w), "Returning H unchanged")
  }
})

test_that("PSD and negligible negative eigenvalues do not warn", {
  # Centered design with unit covariance, used to generate known covariances.
  base <- sqrt(3 / 4) * rbind(c(1, 1, 1), c(1, -1, -1),
                             c(-1, 1, -1), c(-1, -1, 1))
  S <- matrix(c(1, .5, 0, .5, 1, .5, 0, .5, 1), 3)
  expect_silent(H <- bl_bdiag_H(base %*% chol(S), c(1, 2), c(1, 3, 6)))
  expect_equal(H, S, tolerance = 1e-14)
  expect_silent(H <- bl_bdiag_H(matrix(0, 4, 3), c(1, 2), c(1, 3, 6)))
  expect_identical(H, matrix(0, 3, 3))

  S[1, 3] <- S[3, 1] <- 1e-5
  S[1, 2] <- S[2, 1] <- S[2, 3] <- S[3, 2] <- (1 + 1e-10) / sqrt(2)
  for (scale in c(1e-6, 1, 1e6)) {
    expect_silent(H <- bl_bdiag_H(scale * base %*% chol(S), c(1, 2), c(1, 3, 6)))
    values <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    expect_lt(min(values), 0)
    expect_gte(min(values), -1e-8 * max(abs(values)))
  }
})

test_that("unassessable H warns and retains nonfinite entries", {
  result <- capture_h_warnings(bl_bdiag_H(matrix(1, 1, 3), c(1, 2), c(1, 3, 6)))
  expect_length(result$warnings, 1)
  expect_false(result$warnings[[1]]$psd_assessed)
  expect_true(is.na(result$warnings[[1]]$minimum_eigenvalue))
  expect_true(all(is.na(diag(result$value))))
  expect_match(conditionMessage(result$warnings[[1]]), "PSD could not be assessed")
})

test_that("multi-curve H warnings identify the curve without duplication", {
  mats <- c(1, 3, 6, 12, 24, 60)
  phi <- t(NS_loadings(.06, mats))
  v <- qr.resid(qr(phi), c(.2, -.4, .6, -.8, .3, 1))
  Y <- rbind(c(2, -.8, .4), c(2.1, -.7, .3), c(2.2, -.6, .2)) %*% t(phi) +
    outer(c(-1, 0, 1), v)
  panel <- data.frame(time = as.Date("2020-01-01") + 0:2, Y)
  names(panel) <- c("time", vapply(mats, numeric_to_matname, character(1)))
  curves <- c("usa", "gbr")
  result <- capture_h_warnings(.bl_prepare_single_curve_state(
    list(usa = panel, gbr = panel), .06, c(2, 4), curves, mats))
  expect_length(result$warnings, 2)
  expect_identical(vapply(result$warnings, `[[`, character(1), "curve"), curves)
  expect_match(conditionMessage(result$warnings[[1]]), "^Curve 'usa':")
  expect_match(conditionMessage(result$warnings[[2]]), "^Curve 'gbr':")
  expect_true(all(vapply(result$warnings, function(w) w$minimum_eigenvalue < 0, logical(1))))
  expect_identical(names(result$value$bl_yields$H), curves)
})
