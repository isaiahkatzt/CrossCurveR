test_that("yield embeddings match curve and maturity identities rather than positions", {
  yields <- make_pair_yields()
  expected <- vecY(yields, fixture_curves, fixture_mats)
  shuffled <- yields[rev(fixture_curves)]
  shuffled <- lapply(shuffled, function(panel) panel[c(1, 7:2)])
  names(shuffled$usa)[names(shuffled$usa) == "01Y"] <- "12M"
  actual <- vecY(shuffled, fixture_curves, fixture_mats)
  expect_identical(actual, expected)
  expect_identical(vecY(unname(yields), fixture_curves, fixture_mats), expected)
  expect_identical(st_validate_yields_input(unname(yields), fixture_curves, "usa"), yields)
  expect_error(mc_fit_end(yields, .35, c(2, 4), NULL,
                          curves = fixture_curves, mats = fixture_mats), "reference")
  expect_error(dly_fit(yields, fixture_mats, .35, reference = character()), "reference")
  expect_identical(vecY(yields, rev(fixture_curves), fixture_mats)$tY[, 1:6],
                   stats::setNames(expected$tY[, 7:12], names(expected$tY)[7:12]))

  matrices <- lapply(shuffled, function(x) as.matrix(x[-1]))
  expect_identical(vecY(matrices, fixture_curves, fixture_mats, ts = FALSE)$tY,
                   vecY(lapply(yields, function(x) as.matrix(x[-1])),
                        fixture_curves, fixture_mats, ts = FALSE)$tY)
})

test_that("misaligned, duplicate, missing and unordered observations fail before estimation", {
  for (change in list(
    function(x) { x$time <- x$time + 1; x },
    function(x) x[rev(seq_len(nrow(x))), ],
    function(x) { x$time[2] <- x$time[1]; x },
    function(x) { x$time[2] <- NA; x },
    function(x) x[-1, ]
  )) {
    yields <- make_pair_yields()
    yields$gbr <- change(yields$gbr)
    expect_error(vecY(yields, fixture_curves, fixture_mats), "index")
    expect_error(dly_fit(yields, fixture_mats, lambdas = .35), "index")
    expect_error(mc_fit_end(yields, .35, c(2, 4), "usa",
                            curves = fixture_curves, mats = fixture_mats), "index")
  }
  yields <- make_pair_yields()
  yields$gbr$time <- as.character(yields$gbr$time)
  expect_identical(vecY(yields, fixture_curves, fixture_mats)$time, fixture_time)
  yields$gbr$time <- as.POSIXct(fixture_time, tz = "UTC")
  expect_identical(vecY(yields, fixture_curves, fixture_mats)$time, fixture_time)
})

test_that("curve identities and economic tenors cannot be duplicated or guessed", {
  yields <- make_pair_yields()
  expect_error(vecY(yields, c("usa", "usa"), fixture_mats), "unique")
  expect_error(vecY(yields, c("usa", "jpn"), fixture_mats), "match")
  names(yields) <- c("usa", "usa")
  expect_error(vecY(yields, fixture_curves, fixture_mats), "unique")
  yields <- make_pair_yields()
  names(yields$gbr)[2] <- "12M"
  expect_error(vecY(yields, fixture_curves, fixture_mats), "match")
  names(yields$gbr)[2] <- "unknown"
  expect_error(vecY(yields, fixture_curves, fixture_mats), "Unrecognized maturity")
  yields <- make_pair_yields()
  yields$gbr <- yields$gbr[-2]
  expect_error(vecY(yields, fixture_curves, fixture_mats), "one column per maturity")
})

test_that("finite numeric yields allow zero and negative rates without rescaling", {
  yields <- make_pair_yields()
  yields$usa[, -1] <- yields$usa[, -1] - 3
  yields$usa[1, 2] <- 0
  expect_identical(as.numeric(vecY(yields, fixture_curves, fixture_mats)$tY[1, 1:6]),
                   as.numeric(yields$usa[1, -1]))
  for (bad in list(NA_real_, NaN, Inf, "2.5", factor("2.5"))) {
    invalid <- make_pair_yields()
    invalid$usa[[2]] <- rep(bad, nrow(invalid$usa))
    expect_error(dly_sc_fit_all(.35, invalid, fixture_mats), "numeric")
  }
})

test_that("maturity and fixed-decay validation does not change parameter conventions", {
  curve <- make_ns_curve_df()
  for (bad in list(c(1, 1, 3), c(3, 1), c(0, 1), c(1, Inf), numeric())) {
    expect_error(sc_fit(curve, bad, ts = FALSE), "mats")
  }
  for (bad in list(0, -.1, NA_real_, Inf, numeric())) {
    expect_error(lambda_grid_search(bad, curve, fixture_mats), "positive")
    expect_error(NS_loadings(bad, fixture_mats), "positive")
  }
  expect_error(bl_bdiag_H(as.matrix(curve[-1]), c(1.5, 4), fixture_mats), "cutoffs")
  expect_error(bl_bdiag_H(as.matrix(curve[-1]), c(1, Inf), fixture_mats), "cutoffs")
  expect_equal(NS_loadings(.0609, fixture_mats),
               NS_loadings(12 * .0609, fixture_mats / 12), tolerance = 1e-14)
  expect_equal(Svensson_loadings(18, 36, fixture_mats),
               Svensson_loadings(18 / 12, 36 / 12, fixture_mats / 12), tolerance = 1e-14)
  expect_identical(numeric_to_matname(12), "01Y")
  expect_equal(matname_to_numeric("01Y"), 1)
})

test_that("precomputed DLY fits reconcile complete maturity layouts and observed dates", {
  yields <- make_pair_yields()
  sc <- dly_sc_fit_all(.35, yields, fixture_mats)
  expected <- dly_fit_from_sc(sc, yields, mats = fixture_mats)
  shuffled <- sc[rev(fixture_curves)]
  for (curve in fixture_curves) {
    order <- rev(seq_along(fixture_mats))
    shuffled[[curve]]$yield_names <- shuffled[[curve]]$yield_names[order]
    shuffled[[curve]]$phi <- lapply(shuffled[[curve]]$phi, function(x) x[, order])
    shuffled[[curve]]$betas$m_yield <- shuffled[[curve]]$betas$m_yield[, order]
    shuffled[[curve]]$nsfit <- lapply(shuffled[[curve]]$nsfit, function(x) x[, order])
  }
  expect_identical(dly_fit_from_sc(shuffled, yields, curves = fixture_curves), expected)
  invalid <- sc
  invalid$gbr$time <- invalid$gbr$time + 1
  expect_error(dly_fit_from_sc(invalid), "index")
  invalid <- sc
  invalid$gbr$phi$phi <- invalid$gbr$phi$phi[, -1]
  expect_error(dly_fit_from_sc(invalid), "loadings")
  invalid_yields <- yields
  invalid_yields <- lapply(invalid_yields, function(x) { x$time <- x$time + 1; x })
  expect_error(dly_fit_from_sc(sc, invalid_yields), "observation index")
  unnamed <- dly_fit(unname(yields), fixture_mats, .35, curves = fixture_curves)
  expect_identical(unnamed$yields, dly_fit(yields, fixture_mats, .35)$yields)
})

test_that("MCX verifies indexed feature rows and preserves its supplied design", {
  yields <- make_pair_yields()
  X <- cbind(intercept = 1, state = seq(-2, 2, length.out = length(fixture_time)))
  rownames(X) <- as.character(fixture_time)
  fit <- function(X, normalize = TRUE) mc_fit_exo(yields, .35, c(2, 4), "usa", X,
    X_normalize = normalize, CR_algo = zero_mean_covreg_em, CR_maxiter = 2,
    CR_Binit = matrix(.01, 2, 2), CR_S0init = diag(.01, 2),
    curves = fixture_curves, mats = fixture_mats)
  baseline <- fit(X)
  expect_identical(baseline$Xt, X)
  expect_identical(fit(X, FALSE), baseline)
  expect_identical(baseline$input_spec$features$columns, colnames(X))
  expect_identical(baseline$input_spec$units$yield, "percent")
  expect_identical(baseline$input_spec$time, fixture_time)
  expect_equal(baseline$input_spec$decay$selected, c(usa = .35, gbr = .35))
  yields <- lapply(yields[rev(fixture_curves)], function(x) x[c(1, 7:2)])
  expect_identical(fit(X), baseline)
  positional <- X
  rownames(positional) <- as.character(101:120)
  positional_fit <- fit(positional)
  expect_identical(positional_fit$Xt, positional)
  expect_identical(positional_fit$input_spec$features$row_alignment, "positional")
  wrong <- X
  rownames(wrong) <- as.character(fixture_time + 1)
  expect_error(fit(wrong), "index")
  expect_error(fit(X[-1, ]), "one row")
  expect_error(fit(array(1, c(2, 2, 2))), "two dimensional")
  wrong <- X
  wrong[1, 2] <- Inf
  expect_error(fit(wrong), "finite")
  colnames(wrong) <- c("state", "state")
  wrong[1, 2] <- 0
  expect_error(fit(wrong), "unique")
})

test_that("covariance initializers and prediction features use named identities", {
  W <- cbind(usa = sin(1:20), gbr = cos(1:20))
  X <- cbind(intercept = 1, state = seq(-1, 1, length.out = 20))
  B <- matrix(c(.1, .2, .3, -.1), 2, dimnames = list(colnames(W), colnames(X)))
  S0 <- matrix(c(.4, .1, .1, .5), 2, dimnames = list(colnames(W), colnames(W)))
  baseline <- zero_mean_covreg_em(W, X, B = B, S0 = S0, max_iter = 2)
  reordered <- zero_mean_covreg_em(W, X, B = B[2:1, 2:1], S0 = S0[2:1, 2:1], max_iter = 2)
  expect_identical(reordered, baseline)
  expect_identical(cr_build(baseline, X), cr_build(baseline, X[, 2:1]))
  invalid <- X
  colnames(invalid) <- c("different", "state")
  expect_error(cr_build(baseline, invalid), "feature names")
  expect_error(zero_mean_covreg_em(W, X, S0 = diag(3)), "one row and column")
  expect_error(zero_mean_covreg_em(W, X, S0 = -diag(2)), "positive semidefinite")
  expect_error(zero_mean_covreg_em(W, X, B = matrix(1, 3, 2)), "dimensions")
  rownames(W) <- as.character(fixture_time)
  rownames(X) <- as.character(fixture_time + 1)
  expect_error(zero_mean_covreg_em(W, X), "row indices")
})

test_that("feature specifications retain the fitted basis and normalization parameters", {
  beta <- cbind(reference = 1:20, country = 21:40)
  basis <- matrix(c(1, -.5), 2, 1)
  raw <- as.numeric(beta %*% basis)
  features <- list(cc_beta = list(gbr = beta),
    cc_ecm = list(gbr = list(rank = 1, Gmatrix = basis)),
    cc_cspread = list(gbr = matrix(raw, ncol = 1, dimnames = list(NULL, "gbr.r1"))))
  spec <- .bl_feature_spec(features, TRUE, FALSE, FALSE, NULL, 5, 9, 3, 1600, 13)
  expect_identical(spec$relationships$gbr$basis, basis)
  expect_equal(unname(spec$relationships$gbr$center), mean(raw))
  expect_equal(unname(spec$relationships$gbr$scale), stats::sd(raw))
  expect_identical(spec$smoothing$enabled, FALSE)
})
