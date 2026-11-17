test_that("shared projection agrees with independent weighted least squares", {
  time <- fixture_time[1]
  yields <- stats::setNames(lapply(seq_along(fixture_curves), function(i) {
    panel <- make_ns_curve_df(curve_shift = .2 * i, time = time)
    panel[, -1] <- panel[, -1] + .05 * sin(seq_along(fixture_mats)^2)
    panel
  }), fixture_curves)
  Y <- vecY(yields, fixture_curves, fixture_mats)
  P <- bln_build_P(fixture_curves, fixture_mats)
  phi <- t(NS_loadings(.35, fixture_mats))
  phi_hat <- do.call(rbind, rep(list(phi), length(fixture_curves)))
  rownames(phi_hat) <- colnames(Y$tY)
  variance <- seq_len(nrow(P))
  result <- .bl_project_rescale(PY_full(Y$tY, P)$py, P, phi_hat,
    list(diag(sqrt(variance))), list(diag(1 / sqrt(variance))),
    time, fixture_curves, fixture_mats)
  expected <- unlist(lapply(seq_along(fixture_curves), function(i) {
    weights <- 1 / variance[seq.int(i, length(variance), by = length(fixture_curves))]
    stats::lm.wfit(phi, as.numeric(yields[[i]][1, -1]), w = weights)$fitted.values
  }), use.names = FALSE)
  expect_equal(as.numeric(result$curve[1, -1]), expected, tolerance = 1e-10)
  expect_equal(as.numeric(result$tenor[1, -1]), as.numeric(P %*% expected), tolerance = 1e-10)
  expect_identical(result$curve$time, time)
  expect_equal(names(result$curve), c("time", colnames(Y$tY)))
  expect_equal(names(result$tenor), c("time", colnames(Y$tY)[apply(P, 1, which.max)]))
})

test_that("single-curve preparation recovers distinct NS curves in curve and tenor order", {
  yields <- list(
    usa = make_ns_curve_df(lambda = .2),
    gbr = make_ns_curve_df(lambda = .5, curve_shift = .2)
  )
  state <- quiet_h(.bl_prepare_single_curve_state(yields, c(.2, .5), c(2, 4),
                                           fixture_curves, fixture_mats))
  expect_equal(state$lambdas, c(usa = .2, gbr = .5))
  for (curve in fixture_curves) {
    expect_equal(unname(as.matrix(state$bl_yields$nsfit[[curve]]$NSfit)),
                 unname(as.matrix(yields[[curve]][, -1])), tolerance = 1e-10)
  }
  expect_identical(rownames(state$phi_hat), colnames(vecY(yields, fixture_curves, fixture_mats)$tY))
  tenors <- paste0("X", vapply(fixture_mats, numeric_to_matname, character(1)))
  zero <- matrix(0, length(fixture_time), length(fixture_curves),
                 dimnames = list(NULL, fixture_curves))
  expect_equal(state$W, setNames(rep(list(zero), length(tenors)), tenors), tolerance = 1e-10)
})

test_that("smoothing dispatch preserves method-specific controls", {
  Xt <- cbind(x = sin(seq_along(fixture_time)^2))
  expect_equal(blce_smooth_Xt(Xt, fixture_time, "rm", k_count = 3, k_pass = 1),
               fe_xt_rm(Xt, k = 3, passes = 1))
  expect_equal(blce_smooth_Xt(Xt, fixture_time, "hp", hp_lambda = 20),
               fe_xt_hp(Xt, lambda = 20))
  expect_equal(blce_smooth_Xt(Xt, fixture_time, "henderson", henderson_k = 7),
               fe_xt_henderson(Xt, k = 7))
})

test_that("MCE and MCX spline fits preserve numeric observation indices", {
  time <- -5 + cumsum(seq_along(fixture_time)^1.3) / 10
  yields <- lapply(make_pair_yields(), function(panel) { panel$time <- time; panel })
  Xt <- cbind(spread = sin(seq_along(time)))
  covariance <- function(W, X, ...) {
    list(S0 = diag(ncol(W)), B = matrix(0, ncol(W), ncol(X)))
  }
  for (kind in c("ns", "bs")) {
    args <- list(yields = yields, lambdas = .35, cutoffs = c(2, 4), reference = "usa",
                 curves = fixture_curves, mats = fixture_mats, Xt_smooth = TRUE,
                 smoother = kind, knot_count = 3, CR_algo = covariance)
    exo <- quiet_h(do.call(mc_fit_exo, c(args, list(Xt = Xt))))
    end <- with_temp_bindings(list(blcc_cointegration = function(...) {
      list(cc_cspread = list(gbr = Xt), cc_ecm = list(gbr = list(rank = 1)), cc_beta = list())
    }), quiet_h(do.call(mc_fit_end, args)))
    for (fit in list(exo, end)) {
      expect_identical(fit$curve$time, time)
      expect_identical(fit$tenor$time, time)
      expect_equal(fit$Xt, blce_smooth_Xt(Xt, time, kind, knot_count = 3))
    }
  }
})

test_that("MCE passes the requested shrinkage to the tenor covariance estimator", {
  baseline <- make_feature_baseline()
  observed_shrinkage <- numeric()
  covariance <- function(W, X, S0_shrink_diag, ...) {
    observed_shrinkage <<- c(observed_shrinkage, S0_shrink_diag)
    list(S0 = diag(ncol(W)), B = matrix(0, ncol(W), ncol(X)))
  }
  ecm <- baseline$xt_state$cc_ci_features$cc_ecm$gbr
  with_temp_bindings(list(blmc_ECM = function(...) ecm), {
    quiet_h(mc_fit_end(baseline$yields, .35, c(2, 4), "usa", CR_algo = covariance,
                       CR_S0_shrink_diag = .6, curves = baseline$curves, mats = baseline$mats))
  })
  expect_equal(observed_shrinkage, rep(.6, length(baseline$mats)))
})
