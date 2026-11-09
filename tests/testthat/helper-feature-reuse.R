# Known feature bases and covariance coefficients isolate reuse from VECM/EM
# estimation while exercising the actual NS, feature and yield-projection code.
make_feature_baseline <- function(normalize = TRUE, trunc = FALSE, smoother = NULL) {
  time <- as.Date("2020-01-01") + 0:79
  t <- seq_along(time)
  mats <- c(1, 3, 6, 12, 24, 60)
  curves <- c("usa", "gbr")
  betas <- list(
    usa = cbind(2 + .002 * t, -.5 + .04 * sin(t / 5), .2 + .03 * cos(t / 3)),
    gbr = cbind(2.3 + .003 * t + .15 * sin(t / 4), -.7 + .05 * cos(t / 6),
                .3 + .06 * sin(t / 7)))
  yields <- lapply(betas, function(beta) {
    Y <- beta %*% NS_loadings(.35, mats) + .001 * sin(outer(t, seq_along(mats)))
    panel <- data.frame(time, Y)
    names(panel) <- c("time", vapply(mats, numeric_to_matname, character(1)))
    panel
  })
  single <- st_prepare_single_curve_state(yields, .35, c(2, 4), curves, mats)
  G <- matrix(0, 6, 2, dimnames = list(NULL, c("level", "slope")))
  G[c(1, 4), 1] <- c(-1, 1)
  G[c(2, 5), 2] <- c(-1, 1)
  ecm <- list(gbr = list(rank = 2, Gmatrix = G, Fmatrix = cbind(rep(.1, 6), rep(.2, 6))))
  config <- list(lambdas = .35, cutoffs = c(2, 4), ECM_estim = "ML", ECM_type = "eigen",
    ECM_alpha = .2, X_normalize = normalize, X_trunc = trunc,
    CR_algo = zero_mean_covreg_em, CR_init = "adaptive", Xt_smooth = !is.null(smoother),
    smoother = smoother, knot_count = 3, k_count = 9, k_pass = 2,
    hp_lambda = 10, henderson_k = 13, CR_maxiter = 4, CR_tol = 1e-8,
    CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = .2,
    CR_check_every = 1, CR_use_loglik = FALSE, CR_store_path = FALSE,
    CR_verb = FALSE, CR_term = TRUE)
  xt <- st_build_xt_state(single$bl_yields, curves, "usa", reuse_ECM = TRUE,
    baseline_ecm = ecm, X_normalize = normalize, X_trunc = trunc,
    Xt_smooth = config$Xt_smooth, smoother = smoother, knot_count = 3,
    k_count = 9, k_pass = 2, hp_lambda = 10, henderson_k = 13, ytime = time)
  BS0 <- lapply(mats, function(m) list(
    B = matrix(seq_len(2 * ncol(xt$Xt)) / 20, 2,
               dimnames = list(curves, colnames(xt$Xt))),
    S0 = matrix(c(.04, .01, .01, .09), 2, dimnames = list(curves, curves))))
  names(BS0) <- paste0("X", vapply(mats, numeric_to_matname, character(1)))
  fitted <- st_finalize_shocked_fit(yields, time, single$phi_hat, xt$Xt, BS0, curves, mats, single$W)
  spec <- getFromNamespace(".cc_input_spec", "CrossCurveR")(
    yields, mats, curves, xt$Xt, xt$feature_spec)
  structure(list(curves = curves, mats = mats, reference = "usa", yields = yields,
    single_curve_state = single, xt_state = xt, BS0 = BS0,
    baseline_fit = fitted, input_spec = spec, config = config), class = "mc_shock_baseline")
}
