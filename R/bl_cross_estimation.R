#############################################
#         Full Cross-Curve Estimation       #
#############################################

blce_smooth_Xt <- function(Xt, ytime, smoother = c("ns", "bs", "rm", "hp", "henderson"),
                           knot_count = 5, k_count = 9, k_pass = 3,
                           hp_lambda = 1600, henderson_k = 13) {
  smoother <- match.arg(smoother)

  if (smoother == "rm") {
    return(fe_xt_rm(Xt, k_count, k_pass))
  }

  if (smoother == "hp") {
    return(fe_xt_hp(Xt, lambda = hp_lambda))
  }

  if (smoother == "henderson") {
    return(fe_xt_henderson(Xt, k = henderson_k))
  }

  kps <- fe_kp_quantile(time = ytime, knot_count = knot_count)
  apply(Xt, 2, fe_xt_spline, time = ytime, spline_fit = smoother, knot_points = kps)
}

blce_rescale_yield <- function(normalized, sqrt_sigma_t, P, ptype = c("tenor", "curve"), time, curves, mats) {
  ptype = match.arg(ptype)
  vec_normalized <- vecY(normalized, curves, mats, ts=FALSE)[['tY']]
  N <- nrow(vec_normalized)
  rescale_yield_matrix <- matrix(NA, nrow = N, ncol = ncol(vec_normalized))
  ## rescaling and curve / tenor permutation

  if (ptype == "tenor") {
    permuted_colnames <- c("time", colnames(vec_normalized)[apply(P, 1, which.max)])

    for (dt in seq_len(N)) {
      rescale_yield_matrix[dt, ] <-
        as.numeric(Matrix::t(sqrt_sigma_t[[dt]] %*% P %*% as.numeric(vec_normalized[dt, ])))
    }
  }

  else if (ptype == "curve") {
    permuted_colnames <- c("time", colnames(vec_normalized))

    for (dt in seq_len(N)) {
      rescale_yield_matrix[dt, ] <-
        as.numeric(Matrix::t(Matrix::t(P) %*% sqrt_sigma_t[[dt]] %*% P %*% as.numeric(vec_normalized[dt, ])))
    }
  }

  rescale_yield <- data.frame(time, rescale_yield_matrix)
  colnames(rescale_yield) <- permuted_colnames
  return(rescale_yield)
}

.bl_prepare_single_curve_state <- function(yields, lambdas, cutoffs, curves, mats) {
  mat_str <- sapply(mats, numeric_to_matname)
  mds <- as.vector(Matrix::t(outer(curves, mat_str, paste, sep = ".")))
  curve_lambdas <- vapply(yields, function(curve_yield) {
    lambda_grid_search(lambdas, curve_yield, mats)
  }, numeric(1))
  curve_fits <- stats::setNames(lapply(seq_along(curves), function(i) {
    bl_single_curve(curve_lambdas[[i]], yields[[i]], cutoffs, mats)
  }), curves)
  bl_yields <- purrr::transpose(curve_fits)

  list(
    mds = mds,
    lambdas = curve_lambdas,
    bl_yields = bl_yields,
    phi_hat = bln_phi_hat(bl_yields$phi, mds),
    W = blcc_build_Wj(bl_yields$nsfit, mats),
    H = Matrix::bdiag(bl_yields$H)
  )
}

.bl_project_rescale <- function(py, P, phi_hat, sqrt_sigma_t, sqrt_inv_sigma_t,
                                time, curves, mats) {
  yb <- bln_YB(py, sqrt_inv_sigma_t)
  phib <- bln_phi_BT(phi_hat, P, sqrt_inv_sigma_t)
  phic <- bln_phi_CT(phib, P)
  yc <- PY_full(yb, Matrix::t(P))

  ycut <- cf_curve_cut(yc$py, curves, mats)
  phicut <- cf_phi_cut(phic, curves, mats)
  normalized <- bln_normalized_fit(ycut, phicut, curves, mats)

  tenor <- blce_rescale_yield(normalized$n_yield, sqrt_sigma_t, P, "tenor",
                            time = time, curves = curves, mats = mats)
  curve <- blce_rescale_yield(normalized$n_yield, sqrt_sigma_t, P, "curve",
                            time = time, curves = curves, mats = mats)
  list(tenor = tenor, curve = curve)
}

mc_fit_end <- function(yields, lambdas, cutoffs, reference,
                       ECM_estim="ML", ECM_type="eigen", ECM_alpha=0.1, X_normalize=TRUE, X_trunc=FALSE,
                       CR_algo=zero_mean_covreg_em, CR_init="adaptive", Xt_smooth=FALSE, smoother=c("ns", "bs", "rm", "hp", "henderson"),
                       knot_count=5, k_count=9, k_pass=3, hp_lambda=1600, henderson_k=13,
                       CR_maxiter=1000, CR_tol=1e-8, CR_Binit=NULL, CR_S0init=NULL, CR_S0_shrink_diag=0,
                       CR_check_every=1, CR_use_loglik=TRUE, CR_store_path=TRUE,
                       CR_verb=FALSE, CR_term=TRUE,
                       curves, mats) {
  ## full endogenous covariate fit

  ## build vector yield embedding and permutation matrix
  Y <- vecY(yields, curves, mats)
  ytime <- Y[['time']]

  P <- bln_build_P(curves, mats)
  pY <- PY_full(Y[['tY']], P)

  ## single-curve NS fits, loadings, and residual covariance
  single_curve <- .bl_prepare_single_curve_state(yields, lambdas, cutoffs, curves, mats)
  bl_yields <- single_curve$bl_yields
  phi_hat <- single_curve$phi_hat
  W <- single_curve$W
  H <- single_curve$H

  ## feature extraction
  cc_ci_features <- blcc_cointegration(blsc_list=bl_yields, curves=curves, reference=reference,
                                       estim=ECM_estim, type=ECM_type, alpha=ECM_alpha, normalize=X_normalize)
  Xt <- blcc_build_Xt(cc_ci_features$cc_cspread, cc_ci_features$cc_ecm, trunc=X_trunc)

  if (Xt_smooth){
    Xt <- blce_smooth_Xt(
      Xt = Xt, ytime = ytime, smoother = smoother, knot_count = knot_count,
      k_count = k_count, k_pass = k_pass, hp_lambda = hp_lambda, henderson_k = henderson_k
    )
  }

  ## covariance regression and sigma estimation
  cc_BS0_feature <- blcc_fe(Xt=Xt, Wt=W, mats=mats, covreg=CR_algo, init=CR_init,
                            max_iter=CR_maxiter, tol=CR_tol, S0=CR_S0init, B=CR_Binit,
                            S0_shrink_diag=CR_S0_shrink_diag, verb=CR_verb, term=CR_term,
                            check_every=CR_check_every, use_loglik=CR_use_loglik,
                            store_path=CR_store_path)

  cc_Sigma <- tSigma_optim(BS0=cc_BS0_feature, X=Xt)

  full_sigma_t <- sqrt_inv_build(cc_Sigma$tSigma)

  sqrt_sigma_t <- full_sigma_t$sqrt
  sqrt_inv_sigma_t <- full_sigma_t$inverse

  # Retain the existing H transformations pending their separate mathematical audit.
  hb <- lapply(seq_along(sqrt_inv_sigma_t), function(dt) {
    bln_Hb(H, P, sqrt_inv_sigma_t[[dt]])
  })
  hc <- lapply(hb, function(hb_dt) Matrix::t(P) %*% hb_dt %*% P)

  projected <- .bl_project_rescale(pY$py, P, phi_hat, sqrt_sigma_t, sqrt_inv_sigma_t,
                                   time = ytime, curves = curves, mats = mats)

  return(list(
    tenor=projected$tenor,
    curve=projected$curve,
    sigma_JT=cc_Sigma,
    curveECM=cc_ci_features[['cc_ecm']],
    BS0=cc_BS0_feature,
    Xt=Xt,
    W=W,
    ns_factor=bl_yields[['betas']])
    )
}

mc_fit_exo <- function(yields, lambdas, cutoffs, reference, Xt, X_normalize=TRUE,
                       CR_algo=zero_mean_covreg, CR_init="adaptive", Xt_smooth=FALSE, smoother=c("ns", "bs", "rm", "hp", "henderson"),
                       knot_count=5, k_count=9, k_pass=3, hp_lambda=1600, henderson_k=13,
                       CR_maxiter=1000, CR_tol=1e-8, CR_Binit=NULL, CR_S0init=NULL, CR_S0_shrink_diag=0,
                       CR_check_every=1, CR_use_loglik=TRUE, CR_store_path=TRUE,
                       CR_verb=FALSE, CR_term=TRUE,
                       curves, mats) {
  ## full exogenous covariate fit

  ## build vector yield embedding and permutation matrix
  Y <- vecY(yields, curves, mats)
  ytime <- Y[['time']]

  P <- bln_build_P(curves, mats)
  pY <- PY_full(Y[['tY']], P)

  ## single-curve NS fits, loadings, and residual covariance
  single_curve <- .bl_prepare_single_curve_state(yields, lambdas, cutoffs, curves, mats)
  phi_hat <- single_curve$phi_hat
  W <- single_curve$W
  H <- single_curve$H

  if (Xt_smooth){
    Xt <- blce_smooth_Xt(
      Xt = Xt, ytime = ytime, smoother = smoother, knot_count = knot_count,
      k_count = k_count, k_pass = k_pass, hp_lambda = hp_lambda, henderson_k = henderson_k
    )
  }

  ## covariance regression and sigma estimation
  cc_BS0_feature <- blcc_fe(Xt=Xt, Wt=W, mats=mats, covreg=CR_algo, init=CR_init,
                            max_iter=CR_maxiter, tol=CR_tol, S0=CR_S0init, B=CR_Binit,
                            S0_shrink_diag=CR_S0_shrink_diag, verb=CR_verb, term=CR_term,
                            check_every=CR_check_every, use_loglik=CR_use_loglik,
                            store_path=CR_store_path)

  cc_Sigma <- tSigma_optim(BS0=cc_BS0_feature, X=Xt)

  full_sigma_t <- sqrt_inv_build(cc_Sigma$tSigma)
  sqrt_sigma_t <- full_sigma_t$sqrt
  sqrt_inv_sigma_t <- full_sigma_t$inverse

  # Retain the existing H transformations pending their separate mathematical audit.
  hb <- lapply(seq_along(sqrt_inv_sigma_t), function(dt) {
    bln_Hb(H, P, sqrt_inv_sigma_t[[dt]])
  })
  hc <- lapply(hb, function(hb_dt) Matrix::t(P) %*% hb_dt %*% P)

  projected <- .bl_project_rescale(pY$py, P, phi_hat, sqrt_sigma_t, sqrt_inv_sigma_t,
                                   time = ytime, curves = curves, mats = mats)

  return(list(
    tenor=projected$tenor,
    curve=projected$curve,
    sigma_JT=cc_Sigma,
    BS0=cc_BS0_feature,
    Xt=Xt,
    W=W))
}
