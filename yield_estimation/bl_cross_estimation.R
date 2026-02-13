source("yield_estimation/curve_reformat.R")
source("yield_estimation/bl_single_curve.R")
source("yield_estimation/bl_feature_extract.R")
source("yield_estimation/cr_core_estim.R")
source("yield_estimation/matrix_helpers.R")
source("yield_estimation/bl_normalize.R")

#############################################
#         Full Cross-Curve Estimation       #
#############################################

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
        as.numeric(t(sqrt_sigma_t[[dt]] %*% P %*% as.numeric(vec_normalized[dt, ])))
    }
  }
  
  else if (ptype == "curve") {
    permuted_colnames <- c("time", colnames(vec_normalized))
    
    for (dt in seq_len(N)) {
      rescale_yield_matrix[dt, ] <- 
        as.numeric(t(t(P) %*% sqrt_sigma_t[[dt]] %*% P %*% as.numeric(vec_normalized[dt, ])))
    }
  }
  
  rescale_yield <- data.frame(time, rescale_yield_matrix)
  colnames(rescale_yield) <- permuted_colnames
  return(rescale_yield) 
}

mc_fit_end <- function(yields, lambdas, cutoffs, reference, 
                       ECM_estim="ML", ECM_type="eigen", ECM_alpha=0.1, X_normalize=TRUE, X_trunc=FALSE,
                       CR_algo=zero_mean_covreg, CR_init="adaptive", Xt_smooth=FALSE, smoother=c("ns", "bs"), knot_count=5,
                       CR_maxiter=1000, CR_tol=1e-8, CR_Binit=NULL, CR_S0init=NULL, CR_verb=FALSE, CR_term=TRUE, 
                       curves, mats) {
  ## full endogenous covariate fit 
  
  ## build maturity x tenor strings 
  mat_str <- sapply(mats, numeric_to_matname)
  mds <- as.vector(t(outer(curves, mat_str, paste, sep='.')))
  
  ## build vector yield embedding and permutation matrix  
  Y <- vecY(yields, curves, mats)   
  ytime <- Y[['time']]
  
  P <- bln_build_P(curves, mats)  
  pY <- PY_full(Y[['tY']], P)
  
  ## grid search lambda construction 
  curve_lambdas <- vapply(yields, function(x) { 
    lambda_grid_search(lambdas, x, mats)
  }, numeric(1) 
  )
  
  ## NS baseline estimation (no cross-curve) 
  bl_yields <- transpose(setNames(lapply(seq_along(yields), function(x) {
    bl_single_curve(curve_lambdas[x], yields[[x]], cutoffs, mats)
  }
  ), curves)) 

  ## residuals by tenor; within-curve full covariance  
  phi_hat <- bln_phi_hat(bl_yields$phi, mds) 
  W <- blcc_build_Wj(bl_yields$nsfit, mats) 
  H <- bdiag(bl_yields$H) 

  ## feature extraction 
  cc_ci_features <- blcc_cointegration(blsc_list=bl_yields, curves=curves, reference=reference, 
                                       estim=ECM_estim, type=ECM_type, alpha=ECM_alpha, normalize=X_normalize)
  Xt <- blcc_build_Xt(cc_ci_features$cc_cspread, cc_ci_features$cc_ecm, trunc=X_trunc) 
  
  if (Xt_smooth){
    kps <- fe_kp_quantile(time=ytime, knot_count=knot_count)
    Xt <- apply(Xt, 2, fe_xt_spline, time=ytime, spline_fit=smoother, knot_points=kps)
  }
  
  ## covariance regression and sigma estimation 
  cc_BS0_feature <- blcc_fe(Xt=Xt, Wt=W, mats=mats, covreg=CR_algo, init=CR_init, 
                            max_iter=CR_maxiter, tol=CR_tol, S0=CR_S0init, B=CR_Binit, verb=CR_verb, term=CR_term) 
  
  cc_Sigma <- tSigma_optim(BS0=cc_BS0_feature, X=Xt)
  
  full_sigma_t <- sqrt_inv_build(cc_Sigma$tSigma)
  sqrt_sigma_t <- full_sigma_t$sqrt
  sqrt_inv_sigma_t <- full_sigma_t$inverse
  
  ## breve components 
  yb <- bln_YB(pY[['py']], sqrt_inv_sigma_t) 
  phib <- bln_phi_BT(phi_hat, P, sqrt_inv_sigma_t)
  hb <- lapply(seq_along(sqrt_inv_sigma_t), function(dt) {
    bln_Hb(H, P, sqrt_inv_sigma_t[[dt]])
  })
  
  ## check components 
  phic <- bln_phi_CT(phib, P) 
  yc <- PY_full(yb, t(P))
  hc <- lapply(hb, function(hb_dt) t(P) %*% hb_dt %*% P)
  
  ## curve cuts  
  ycut <- cf_curve_cut(yc[['py']], curves, mats)
  phicut <- cf_phi_cut(phic, curves, mats)
  
  ## normalized 
  normalized <- bln_normalized_fit(ycut, phicut, curves, mats)
  ystar_tenor <- blce_rescale_yield(normalized[['n_yield']], sqrt_sigma_t, P, "tenor", time=ytime, curves=curves, mats=mats)
  ystar_curve <- blce_rescale_yield(normalized$n_yield, sqrt_sigma_t, P, "curve", time=ytime, curves=curves, mats=mats)
  
  return(list(
    tenor=ystar_tenor, 
    curve=ystar_curve, 
    sigma_JT=cc_Sigma,
    curveECM=cc_ci_features[['cc_ecm']], 
    BS0=cc_BS0_feature, 
    Xt=Xt))  
}

mc_fit_exo <- function(yields, lambdas, cutoffs, reference, Xt, X_normalize=TRUE,
                       CR_algo=zero_mean_covreg, CR_init="adaptive", Xt_smooth=FALSE, smoother=c("ns", "bs"), knot_count=5,
                       CR_maxiter=1000, CR_tol=1e-8, CR_Binit=NULL, CR_S0init=NULL, CR_verb=FALSE, CR_term=TRUE, 
                       curves, mats) {
  ## full exogenous covariate fit 
  
  ## build maturity x tenor strings 
  mat_str <- sapply(mats, numeric_to_matname)
  mds <- as.vector(t(outer(curves, mat_str, paste, sep='.')))
  
  ## build vector yield embedding and permutation matrix  
  Y <- vecY(yields, curves, mats)   
  ytime <- Y[['time']]
  
  P <- bln_build_P(curves, mats)  
  pY <- PY_full(Y[['tY']], P)
  
  ## grid search lambda construction 
  curve_lambdas <- vapply(yields, function(x) { 
    lambda_grid_search(lambdas, x, mats)
  }, numeric(1) 
  )
  
  ## NS baseline estimation (no cross-curve) 
  bl_yields <- transpose(setNames(lapply(seq_along(yields), function(x) {
    bl_single_curve(curve_lambdas[x], yields[[x]], cutoffs, mats)
  }
  ), curves)) 
  
  ## residuals by tenor; within-curve full covariance  
  phi_hat <- bln_phi_hat(bl_yields$phi, mds) 
  W <- blcc_build_Wj(bl_yields$nsfit, mats) 
  H <- bdiag(bl_yields$H) 
  
  if (Xt_smooth){
    kps <- fe_kp_quantile(time=ytime, knot_count=knot_count)
    Xt <- apply(Xt, 2, fe_xt_spline, time=ytime, spline_fit=smoother, knot_points=kps)
  }
  
  ## covariance regression and sigma estimation 
  cc_BS0_feature <- blcc_fe(Xt=Xt, Wt=W, mats=mats, covreg=CR_algo, init=CR_init, 
                            max_iter=CR_maxiter, tol=CR_tol, S0=CR_S0init, B=CR_Binit, verb=CR_verb, term=CR_term) 
  
  cc_Sigma <- tSigma_optim(BS0=cc_BS0_feature, X=Xt)
  
  full_sigma_t <- sqrt_inv_build(cc_Sigma$tSigma)
  sqrt_sigma_t <- full_sigma_t$sqrt
  sqrt_inv_sigma_t <- full_sigma_t$inverse
  
  ## breve components 
  yb <- bln_YB(pY[['py']], sqrt_inv_sigma_t) 
  phib <- bln_phi_BT(phi_hat, P, sqrt_inv_sigma_t)
  hb <- lapply(seq_along(sqrt_inv_sigma_t), function(dt) {
    bln_Hb(H, P, sqrt_inv_sigma_t[[dt]])
  })
  
  ## check components 
  phic <- bln_phi_CT(phib, P) 
  yc <- PY_full(yb, t(P))
  hc <- lapply(hb, function(hb_dt) t(P) %*% hb_dt %*% P)
  
  ## curve cuts  
  ycut <- cf_curve_cut(yc[['py']], curves, mats)
  phicut <- cf_phi_cut(phic, curves, mats)
  
  ## normalized 
  normalized <- bln_normalized_fit(ycut, phicut, curves, mats)
  ystar_tenor <- blce_rescale_yield(normalized[['n_yield']], sqrt_sigma_t, P, "tenor", time=ytime, curves=curves, mats=mats)
  ystar_curve <- blce_rescale_yield(normalized$n_yield, sqrt_sigma_t, P, "curve", time=ytime, curves=curves, mats=mats)
  
  return(list(
    tenor=ystar_tenor, 
    curve=ystar_curve, 
    sigma_JT=cc_Sigma,
    BS0=cc_BS0_feature, 
    Xt=Xt))  
}