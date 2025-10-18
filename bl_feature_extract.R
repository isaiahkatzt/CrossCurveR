source("packages.R")
source("curve_reformat.R")

 #############################################
#            Feature Extraction             #
#############################################

blmc_phi <- function(bl_sc_list){
  return(t(do.call(cbind, lapply(tsc_list$phi, `[[`, "phi"))))
}

blmc_VAR <- function(blsc_list, reference, curves) {
  marginals <- setdiff(curves, reference) 
  reference_beta <- blsc_list$betas[[reference]]
  cross_curve_beta <- setNames(vector("list", length(marginals)), marginals)
  
  for (m in marginals) {
    m_beta <- blsc_list$betas[[m]]
    cc_beta <- cbind(reference_beta$NSbetas, m_beta$NSbetas)
    colnames(cc_beta) <- c(
      paste0(reference, c(".L", ".S", ".C")),
      paste0(m, c(".L", ".S", ".C"))
    )
    cross_curve_beta[[m]] <- cc_beta
  }
  return(cross_curve_beta)
}

blmc_ECM <- function(cc_betas, estim = "ML", type = "eigen", alpha = 0.1) {
  lag_opt <- max(VARselect(cc_betas, type = "const")$selection["SC(n)"], 2)
  vecm_unrestricted <- VECM(cc_betas, lag = lag_opt - 1, estim = estim, include = "const")
  r_test <- rank.test(vecm_unrestricted, type = type, cval = alpha)
  r <- r_test$r
  
  ## rank zero condition 
  if (r == 0) {
    return(list(Gmatrix = 0, Fmatrix = 0, rank = 0))
  }
  
  vecm_restricted <- VECM(cc_betas, lag = lag_opt - 1, r = r, estim = estim, include = "const")
  
  ## round for floating point precision 
  Fmat <- round(coefA(vecm_restricted), 14)
  Gmat <- round(coefB(vecm_restricted), 14)
  return(list(Fmatrix = Fmat, Gmatrix = Gmat, rank = r))
}

blmc_cspread <- function(cc_betas, Gmatrix, normalize = TRUE) {
  Xt <- cc_betas %*% Gmatrix 
  if (normalize) {
    Xt <- sweep(Xt, 2, colMeans(Xt), FUN = "-")
    Xt <- sweep(Xt, 2, apply(Xt, 2, sd), FUN = "/")
  }
  return(Xt)
}

blmc_Fselect <- function(Fmatrix) {
  reversion_norm <- sqrt(colSums(Fmatrix^2))
  
}

blcc_cointegration <- function(blsc_list, curves, reference, estim = "ML", type = "eigen", alpha = 0.1, normalize = TRUE) {
  ## cross curve beta list 
  cc_betas <- blmc_VAR(blsc_list, reference, curves) 
  
  ## curve-specific error correct model construction 
  full_ecm <- lapply(cc_betas, function(cc_var) blmc_ECM(cc_betas = cc_var, estim = estim, type = type, alpha = alpha))
  
  ## cointegration spread construction 
  full_cspread <- setNames(lapply(seq_along(cc_betas), function(curve) {
    if (full_ecm[[curve]]$rank != 0) {
      blmc_cspread(cc_betas = cc_betas[[curve]], Gmatrix = full_ecm[[curve]]$Gmatrix, normalize = normalize)
    }
  }), names(cc_betas))
  
  return(list(
    cc_cspread = full_cspread, 
    cc_ecm = full_ecm, 
    cc_beta = cc_betas
  ))
}

blcc_build_Xt <- function(cc_cspread, cc_ecm, trunc = FALSE) {
  ## construct full Xt covariate matrix 
  if (!trunc) {
    Xt <- do.call(cbind, cc_cspread) 
  } else {
  ## placeholder condition for testing  
    Xt <- do.call(cbind, lapply(cc_cspread, function(x) x[, 1]))
  }
  return(Xt) 
}

blcc_build_Wj <- function(curve_nsfit, mats){
  ## construct tenor-specific Wj list 
  wc_list <- lapply(curve_nsfit, `[[`, 1) 
  mat_names <- paste0("X", sapply(mats, numeric_to_matname))
  wj_full <- setNames(lapply(seq_along(mats), function(i) {
    do.call(cbind, lapply(wc_list, function(x) x[, i])) 
  }), mat_names)
  
  return(Wj = wj_full)
}

blcc_fe <- function(Xt, Wt, mats, covreg=zero_mean_covreg,
                    init = "adaptive", max_iter = 1000, tol = 1e-10, S0 = NULL, B = NULL, verb = FALSE, term = FALSE) {
  ## baseline cross-curve feature extraction
  BS0_list <- lapply(Wt, function(Wj) {
    covreg(W=Wj, X=Xt, init = init, max_iter = max_iter, tol = tol, S0 = S0, B = B, 
           verb = verb, term = term)
  })
  return(BS0_list) 
}

sigma_jt_optim <- function(BS0j, x) {
  ## compute single-day sigma_{jt}   
  B <- BS0j$B
  S0 <- BS0j$S0
  Bx <- B %*% x
  Sigma <- S0 + tcrossprod(Bx)  
  return(Sigma)
}

jSigma_optim <- function(BS0j, X) {
  # compute full period sigma_{jt}
  N <- nrow(X)
  day_list <- vector("list", N)
  
  jSigma <- lapply(seq_len(N), function(i) {
    sigma_jt_optim(BS0j = BS0j, x = X[i, ])
  })
  return(jSigma)
}

tSigma_optim <- function(BS0, X) {
  # compute full period, all sigma_{jt} 
  N <- nrow(X) 
  sigma_j_full <- lapply(BS0, function(j) jSigma_optim(j, X))
  sigma_t_full <- lapply(seq_len(N), function(t) {
    bdiag(lapply(sigma_j_full, function(j) j[[t]]))
  })
  return(list(tSigma = sigma_t_full, jSigma = sigma_j_full))
}

