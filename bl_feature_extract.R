source("packages.R")

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
  if (!trunc) {
    Xt <- do.call(cbind, cc_cspread) 
  } else {
    ## REVERSION RATE CONDITION 
  }
}

