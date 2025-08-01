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

cc_VAR_optim <- function(cDNS, reference, curves) {
  ## combine reference curve and marginal curve DNS parameters as VARs  
  marginals <- setdiff(curves, reference)
  reference_beta <- cDNS[[reference]][[2]]
  marginal_betas <- lapply(marginals, function(m) cDNS[[m]][[2]])
  
  cross_curve_beta <- vector("list", length(marginals))
  names(cross_curve_beta) <- marginals
  
  for (i in seq_along(marginals)) {
    m <- marginals[i]
    cc_beta <- cbind(reference_beta, marginal_betas[[i]])
    colnames(cc_beta) <- c(
      paste0(reference, ".L"), paste0(reference, ".S"), paste0(reference, ".C"),
      paste0(m, ".L"), paste0(m, ".S"), paste0(m, ".C")
    )
    cross_curve_beta[[i]] <- cc_beta
  }
  return(cross_curve_beta)
}
