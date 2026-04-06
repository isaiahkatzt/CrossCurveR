source("packages.R")
source("yield_estimation/curve_reformat.R")
source("yield_estimation/bl_single_curve.R")

#############################################
#     Diebold-Li-Yue Multi-Curve Model      #
#############################################

dly_phi <- function(lambda_optim, mats) {
  phi <- NS_loadings_ls(lambda_optim, mats) 
  cross_phi <- solve(tcrossprod(phi)) %*% phi 
  return(list(phi = phi, cross_phi = cross_phi))
}

dly_single_curve <- function(lambda_select, yields, mats) {
  return(setNames(list(
    curve_phi <- dly_phi(lambda_select, mats), 
    curve_betas <- bl_NSbetas(yields, lambda_select, curve_phi$cross_phi, mats),
    curve_nsfit <- bl_NSfit(curve_betas$m_yield, curve_phi$phi, curve_betas$NSbetas, lambda_select, mats)), c("phi", "betas", "nsfit")
  ))
}

dly_all_sc <- function(lambdas, yields, mats) {
  curve_lambdas <- vapply(yields, function(x) {lambda_grid_search(lambdas, x, mats, l2_loss_ls)}, numeric(1))
  curve_NSparam <- setNames(lapply(seq_along(yields), function(x) {
    dly_single_curve(curve_lambdas[x], yields[[x]], mats)
  }), curves) 
  return(curve_NSparam) 
}

dly_pull_beta <- function(dly_sclist, curves){
  tdly_sc <- transpose(dly_sclist) 
  lapply(tdly_sc[['betas']], function(x) x[['NSbetas']])
}

dly_group_factor <- function(sc_full, curves, latent=c("L", "S")) {
  latent = match.arg(latent, c("L", "S"))
  latent <- ifelse(latent=="L", 1, 2) 
  curve_beta <- dly_pull_beta(sc_full, curves) 
  factor_draw <- do.call(cbind, lapply(curve_beta, function(d) d[, latent]))
  return(factor_draw) 
}

