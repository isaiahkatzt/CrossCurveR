source("packages.R")
source("yield_estimation/curve_reformat.R")

#############################################
#     Normalized Component Construction     #
#############################################

bln_build_P <- function(curves, mats) {
  M <- length(mats)    
  D <- length(curves) 
  
  total_size <- M * D
  target_indices <- as.vector(
    matrix(seq_len(total_size), nrow = D, ncol = M, byrow = TRUE)
  )
  
  original_indices <- seq_len(total_size)
  
  P <- Matrix::sparseMatrix(
    i = target_indices,
    j = original_indices,
    x = 1,
    dims = c(total_size, total_size)
  )
  return(P)
}

bln_Hb <- function(HT, P, sqrt_inv_sigmaT) {
  ## single-day \breve{H} computation 
  tmp <- P %*% HT
  Hb <- sqrt_inv_sigmaT %*% tmp %*% t(P) %*% t(sqrt_inv_sigmaT)
  return(Hb)
}

bln_YB <- function(PY, sqrt_inv_sigmaT) {
  ## compute \breve{Y} across all days  
  N <- nrow(PY) 
  YB_matrix <- matrix(NA, nrow = N, ncol = ncol(PY))
  
  ## compute Sigma %*% t(PY[t, ]) <-> PY[t, ] %*% t(Sigma) 
  for (dt in seq_len(N)) {
    YB_matrix[dt, ] <- as.numeric(PY[dt, ]) %*% t(sqrt_inv_sigmaT[[dt]])
  }
  
  YB <- data.frame(YB_matrix)
  colnames(YB) <- colnames(PY) 
  return(YB)
}

bln_phi_hat <- function(curve_full_phi, mds) {
  ## stack phi_hat 
  curve_phi <- do.call(rbind, lapply(curve_full_phi, function(curve) t(curve[['phi']])))
  rownames(curve_phi) <- mds 
  return(curve_phi)
}

bln_phi_BT <- function(phi_hat, P, sqrt_inv_sigmaT) {
  ## compute \breve{Phi} across all days 
  P_Phi <- P %*% phi_hat 
  new_indices <- as.numeric(P %*% seq_len(nrow(phi_hat)))
  permuted_rn <- rownames(phi_hat)[new_indices] 
  
  lapply(sqrt_inv_sigmaT, function(sqrt_inv_sigma_t){
    phi_B_t <- sqrt_inv_sigma_t %*% P_Phi
    rownames(phi_B_t) <- permuted_rn
    phi_B_t
  })
}

bln_phi_CT <- function(phi_BT, P) {
  ## compute \check{Phi} across all days 
  tP <- t(P)
  
  n <- nrow(phi_BT[[1]])
  rn <- rownames(phi_BT[[1]])
  
  new_indices <- as.numeric(tP %*% seq_len(n))
  permuted_rn <- rn[new_indices]
  
  lapply(phi_BT, function(phi_Bt) {
    phi_Ct <- tP %*% phi_Bt
    rownames(phi_Ct) <- permuted_rn
    phi_Ct
  })
}

