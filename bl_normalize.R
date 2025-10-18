source("packages.R")
source("curve_reformat.R")

#############################################
#     Normalized Component Construction     #
#############################################

bln_build_P <- function(curves, mats) {
  ## build permutation matrix P swapping curve <-> tenor ordering
  M <- length(mats)
  D <- length(curves)
  
  total_size <- M * D
  original_indices <- seq_len(total_size) 
  target_indices <- as.vector(matrix(seq_len(total_size), nrow = M, byrow = TRUE))
  
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








