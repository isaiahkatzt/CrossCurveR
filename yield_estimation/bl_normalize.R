source("packages.R")
source("yield_estimation/curve_reformat.R")

#############################################
#     Normalized Component Construction     #
#############################################

bln_build_P <- function(curves, mats) {
  ## build curve <-> tenor permutation matrix  
  M <- length(mats)
  D <- length(curves)
  
  m <- rep(seq_len(M), each = D)
  d <- rep(seq_len(D), times = M)
  
  j <- (d - 1) * M + m # curve-major
  i <- (m - 1) * D + d # tenor-major
  
  Matrix::sparseMatrix(
    i = i,
    j = j,
    x = 1,
    dims = c(M * D, M * D)
  )
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

bln_rescale_beta_NS <- function(yield, check_phi, mats) {
  ## single curve rescaled estimates 
  N <- dim(yield)[1]
  mat_str <- paste0("X", lapply(mats, numeric_to_matname))
  matrix_yield <- as.matrix(yield) 
  daily_betas <- matrix(data = NA, nrow = N, ncol = 3) 
  NS_fit <- matrix(data = NA, nrow = N, ncol = length(mats))
  
  ## NS computation 
  for (i in 1:N) {
    crosscheck_phi <- solve(crossprod(check_phi[[i]])) %*% t(check_phi[[i]])
    daily_betas[i, ] <- matrix(t(crosscheck_phi %*% matrix_yield[i, ]))
    NS_fit[i, ] <- matrix(check_phi[[i]] %*% daily_betas[i, ])
  }
  
  colnames(NS_fit) <- mat_str
  return(list(
    yhat = NS_fit, 
    betahat = daily_betas))
}

bln_normalized_fit <- function(cut_yield, cut_phi, curves, mats) {
  stable_yield <- lapply(seq_along(cut_yield), function(curve) {
    bln_rescale_beta_NS(cut_yield[[curve]], cut_phi[[curve]], mats) 
  })
  names(stable_yield) <- curves 
  
  stable_ystar <- lapply(stable_yield, `[[`, "yhat") 
  stable_beta <- lapply(stable_yield, `[[`, "betahat")
  
  return(list(
    n_yield=stable_ystar, 
    n_beta=stable_beta)) 
}