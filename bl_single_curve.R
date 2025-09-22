source("packages.R")

#############################################
#      Baseline Single Curve Estimation     #
#############################################

l2_loss_lambda <- function(lambda_grid, daily_yield, mats) {
  n <- length(lambda_grid)
  loss_row <- numeric(n)  
  factor_loadings <- lapply(lambda_grid, NS_loadings, mats)
  
  for (i in seq_along(lambda_grid)) {
    X <- t(factor_loadings[[i]])  
    fit <- lm.fit(X, daily_yield)  
    loss_row[i] <- sum(fit$residuals^2)  
  }
  loss_row
}

lambda_grid_search <- function(lambda_grid, yields, mats) {  
  ## fast search lambda grid 
  yields <- as.matrix(yields[, -1, drop = FALSE])  # Ensure matrix structure
  N <- nrow(yields)
  n <- length(lambda_grid)
  
  ## vectorized loss computation 
  loss_matrix <- t(vapply(seq_len(N), function(i) l2_loss_lambda(lambda_grid, yields[i, ], mats), numeric(n)))
  
  ## optimal lambda selection; varies across curves! 
  optimal_lambda <- lambda_grid[which.min(colSums(loss_matrix))]
  return(optimal_lambda)
}

bl_phi <- function(lambda_optim, mats) {
  phi <- NS_loadings(lambda_optim, mats) 
  cross_phi <- solve(tcrossprod(phi)) %*% phi 
  return(list(phi = phi, cross_phi = cross_phi))
}

bl_NSbetas <- function(yield, lambda_optim, cross_phi, mats) {
  matrix_yield <- as.matrix(yield[-1]) 
  NSbetas <- t(cross_phi %*% t(matrix_yield))
  return(list(NSbetas = NSbetas, m_yield = matrix_yield))
}

bl_NSfit <- function(matrix_yield, phi, NSbetas, lambda_optim, mats) {
  NSfit <- NSbetas %*% phi 
  residual <- NSfit - matrix_yield 
  return(list(NSfit = NSfit, W = residual))
}

bl_bdiag_H <- function(Wt, rate_cutoffs, mats) {
  # Validate cutoffs
  if (!is.numeric(rate_cutoffs) || length(rate_cutoffs) != 2) {
    stop("Invalid rate cutoffs (requires a numeric vector of length 2)")
  }
  r1 <- rate_cutoffs[1]
  r2 <- rate_cutoffs[2]
  p <- length(mats)
  
  if (r1 < 1 || r2 <= r1 || r2 >= p) {
    stop("Rate cutoffs must be valid indices within the maturity range")
  }
  
  # cutoff slicing 
  short_idx <- 1:r1
  med_idx   <- (r1 + 1):r2
  long_idx  <- (r2 + 1):p
  sm_idx    <- r1:(r1 + 1)
  ml_idx    <- r2:(r2 + 1)
  
  # extract rate-specific views
  short_rates <- Wt[, short_idx, drop = FALSE]
  med_rates   <- Wt[, med_idx, drop = FALSE]
  long_rates  <- Wt[, long_idx, drop = FALSE]
  coupling_SM <- Wt[, sm_idx, drop = FALSE]
  coupling_ML <- Wt[, ml_idx, drop = FALSE]
  
  # preallocate zero matrix and insert covariance blocks
  H_matrix <- matrix(0, nrow = p, ncol = p)
  
  H_matrix[short_idx, short_idx] <- cov(short_rates)
  H_matrix[med_idx,   med_idx]   <- cov(med_rates)
  H_matrix[long_idx,  long_idx]  <- cov(long_rates)
  
  sigma_SM <- cov(coupling_SM)
  sigma_ML <- cov(coupling_ML)
  
  # insert coupling terms 
  H_matrix[sm_idx[1], sm_idx[2]] <- sigma_SM[1, 2]
  H_matrix[sm_idx[2], sm_idx[1]] <- sigma_SM[1, 2]
  H_matrix[ml_idx[1], ml_idx[2]] <- sigma_ML[1, 2]
  H_matrix[ml_idx[2], ml_idx[1]] <- sigma_ML[1, 2]
  
  return(H_matrix)
}

bl_single_curve <- function(lambda_select, yields, cutoffs, mats) {
  return(setNames(list(
    curve_phi <- bl_phi(lambda_select, mats), 
    curve_betas <- bl_NSbetas(yields, lambda_select, curve_phi$cross_phi, mats), 
    curve_nsfit <- bl_NSfit(curve_betas$m_yield, curve_phi$phi, curve_betas$NSbetas, lambda_select, mats), 
    curve_H <- bl_bdiag_H(curve_nsfit$W, cutoffs, mats)), c("phi", "betas", "nsfit", "H")
  ))
}
