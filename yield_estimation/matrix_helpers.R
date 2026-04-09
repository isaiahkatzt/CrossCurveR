#############################################
#         Matrix Arithmetic Helpers         #
#############################################

symmPSD_sqrt <- function(sPSD) {
  ## compute square root matrix of a symmetric PSD matrix 
  eig <- eigen(as.matrix(sPSD), symmetric = TRUE)
  lambda <- eig$values; U <- eig$vectors
  
  ## negative eigenvalue thresholding 
  lambda[lambda < 0] <- 0 
  
  sqrt_lambda <- diag(sqrt(lambda)) 
  sqrt_sPSD <- U %*% sqrt_lambda %*% t(U) 
  return(sqrt_sPSD)
}

sqrt_inv_build <- function(cov_matrix_ts) {
  sqrt_matrix_ts <- lapply(cov_matrix_ts, symmPSD_sqrt) 
  inv_matrix_ts <- lapply(sqrt_matrix_ts, function(matrix) solve(matrix)) 
  return(list(sqrt=sqrt_matrix_ts, inverse=inv_matrix_ts))
}
