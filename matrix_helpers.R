source("packages.R")
source("curve_reformat.R")

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

