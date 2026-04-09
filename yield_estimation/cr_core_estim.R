#############################################
#     Covariance Regression Algorithms      #
#############################################

zero_mean_covreg <- function(W, X, init = c("adaptive", "static"), max_iter = 1000, tol = 1e-10, S0 = NULL, B = NULL, 
                             verb = FALSE, term = FALSE) {
  bnorm = 1e10; snorm = 1e10
  W <- as.matrix(W)
  P <- ncol(W); Q <- ncol(X); N <- nrow(W)
  
  if (is.null(S0)) { S0 <- cov(W) }
  if (is.null(B)) {
    if (init == "adaptive") {
      v0 <- rep(1, N)
      m0 <- rnorm(N)
      s0 <- rep(1, N) 
      X_tilde <- rbind(X*m0, X*s0) 
      W_tilde <-rbind(W, matrix(0, N, P))
      XtX <- t(X_tilde) %*% X_tilde
      B <- t(W_tilde) %*% X_tilde %*% solve(XtX)
    } 
    else {
      B <- matrix(rnorm(P * Q), nrow = P, ncol = Q)*1e-4
    }
  }
  
  b_est <- vector("list", length = max_iter) 
  s_est <- vector("list", length = max_iter)
  
  iS0 <- solve(S0) 
  iter <- 0 
  
  b_est[[1]] <- B; s_est[[1]] <- S0
  
  while(iter < max_iter) {
    iter <- iter + 1 

    ## E STEP 
    XtB <- X %*% t(B) 
    vz <- 1/(1 + apply((XtB %*% iS0) * XtB, 1, sum))
    mz <- vz * apply((XtB %*% iS0) * (W), 1, sum)
    sz <- sqrt(vz) 
    
    ## M STEP 
    X_tilde <- rbind(X*mz, X*sz) 
    W_tilde <-rbind(W, matrix(0, N, P))
    C <- t(W_tilde) %*% X_tilde %*% solve(t(X_tilde)%*%X_tilde)
    E <- W_tilde - (X_tilde %*% t(C))
    S0 <- (t(E) %*% E) / N
    iS0 <- solve(S0) 
    B <- C 
    
    ## STORE ESTIMATES 
    b_est[[iter+1]] <- B
    s_est[[iter+1]] <- S0 
    
    if (term){
      if (iter %% 10 == 0){
        bnorm <- norm(b_est[[iter]] - b_est[[iter-9]], type = 'F')
        snorm <- norm(s_est[[iter]] - s_est[[iter-9]], type = 'F')
      }
      if (bnorm < tol & snorm < tol) {return(list(S0=S0, B=B, b_est=b_est, s_est=s_est, iter=iter))}
    }
  }
  return(list(S0=S0, B=B, b_est=b_est, s_est=s_est))
}

dynamic_mean_covreg <- function() {
  ## TO WRITE 
}
