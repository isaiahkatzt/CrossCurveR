#############################################
#      Covariance Regression Testing        #
#############################################

zero_mean_covreg_em <- function(W, X, init = c("adaptive", "static"), max_iter = 500, tol = 1e-8,
    S0 = NULL, B = NULL, ridge_S0 = 1e-8, ridge_X = 1e-10, S0_shrink_diag = 0,
    check_every = 1, use_loglik = TRUE, verbose = FALSE, store_path = TRUE) {
  
  init <- match.arg(init)
  
  W <- as.matrix(W)
  X <- as.matrix(X)
  
  if (nrow(W) != nrow(X)) {
    stop("W and X must have the same number of rows.")
  }
  
  N <- nrow(W)
  P <- ncol(W)
  Q <- ncol(X)
  
  if (N < 2) {
    stop("Need at least two observations.")
  }
  
  if (qr(X)$rank < Q) {
    stop("X is rank deficient. Remove collinear columns or regularize more heavily.")
  }
  
  safe_solve <- function(A, ridge = 0) {
    A2 <- A
    if (ridge > 0) {
      A2 <- A2 + diag(ridge, nrow(A2))
    }
    cholA <- tryCatch(chol(A2), error = function(e) NULL)
    if (is.null(cholA)) {
      solve(A2)
    } else {
      chol2inv(cholA)
    }
  }
  
  symmetrize <- function(A) {
    0.5 * (A + t(A))
  }

  shrink_to_diagonal <- function(A, shrink) {
    if (!is.numeric(shrink) || length(shrink) != 1 || is.na(shrink) || shrink < 0 || shrink > 1) {
      stop("S0_shrink_diag must be a numeric scalar between 0 and 1.")
    }
    diag_A <- diag(diag(A))
    (1 - shrink) * A + shrink * diag_A
  }

  regularize_S0 <- function(A) {
    shrink_to_diagonal(symmetrize(A), S0_shrink_diag) + diag(ridge_S0, nrow(A))
  }
  
  compute_loglik <- function(W, X, B, S0) {
    iS0 <- safe_solve(S0, ridge = ridge_S0)
    logdet_S0 <- as.numeric(determinant(S0, logarithm = TRUE)$modulus)
    
    ll <- 0
    for (i in seq_len(N)) {
      xi <- matrix(X[i, ], ncol = 1)
      wi <- matrix(W[i, ], ncol = 1)
      bxi <- B %*% xi
      
      quad_latent <- as.numeric(t(bxi) %*% iS0 %*% bxi)
      logdet_Sigma_i <- logdet_S0 + log1p(quad_latent)
      
      iSigma_i <- iS0 - (iS0 %*% bxi %*% t(bxi) %*% iS0) / (1 + quad_latent)
      quad_i <- as.numeric(t(wi) %*% iSigma_i %*% wi)
      
      ll <- ll - 0.5 * (P * log(2 * pi) + logdet_Sigma_i + quad_i)
    }
    ll
  }
  
  # initialization 
  if (is.null(S0)) {
    S0 <- cov(W)
  }
  S0 <- regularize_S0(S0)
  
  if (is.null(B)) {
    if (init == "static") {
      B <- matrix(rnorm(P * Q, sd = 1e-4), nrow = P, ncol = Q)
    } else {
      # adaptive start 
      m0 <- rnorm(N)
      s0 <- rep(1, N)
      
      X_tilde <- rbind(X * m0, X * s0)
      W_tilde <- rbind(W, matrix(0, N, P))
      
      XtX <- crossprod(X_tilde) + diag(ridge_X, Q)
      B <- t(W_tilde) %*% X_tilde %*% safe_solve(XtX)
    }
  } else {
    B <- as.matrix(B)
    if (!all(dim(B) == c(P, Q))) {
      stop("B must have dimensions ncol(W) x ncol(X).")
    }
  }
  
  # path storage 
  if (store_path) {
    B_path <- vector("list", max_iter + 1)
    S0_path <- vector("list", max_iter + 1)
    ll_path <- rep(NA_real_, max_iter + 1)
    B_path[[1]] <- B
    S0_path[[1]] <- S0
    if (use_loglik) {
      ll_path[1] <- compute_loglik(W, X, B, S0)
    }
  } else {
    B_path <- S0_path <- ll_path <- NULL
  }
  
  converged <- FALSE
  iter <- 0
  
  # EM algorithm  
  while (iter < max_iter) {
    iter <- iter + 1
    
    B_old <- B
    S0_old <- S0
    
    iS0 <- safe_solve(S0, ridge = ridge_S0)
    
    # E step 
    XtB <- X %*% t(B)  # N x P
    
    quad <- rowSums((XtB %*% iS0) * XtB)
    vz <- 1 / (1 + quad)
    mz <- vz * rowSums((XtB %*% iS0) * W)
    sz <- sqrt(vz)
    
    # M step 
    X_tilde <- rbind(X * mz, X * sz)          
    W_tilde <- rbind(W, matrix(0, N, P))        
    
    XtX <- crossprod(X_tilde) + diag(ridge_X, Q)
    B <- t(W_tilde) %*% X_tilde %*% safe_solve(XtX)
    
    E <- W_tilde - X_tilde %*% t(B)
    S0 <- regularize_S0(crossprod(E) / N)
    
    # sign normalization 
    idx <- which.max(abs(B))
    if (length(idx) == 1 && B[idx] < 0) {
      B <- -B
    }
    
    # diagnostic checking 
    dB <- norm(B - B_old, type = "F")
    dS <- norm(S0 - S0_old, type = "F")
    
    ll_now <- NA_real_
    if (use_loglik) {
      ll_now <- compute_loglik(W, X, B, S0)
    }
    
    if (store_path) {
      B_path[[iter + 1]] <- B
      S0_path[[iter + 1]] <- S0
      ll_path[iter + 1] <- ll_now
    }
    
    if (verbose && (iter %% check_every == 0 || iter == 1)) {
      if (use_loglik) {
        cat(sprintf(
          "iter = %d | dB = %.3e | dS0 = %.3e | loglik = %.8f\n",
          iter, dB, dS, ll_now
        ))
      } else {
        cat(sprintf(
          "iter = %d | dB = %.3e | dS0 = %.3e\n",
          iter, dB, dS
        ))
      }
    }
    
    # convergence check 
    if (iter %% check_every == 0) {
      if (use_loglik && store_path && iter >= 1) {
        ll_prev <- ll_path[iter]
        if (is.finite(ll_now) && is.finite(ll_prev)) {
          if (abs(ll_now - ll_prev) < tol && dB < sqrt(tol) && dS < sqrt(tol)) {
            converged <- TRUE
            break
          }
        }
      } else {
        if (dB < tol && dS < tol) {
          converged <- TRUE
          break
        }
      }
    }
  }
  
  out <- list(
    B = B,
    S0 = S0,
    mz = mz,
    vz = vz,
    sz = sz,
    iter = iter,
    converged = converged,
    loglik = if (use_loglik) ll_now else NULL,
    S0_shrink_diag = S0_shrink_diag
  )
  
  if (store_path) {
    out$B_path <- B_path[seq_len(iter + 1)]
    out$S0_path <- S0_path[seq_len(iter + 1)]
    if (use_loglik) {
      out$loglik_path <- ll_path[seq_len(iter + 1)]
    }
  }
  
  return(out) 
}

cr_build <- function(object, X_new = NULL, return_array = TRUE, drop_single = TRUE) {
  # Accept either the fitted object from zero_mean_covreg_em()
  # or a list containing $B and $S0
  if (is.null(object$B) || is.null(object$S0)) {
    stop("object must contain components $B and $S0.")
  }
  
  B  <- as.matrix(object$B)
  S0 <- as.matrix(object$S0)
  
  P <- nrow(B)
  Q <- ncol(B)
  
  if (!all(dim(S0) == c(P, P))) {
    stop("Dimensions of S0 are incompatible with B.")
  }
  
  # If no new design matrix is supplied, try to use the original X if stored
  if (is.null(X_new)) {
    if (!is.null(object$X)) {
      X_new <- object$X
    } else {
      stop("X_new is NULL and object does not contain a stored design matrix X.")
    }
  }
  
  X_new <- as.matrix(X_new)
  
  if (ncol(X_new) != Q) {
    stop("X_new must have ncol(X_new) = ncol(B).")
  }
  
  N_new <- nrow(X_new)
  
  # 3D array: P x P x N_new
  Sigma_array <- array(0, dim = c(P, P, N_new))
  
  for (i in seq_len(N_new)) {
    x_i <- matrix(X_new[i, ], ncol = 1)
    bx_i <- B %*% x_i
    Sigma_array[, , i] <- S0 + bx_i %*% t(bx_i)
  }
  
  if (return_array) {
    if (N_new == 1 && drop_single) {
      return(Sigma_array[, , 1, drop = FALSE])
    }
    return(Sigma_array)
  }
  
  # Alternative: return a list of covariance matrices
  Sigma_list <- lapply(seq_len(N_new), function(i) Sigma_array[, , i, drop = FALSE])
  
  if (N_new == 1 && drop_single) {
    return(Sigma_list[[1]])
  }
  
  return(Sigma_list)
}

#############################################
#         Synthetic CR Data Checker         #
#############################################

cr_synthetic <- function(N = 500, B, S0, X = NULL, intercept = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  B  <- as.matrix(B)
  S0 <- as.matrix(S0)
  
  P <- nrow(B)
  Q <- ncol(B)
  
  if (!all(dim(S0) == c(P, P))) {
    stop("S0 must be a P x P matrix where P = nrow(B).")
  }
  
  if (is.null(X)) {
    if (intercept) {
      if (Q < 1) stop("Q must be at least 1 when intercept = TRUE.")
      X <- cbind(1, matrix(rnorm(N * (Q - 1)), nrow = N, ncol = Q - 1))
    } else {
      X <- matrix(rnorm(N * Q), nrow = N, ncol = Q)
    }
  } else {
    X <- as.matrix(X)
    if (ncol(X) != Q) stop("X must have ncol(X) = ncol(B).")
    if (nrow(X) != N) stop("X must have N rows.")
  }
  
  gamma <- rnorm(N)
  
  # Cholesky draw for epsilon
  R <- chol(S0)
  eps <- matrix(rnorm(N * P), nrow = N, ncol = P) %*% R
  
  BX <- X %*% t(B)   # N x P ; row i is (B x_i)'
  W  <- gamma * BX + eps
  
  Sigma_true <- array(0, dim = c(P, P, N))
  for (i in seq_len(N)) {
    bxi <- matrix(BX[i, ], ncol = 1)
    Sigma_true[, , i] <- S0 + bxi %*% t(bxi)
  }
  
  list(
    W = W,
    X = X,
    gamma = gamma,
    eps = eps,
    B_true = B,
    S0_true = S0,
    Sigma_true = Sigma_true
  )
}

#############################################
#            Evaluation Helpers             #
#############################################

frobenius_error <- function(A, B) {
  norm(A - B, type = "F")
}

cov_array_mse <- function(S_hat, S_true) {
  n <- dim(S_hat)[3]
  mean(vapply(seq_len(n), function(i) {
    mean((S_hat[, , i] - S_true[, , i])^2)
  }, numeric(1)))
}

sign_align_B <- function(B_hat, B_true) {
  err_plus  <- norm(B_hat - B_true, type = "F")
  err_minus <- norm(-B_hat - B_true, type = "F")
  if (err_plus <= err_minus) B_hat else -B_hat
}

is_psd <- function(S, tol = 1e-8) {
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values
  min(ev) > -tol
}
