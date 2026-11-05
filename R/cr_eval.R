#############################################
#      Covariance Regression Testing        #
#############################################

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
    Sigma_array[, , i] <- .cr_covariance(S0, bx_i)
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
      X <- cbind(1, matrix(stats::rnorm(N * (Q - 1)), nrow = N, ncol = Q - 1))
    } else {
      X <- matrix(stats::rnorm(N * Q), nrow = N, ncol = Q)
    }
  } else {
    X <- as.matrix(X)
    if (ncol(X) != Q) stop("X must have ncol(X) = ncol(B).")
    if (nrow(X) != N) stop("X must have N rows.")
  }

  gamma <- stats::rnorm(N)

  # Cholesky draw for epsilon
  R <- chol(S0)
  eps <- matrix(stats::rnorm(N * P), nrow = N, ncol = P) %*% R

  BX <- X %*% Matrix::t(B)   # N x P ; row i is (B x_i)'
  W  <- gamma * BX + eps

  Sigma_true <- array(0, dim = c(P, P, N))
  for (i in seq_len(N)) {
    bxi <- matrix(BX[i, ], ncol = 1)
    Sigma_true[, , i] <- .cr_covariance(S0, bxi)
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
  Matrix::norm(A - B, type = "F")
}

cov_array_mse <- function(S_hat, S_true) {
  n <- dim(S_hat)[3]
  mean(vapply(seq_len(n), function(i) {
    mean((S_hat[, , i] - S_true[, , i])^2)
  }, numeric(1)))
}

sign_align_B <- function(B_hat, B_true) {
  err_plus  <- Matrix::norm(B_hat - B_true, type = "F")
  err_minus <- Matrix::norm(-B_hat - B_true, type = "F")
  if (err_plus <= err_minus) B_hat else -B_hat
}

is_psd <- function(S, tol = 1e-8) {
  ev <- eigen((S + Matrix::t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values
  min(ev) > -tol
}
