#############################################
#         Matrix Arithmetic Helpers         #
#############################################

symmPSD_sqrt <- function(sPSD) {
  ## compute square root matrix of a symmetric PSD matrix
  sPSD <- as.matrix(sPSD)
  if (!is.numeric(sPSD) || nrow(sPSD) == 0 || nrow(sPSD) != ncol(sPSD) || !all(is.finite(sPSD))) {
    stop("Covariance must be a finite, nonempty numeric square matrix.", call. = FALSE)
  }
  # Relative tolerances preserve behavior when covariance units change.
  if (max(abs(sPSD - Matrix::t(sPSD))) > 1e-8 * max(abs(sPSD))) {
    stop("Covariance matrix is materially asymmetric.", call. = FALSE)
  }
  sPSD <- 0.5 * sPSD + 0.5 * Matrix::t(sPSD)
  eig <- eigen(sPSD, symmetric = TRUE)
  lambda <- eig$values; U <- eig$vectors

  tolerance <- 1e-8 * max(abs(lambda))
  if (min(lambda) < -tolerance) {
    stop(sprintf("Covariance is not positive semidefinite: minimum eigenvalue %.6g; tolerance %.6g.",
                 min(lambda), tolerance), call. = FALSE)
  }
  lambda[lambda < 0] <- 0

  if (length(lambda) == 1L) return(matrix(sqrt(lambda), 1, 1))
  sqrt_lambda <- Matrix::diag(sqrt(lambda))
  sqrt_sPSD <- U %*% sqrt_lambda %*% Matrix::t(U)
  return(sqrt_sPSD)
}

sqrt_inv_build <- function(cov_matrix_ts) {
  sqrt_matrix_ts <- lapply(cov_matrix_ts, symmPSD_sqrt)
  inv_matrix_ts <- lapply(seq_along(sqrt_matrix_ts), function(i) {
    root <- sqrt_matrix_ts[[i]]
    # A square-root condition limit of eps^(-1/2) bounds covariance conditioning.
    if (rcond(root) <= sqrt(.Machine$double.eps)) {
      stop(sprintf("Covariance at list element %d is singular or numerically singular; whitening requires positive definiteness.", i),
           call. = FALSE)
    }
    Matrix::solve(root)
  })
  names(inv_matrix_ts) <- names(sqrt_matrix_ts)
  return(list(sqrt=sqrt_matrix_ts, inverse=inv_matrix_ts))
}
