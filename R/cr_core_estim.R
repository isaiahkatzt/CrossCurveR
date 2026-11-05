#############################################
#     Covariance Regression Algorithms      #
#############################################

.cr_validate_iteration_controls <- function(max_iter, tol, check_every) {
  counts <- list(max_iter = max_iter, check_every = check_every)
  for (name in names(counts)) {
    value <- counts[[name]]
    if (!(is.double(value) || is.integer(value)) || length(value) != 1L ||
        !is.finite(value) || value < 1 || value != floor(value)) {
      stop(sprintf("`%s` must be a finite positive integer (at least 1).", name), call. = FALSE)
    }
  }
  if (!(is.double(tol) || is.integer(tol)) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a finite positive numeric scalar.", call. = FALSE)
  }
}

zero_mean_covreg <- function(W, X, init = c("adaptive", "static"), max_iter = 1000, tol = 1e-10, S0 = NULL, B = NULL,
                             S0_shrink_diag = 0, verb = FALSE, term = FALSE) {
  init <- match.arg(init)
  bnorm = 1e10; snorm = 1e10
  W <- as.matrix(W)
  P <- ncol(W); Q <- ncol(X); N <- nrow(W)

  shrink_to_diagonal <- function(A, shrink) {
    if (!is.numeric(shrink) || length(shrink) != 1 || is.na(shrink) || shrink < 0 || shrink > 1) {
      stop("S0_shrink_diag must be a numeric scalar between 0 and 1.")
    }
    diag_A <- Matrix::diag(Matrix::diag(A))
    (1 - shrink) * A + shrink * diag_A
  }

  if (is.null(S0)) { S0 <- stats::cov(W) }
  S0 <- shrink_to_diagonal(0.5 * (S0 + Matrix::t(S0)), S0_shrink_diag)
  if (is.null(B)) {
    if (init == "adaptive") {
      v0 <- rep(1, N)
      m0 <- stats::rnorm(N)
      s0 <- rep(1, N)
      X_tilde <- rbind(X*m0, X*s0)
      W_tilde <-rbind(W, matrix(0, N, P))
      XtX <- Matrix::t(X_tilde) %*% X_tilde
      B <- Matrix::t(W_tilde) %*% X_tilde %*% Matrix::solve(XtX)
    }
    else {
      B <- matrix(stats::rnorm(P * Q), nrow = P, ncol = Q)*1e-4
    }
  }

  b_est <- vector("list", length = max_iter)
  s_est <- vector("list", length = max_iter)

  iS0 <- Matrix::solve(S0)
  iter <- 0

  b_est[[1]] <- B; s_est[[1]] <- S0

  while(iter < max_iter) {
    iter <- iter + 1

    ## E STEP
    XtB <- X %*% Matrix::t(B)
    vz <- 1/(1 + apply((XtB %*% iS0) * XtB, 1, sum))
    mz <- vz * apply((XtB %*% iS0) * (W), 1, sum)
    sz <- sqrt(vz)

    ## M STEP
    X_tilde <- rbind(X*mz, X*sz)
    W_tilde <-rbind(W, matrix(0, N, P))
    C <- Matrix::t(W_tilde) %*% X_tilde %*% Matrix::solve(Matrix::t(X_tilde)%*%X_tilde)
    E <- W_tilde - (X_tilde %*% Matrix::t(C))
    S0 <- shrink_to_diagonal((Matrix::t(E) %*% E) / N, S0_shrink_diag)
    iS0 <- Matrix::solve(S0)
    B <- C

    ## STORE ESTIMATES
    b_est[[iter+1]] <- B
    s_est[[iter+1]] <- S0

    if (term){
      if (iter %% 10 == 0){
        bnorm <- Matrix::norm(b_est[[iter]] - b_est[[iter-9]], type = 'F')
        snorm <- Matrix::norm(s_est[[iter]] - s_est[[iter-9]], type = 'F')
      }
      if (bnorm < tol & snorm < tol) {return(list(S0=S0, B=B, b_est=b_est, s_est=s_est, iter=iter, S0_shrink_diag=S0_shrink_diag))}
    }
  }
  return(list(S0=S0, B=B, b_est=b_est, s_est=s_est, S0_shrink_diag=S0_shrink_diag))
}

covreg_comp_estimator <- function(W, X, init = c("adaptive", "static"), max_iter = 500, tol = 1e-8,
                                  S0 = NULL, B = NULL, ridge_S0 = 1e-8, ridge_X = 1e-10,
                                  S0_shrink_diag = 0,
                                  check_every = 1, use_loglik = TRUE, verbose = FALSE,
                                  store_path = TRUE) {
  init <- match.arg(init)
  .cr_validate_iteration_controls(max_iter, tol, check_every)

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

  # Covariance and design matrices already include their intended ridge.
  safe_solve <- function(A) {
    cholA <- tryCatch(chol(A), error = function(e) NULL)
    if (is.null(cholA)) {
      Matrix::solve(A)
    } else {
      chol2inv(cholA)
    }
  }

  symmetrize <- function(A) {
    0.5 * (A + Matrix::t(A))
  }

  shrink_to_diagonal <- function(A, shrink) {
    if (!is.numeric(shrink) || length(shrink) != 1 || is.na(shrink) || shrink < 0 || shrink > 1) {
      stop("S0_shrink_diag must be a numeric scalar between 0 and 1.")
    }
    diag_A <- Matrix::diag(Matrix::diag(A))
    (1 - shrink) * A + shrink * diag_A
  }

  regularize_Psi <- function(A) {
    shrink_to_diagonal(symmetrize(A), S0_shrink_diag) + Matrix::diag(ridge_S0, nrow(A))
  }

  compute_loglik <- function(W, X, B, Psi) {
    iPsi <- safe_solve(Psi)
    logdet_Psi <- as.numeric(determinant(Psi, logarithm = TRUE)$modulus)

    ll <- 0
    for (i in seq_len(N)) {
      xi <- matrix(X[i, ], ncol = 1)
      wi <- matrix(W[i, ], ncol = 1)
      bxi <- B %*% xi

      quad_latent <- as.numeric(Matrix::t(bxi) %*% iPsi %*% bxi)
      logdet_Sigma_i <- logdet_Psi + log1p(quad_latent)
      iSigma_i <- iPsi - (iPsi %*% bxi %*% Matrix::t(bxi) %*% iPsi) / (1 + quad_latent)
      quad_i <- as.numeric(Matrix::t(wi) %*% iSigma_i %*% wi)

      ll <- ll - 0.5 * (P * log(2 * pi) + logdet_Sigma_i + quad_i)
    }
    ll
  }

  Psi <- if (is.null(S0)) stats::cov(W) else as.matrix(S0)
  Psi <- regularize_Psi(Psi)
  iPsi <- safe_solve(Psi)

  if (is.null(B)) {
    if (init == "static") {
      B <- matrix(stats::rnorm(P * Q, sd = 1e-4), nrow = P, ncol = Q)
    } else {
      m0 <- stats::rnorm(N)
      s0 <- rep(1, N)

      X_tilde <- rbind(X * m0, X * s0)
      W_tilde <- rbind(W, matrix(0, N, P))

      XtX <- Matrix::crossprod(X_tilde) + Matrix::diag(ridge_X, Q)
      B <- Matrix::t(W_tilde) %*% X_tilde %*% safe_solve(XtX)
    }
  } else {
    B <- as.matrix(B)
    if (!all(dim(B) == c(P, Q))) {
      stop("B must have dimensions ncol(W) x ncol(X).")
    }
  }

  if (store_path) {
    B_path <- vector("list", max_iter + 1)
    S0_path <- vector("list", max_iter + 1)
    ll_path <- rep(NA_real_, max_iter + 1)
    B_path[[1]] <- B
    S0_path[[1]] <- Psi
    if (use_loglik) {
      ll_path[1] <- compute_loglik(W, X, B, Psi)
    }
  } else {
    B_path <- S0_path <- ll_path <- NULL
  }

  converged <- FALSE
  iter <- 0

  while (iter < max_iter) {
    iter <- iter + 1
    B_old <- B

    XtB <- X %*% Matrix::t(B)
    quad <- Matrix::rowSums((XtB %*% iPsi) * XtB)
    vz <- 1 / (1 + quad)
    mz <- vz * Matrix::rowSums((XtB %*% iPsi) * W)
    sz <- sqrt(vz)

    X_tilde <- rbind(X * mz, X * sz)
    W_tilde <- rbind(W, matrix(0, N, P))

    XtX <- Matrix::crossprod(X_tilde) + Matrix::diag(ridge_X, Q)
    B <- Matrix::t(W_tilde) %*% X_tilde %*% safe_solve(XtX)

    idx <- which.max(abs(B))
    if (length(idx) == 1 && B[idx] < 0) {
      B <- -B
    }

    dB <- Matrix::norm(B - B_old, type = "F")
    ll_now <- NA_real_
    if (use_loglik) {
      ll_now <- compute_loglik(W, X, B, Psi)
    }

    if (store_path) {
      B_path[[iter + 1]] <- B
      S0_path[[iter + 1]] <- Psi
      ll_path[iter + 1] <- ll_now
    }

    if (verbose && (iter %% check_every == 0 || iter == 1)) {
      if (use_loglik) {
        cat(sprintf("iter = %d | dB = %.3e | loglik = %.8f\n", iter, dB, ll_now))
      } else {
        cat(sprintf("iter = %d | dB = %.3e\n", iter, dB))
      }
    }

    if (iter %% check_every == 0) {
      if (use_loglik && store_path && iter >= 1) {
        ll_prev <- ll_path[iter]
        if (is.finite(ll_now) && is.finite(ll_prev)) {
          if (abs(ll_now - ll_prev) < tol && dB < sqrt(tol)) {
            converged <- TRUE
            break
          }
        }
      } else {
        if (dB < tol) {
          converged <- TRUE
          break
        }
      }
    }
  }

  out <- list(
    B = B,
    S0 = Psi,
    Psi = Psi,
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
    out$Psi_path <- S0_path[seq_len(iter + 1)]
    if (use_loglik) {
      out$loglik_path <- ll_path[seq_len(iter + 1)]
    }
  }

  return(out)
}

dynamic_mean_covreg <- function() {
  ## TO WRITE
}
