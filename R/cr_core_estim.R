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

.cr_warn_nonconvergence <- function(estimator, iter, max_iter, tol, dimensions,
                                    changes, check_every, last_check, criterion,
                                    loglik = NULL, loglik_change = NULL, note = NULL) {
  details <- c(
    sprintf("iter=%d, tol=%.3e", iter, tol),
    paste(sprintf("%s=%d", names(dimensions), dimensions), collapse = ", "),
    sprintf("check_every=%d, last_check=%d", check_every, last_check),
    sprintf("criterion=%s", criterion),
    sprintf("%s=%.3e", names(changes), changes)
  )
  if (!is.null(loglik)) details <- c(details, sprintf("final_score=%.6g", loglik))
  if (!is.null(loglik_change)) {
    details <- c(details, sprintf("absolute_score_change=%.3e", loglik_change))
    if (!is.finite(loglik_change)) details <- c(details, "likelihood score change is nonfinite")
  }
  if (!is.null(loglik) && !is.finite(loglik)) {
    details <- c(details, "likelihood score is nonfinite")
  }
  if (last_check == 0) {
    details <- c(details, "no convergence check was performed")
  } else if (last_check < iter) {
    details <- c(details, "iteration limit fell between convergence checks")
  }
  details <- c(details, note)
  message <- sprintf(
    "%s reached max_iter=%d without confirmed convergence (%s). Returning final estimates.",
    estimator, max_iter, paste(details, collapse = "; ")
  )
  warning(structure(list(
    message = message, call = NULL, estimator = estimator,
    iter = iter, max_iter = max_iter, tol = tol, dimensions = dimensions,
    parameter_changes = changes, check_every = check_every, last_check = last_check,
    criterion = criterion, loglik = loglik, loglik_change = loglik_change,
    note = note
  ), class = c("crosscurver_nonconvergence", "warning", "condition")))
}

# Shared numerical operations.
.cr_inverse <- function(A, ridge = 0) {
  A2 <- A
  if (ridge > 0) {
    A2 <- A2 + Matrix::diag(ridge, nrow(A2))
  }
  cholA <- tryCatch(chol(A2), error = function(e) NULL)
  if (is.null(cholA)) {
    Matrix::solve(A2)
  } else {
    chol2inv(cholA)
  }
}

.cr_shrink_diag <- function(A, shrink) {
  if (!is.numeric(shrink) || length(shrink) != 1 || is.na(shrink) || shrink < 0 || shrink > 1) {
    stop("S0_shrink_diag must be a numeric scalar between 0 and 1.")
  }
  diag_A <- Matrix::diag(Matrix::diag(A))
  (1 - shrink) * A + shrink * diag_A
}

.cr_regularize <- function(A, shrink, ridge) {
  .cr_shrink_diag(0.5 * (A + Matrix::t(A)), shrink) + Matrix::diag(ridge, nrow(A))
}

.cr_covariance <- function(S0, Bx) {
  S0 + Matrix::tcrossprod(Bx)
}

.cr_loglik <- function(W, X, B, S0, ridge_S0 = 0) {
  N <- nrow(W)
  P <- ncol(W)
  # Preserve the historical score: the inverse has an additional covariance
  # ridge, while the determinant uses stored S0. Changing this alters fitting.
  iS0 <- .cr_inverse(S0, ridge = ridge_S0)
  logdet_S0 <- as.numeric(determinant(S0, logarithm = TRUE)$modulus)

  ll <- 0
  for (i in seq_len(N)) {
    xi <- matrix(X[i, ], ncol = 1)
    wi <- matrix(W[i, ], ncol = 1)
    bxi <- B %*% xi

    quad_latent <- as.numeric(Matrix::t(bxi) %*% iS0 %*% bxi)
    logdet_Sigma_i <- logdet_S0 + log1p(quad_latent)

    iSigma_i <- iS0 - (iS0 %*% bxi %*% Matrix::t(bxi) %*% iS0) / (1 + quad_latent)
    quad_i <- as.numeric(Matrix::t(wi) %*% iSigma_i %*% wi)

    ll <- ll - 0.5 * (P * log(2 * pi) + logdet_Sigma_i + quad_i)
  }
  ll
}

zero_mean_covreg_em <- function(W, X, init = c("adaptive", "static"), max_iter = 500, tol = 1e-8,
    S0 = NULL, B = NULL, ridge_S0 = 1e-8, ridge_X = 1e-10, S0_shrink_diag = 0,
    check_every = 1, use_loglik = TRUE, verbose = FALSE, store_path = TRUE) {

  init <- match.arg(init)
  .cr_validate_iteration_controls(max_iter, tol, check_every)

  inputs <- .cc_cov_inputs(W, X, S0, B)
  W <- inputs$W; X <- inputs$X; S0 <- inputs$S0; B <- inputs$B

  if (nrow(W) != nrow(X)) {
    stop("W and X must have the same number of rows.")
  }

  N <- nrow(W)
  P <- ncol(W)
  Q <- ncol(X)

  if (N < 2) {
    stop("Covariance regression: at least two observations are required; received ", N, ".", call. = FALSE)
  }

  if (qr(X)$rank < Q) {
    stop("X is rank deficient under the current numerical tolerance. ",
         "This estimator requires full-column-rank X before applying ridge stabilization. ",
         "Remove or explicitly reparameterize redundant columns; increasing `ridge_X` will not bypass this check.")
  }

  # initialization
  if (is.null(S0)) {
    S0 <- stats::cov(W)
  }
  S0 <- .cr_regularize(S0, S0_shrink_diag, ridge_S0)

  if (is.null(B)) {
    if (init == "static") {
      B <- matrix(stats::rnorm(P * Q, sd = 1e-4), nrow = P, ncol = Q)
    } else {
      # adaptive start
      m0 <- stats::rnorm(N)
      s0 <- rep(1, N)

      X_tilde <- rbind(X * m0, X * s0)
      W_tilde <- rbind(W, matrix(0, N, P))

      XtX <- Matrix::crossprod(X_tilde) + Matrix::diag(ridge_X, Q)
      B <- Matrix::t(W_tilde) %*% X_tilde %*% .cr_inverse(XtX)
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
      ll_path[1] <- .cr_loglik(W, X, B, S0, ridge_S0)
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

    # Preserve the pre-correction inversion ridge in addition to stored S0's ridge.
    iS0 <- .cr_inverse(S0, ridge = ridge_S0)

    # E step
    XtB <- X %*% Matrix::t(B)  # N x P

    quad <- Matrix::rowSums((XtB %*% iS0) * XtB)
    vz <- 1 / (1 + quad)
    mz <- vz * Matrix::rowSums((XtB %*% iS0) * W)
    sz <- sqrt(vz)

    # M step
    X_tilde <- rbind(X * mz, X * sz)
    W_tilde <- rbind(W, matrix(0, N, P))

    XtX <- Matrix::crossprod(X_tilde) + Matrix::diag(ridge_X, Q)
    B <- Matrix::t(W_tilde) %*% X_tilde %*% .cr_inverse(XtX)

    E <- W_tilde - X_tilde %*% Matrix::t(B)
    S0 <- .cr_regularize(Matrix::crossprod(E) / N, S0_shrink_diag, ridge_S0)

    # sign normalization
    idx <- which.max(abs(B))
    if (length(idx) == 1 && B[idx] < 0) {
      B <- -B
    }

    # diagnostic checking
    dB <- Matrix::norm(B - B_old, type = "F")
    dS <- Matrix::norm(S0 - S0_old, type = "F")

    ll_now <- NA_real_
    if (use_loglik) {
      ll_now <- .cr_loglik(W, X, B, S0, ridge_S0)
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
          if (abs(ll_now - ll_prev) < tol && dB < tol && dS < tol) {
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
    S0_shrink_diag = S0_shrink_diag,
    input_spec = .cc_cov_spec(W, X)
  )

  if (store_path) {
    out$B_path <- B_path[seq_len(iter + 1)]
    out$S0_path <- S0_path[seq_len(iter + 1)]
    if (use_loglik) {
      out$loglik_path <- ll_path[seq_len(iter + 1)]
    }
  }

  if (!converged) {
    .cr_warn_nonconvergence(
      "zero_mean_covreg_em", iter, max_iter, tol, c(N = N, P = P, Q = Q),
      c(dB = dB, dS0 = dS), check_every, (iter %/% check_every) * check_every,
      if (use_loglik && store_path) "parameter and likelihood changes" else "parameter changes",
      loglik = if (use_loglik) ll_now else NULL,
      loglik_change = if (use_loglik && store_path) abs(ll_now - ll_path[iter]) else NULL
    )
  }

  return(out)
}

zero_mean_covreg <- function(W, X, init = c("adaptive", "static"), max_iter = 1000, tol = 1e-10, S0 = NULL, B = NULL,
                             S0_shrink_diag = 0, verb = FALSE, term = FALSE) {
  init <- match.arg(init)
  bnorm = 1e10; snorm = 1e10
  inputs <- .cc_cov_inputs(W, X, S0, B)
  W <- inputs$W; X <- inputs$X; S0 <- inputs$S0; B <- inputs$B
  P <- ncol(W); Q <- ncol(X); N <- nrow(W)

  if (is.null(S0)) { S0 <- stats::cov(W) }
  S0 <- .cr_shrink_diag(0.5 * (S0 + Matrix::t(S0)), S0_shrink_diag)
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
  converged <- FALSE
  mz <- vz <- sz <- NULL

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
    S0 <- .cr_shrink_diag((Matrix::t(E) %*% E) / N, S0_shrink_diag)
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
      if (bnorm < tol & snorm < tol) {
        converged <- TRUE
        break
      }
    }
  }
  if (!converged && iter > 0 && iter >= max_iter) {
    last_check <- if (term) (iter %/% 10) * 10 else 0
    changes <- if (term) {
      if (last_check > 0) c(dB = bnorm, dS0 = snorm) else numeric()
    } else {
      c(dB = Matrix::norm(B - b_est[[iter]], "F"),
        dS0 = Matrix::norm(S0 - s_est[[iter]], "F"))
    }
    .cr_warn_nonconvergence(
      "zero_mean_covreg", iter, max_iter, tol, c(N = N, P = P, Q = Q),
      changes, 10, last_check,
      if (term) "legacy history comparison" else "not assessed (term=FALSE)",
      note = if (term && last_check > 0) {
        "dB and dS0 are from the last scheduled legacy history comparison"
      } else NULL
    )
  }
  list(
    B = B, S0 = S0, mz = mz, vz = vz, sz = sz,
    iter = iter, converged = converged, loglik = NULL,
    S0_shrink_diag = S0_shrink_diag, input_spec = .cc_cov_spec(W, X),
    B_path = b_est[seq_len(iter + 1)], S0_path = s_est[seq_len(iter + 1)],
    b_est = b_est, s_est = s_est,
    termination_reason = if (converged) "converged" else if (!term) "not_assessed" else "max_iter"
  )
}

covreg_comp_estimator <- function(W, X, init = c("adaptive", "static"), max_iter = 500, tol = 1e-8,
                                  S0 = NULL, B = NULL, ridge_S0 = 1e-8, ridge_X = 1e-10,
                                  S0_shrink_diag = 0,
                                  check_every = 1, use_loglik = TRUE, verbose = FALSE,
                                  store_path = TRUE) {
  init <- match.arg(init)
  .cr_validate_iteration_controls(max_iter, tol, check_every)

  inputs <- .cc_cov_inputs(W, X, S0, B)
  W <- inputs$W; X <- inputs$X; S0 <- inputs$S0; B <- inputs$B

  if (nrow(W) != nrow(X)) {
    stop("W and X must have the same number of rows.")
  }

  N <- nrow(W)
  P <- ncol(W)
  Q <- ncol(X)

  if (N < 2) {
    stop("Covariance regression: at least two observations are required; received ", N, ".", call. = FALSE)
  }

  if (qr(X)$rank < Q) {
    stop("X is rank deficient under the current numerical tolerance. ",
         "This estimator requires full-column-rank X before applying ridge stabilization. ",
         "Remove or explicitly reparameterize redundant columns; increasing `ridge_X` will not bypass this check.")
  }

  Psi <- if (is.null(S0)) stats::cov(W) else as.matrix(S0)
  Psi <- .cr_regularize(Psi, S0_shrink_diag, ridge_S0)
  iPsi <- .cr_inverse(Psi, ridge = ridge_S0)

  if (is.null(B)) {
    if (init == "static") {
      B <- matrix(stats::rnorm(P * Q, sd = 1e-4), nrow = P, ncol = Q)
    } else {
      m0 <- stats::rnorm(N)
      s0 <- rep(1, N)

      X_tilde <- rbind(X * m0, X * s0)
      W_tilde <- rbind(W, matrix(0, N, P))

      XtX <- Matrix::crossprod(X_tilde) + Matrix::diag(ridge_X, Q)
      B <- Matrix::t(W_tilde) %*% X_tilde %*% .cr_inverse(XtX)
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
      ll_path[1] <- .cr_loglik(W, X, B, Psi, ridge_S0)
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
    B <- Matrix::t(W_tilde) %*% X_tilde %*% .cr_inverse(XtX)

    idx <- which.max(abs(B))
    if (length(idx) == 1 && B[idx] < 0) {
      B <- -B
    }

    dB <- Matrix::norm(B - B_old, type = "F")
    ll_now <- NA_real_
    if (use_loglik) {
      ll_now <- .cr_loglik(W, X, B, Psi, ridge_S0)
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
          if (abs(ll_now - ll_prev) < tol && dB < tol) {
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
    S0_shrink_diag = S0_shrink_diag,
    input_spec = .cc_cov_spec(W, X)
  )

  if (store_path) {
    out$B_path <- B_path[seq_len(iter + 1)]
    out$S0_path <- S0_path[seq_len(iter + 1)]
    out$Psi_path <- S0_path[seq_len(iter + 1)]
    if (use_loglik) {
      out$loglik_path <- ll_path[seq_len(iter + 1)]
    }
  }

  if (!converged) {
    .cr_warn_nonconvergence(
      "covreg_comp_estimator", iter, max_iter, tol, c(N = N, P = P, Q = Q),
      c(dB = dB), check_every, (iter %/% check_every) * check_every,
      if (use_loglik && store_path) "parameter and likelihood changes" else "parameter changes",
      loglik = if (use_loglik) ll_now else NULL,
      loglik_change = if (use_loglik && store_path) abs(ll_now - ll_path[iter]) else NULL
    )
  }

  return(out)
}
