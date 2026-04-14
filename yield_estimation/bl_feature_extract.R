 #############################################
#            Feature Extraction             #
#############################################

blmc_phi <- function(bl_sc_list){
  return(t(do.call(cbind, lapply(bl_sc_list$phi, `[[`, "phi"))))
}

blmc_VAR <- function(blsc_list, reference, curves) {
  marginals <- setdiff(curves, reference) 
  reference_beta <- blsc_list$betas[[reference]]
  cross_curve_beta <- setNames(vector("list", length(marginals)), marginals)
  
  for (m in marginals) {
    m_beta <- blsc_list$betas[[m]]
    cc_beta <- cbind(reference_beta$NSbetas, m_beta$NSbetas)
    colnames(cc_beta) <- c(
      paste0(reference, c(".L", ".S", ".C")),
      paste0(m, c(".L", ".S", ".C"))
    )
    cross_curve_beta[[m]] <- cc_beta
  }
  return(cross_curve_beta)
}

blmc_identity_sigma <- function(curves, mats) {
  ## test function for methodological verification  
  M = length(mats) 
  D = length(curves) 
  identity_sigma <- diag(M * D)
  return(identity_sigma) 
}

blmc_ECM <- function(cc_betas, estim = "ML", type = "eigen", alpha = 0.1) {
  lag_opt <- max(VARselect(cc_betas, type = "const")$selection["SC(n)"], 2)
  vecm_unrestricted <- VECM(cc_betas, lag = lag_opt - 1, estim = estim, include = "const")
  r_test <- rank.test(vecm_unrestricted, type = type, cval = alpha)
  r <- r_test$r
  K <- ncol(cc_betas)
  
  ## tsDyn::VECM requires 1 <= r <= K - 1. Treat degenerate selections
  ## as "no usable cointegration restriction" for downstream Xt construction.
  if (!is.finite(r) || r <= 0 || r >= K) {
    return(list(Gmatrix = 0, Fmatrix = 0, rank = 0))
  }
  
  vecm_restricted <- VECM(cc_betas, lag = lag_opt - 1, r = r, estim = estim, include = "const")
  
  ## round for floating point precision 
  Fmat <- round(coefA(vecm_restricted), 14)
  Gmat <- round(coefB(vecm_restricted), 14)
  return(list(Fmatrix = Fmat, Gmatrix = Gmat, rank = r))
}

blmc_cspread <- function(cc_betas, Gmatrix, curve_name = NULL, normalize = TRUE) {
  Xt <- cc_betas %*% Gmatrix 

  if (is.null(curve_name)) {
    curve_labels <- unique(sub("\\..*$", "", colnames(cc_betas)))
    if (length(curve_labels) >= 2) {
      curve_name <- curve_labels[[2]]
    }
  }

  spread_names <- colnames(Xt)
  if (is.null(spread_names)) {
    spread_names <- paste0("r", seq_len(ncol(Xt)))
  }

  if (!is.null(curve_name) && nzchar(curve_name)) {
    colnames(Xt) <- paste0(curve_name, ".", spread_names)
  } else {
    colnames(Xt) <- spread_names
  }

  if (normalize) {
    Xt <- sweep(Xt, 2, colMeans(Xt), FUN = "-")
    Xt <- sweep(Xt, 2, apply(Xt, 2, sd), FUN = "/")
  }
  return(Xt)
}

blmc_Fselect <- function(Fmatrix) {
  reversion_norm <- sqrt(colSums(Fmatrix^2))
  return(which.max(reversion_norm))
}

blcc_cointegration <- function(blsc_list, curves, reference, estim = "ML", type = "eigen", alpha = 0.1, normalize = TRUE) {
  ## cross curve beta list 
  cc_betas <- blmc_VAR(blsc_list, reference, curves) 
  
  ## curve-specific error correct model construction 
  full_ecm <- lapply(cc_betas, function(cc_var) blmc_ECM(cc_betas = cc_var, estim = estim, type = type, alpha = alpha))
  
  ## cointegration spread construction 
  full_cspread <- setNames(lapply(seq_along(cc_betas), function(curve) {
    if (full_ecm[[curve]]$rank != 0) {
      blmc_cspread(
        cc_betas = cc_betas[[curve]],
        Gmatrix = full_ecm[[curve]]$Gmatrix,
        curve_name = names(cc_betas)[[curve]],
        normalize = normalize
      )
    }
  }), names(cc_betas))
  
  return(list(
    cc_cspread = full_cspread, 
    cc_ecm = full_ecm, 
    cc_beta = cc_betas
  ))
}

blcc_build_Xt <- function(cc_cspread, cc_ecm, trunc = FALSE) {
  ## construct full Xt covariate matrix 
  valid_idx <- which(vapply(cc_ecm, function(x) {
    is.list(x) && is.numeric(x$rank) && length(x$rank) == 1 && x$rank > 0
  }, logical(1)))
  
  if (length(valid_idx) == 0) {
    stop("No usable cointegration relationships were identified for Xt construction.")
  }
  
  if (!trunc) {
    Xt <- do.call(cbind, cc_cspread[valid_idx]) 
  } else {
    Fmatrix_list <- lapply(cc_ecm[valid_idx], `[[`, "Fmatrix")
    Fnorms <- lapply(Fmatrix_list, blmc_Fselect)
    Xt <- do.call(cbind, lapply(seq_along(valid_idx), function(i) {
      cc_cspread[[valid_idx[[i]]]][, Fnorms[[i]], drop = FALSE]
    }))
  }
  return(Xt) 
}

blcc_build_Wj <- function(curve_nsfit, mats){
  ## construct tenor-specific Wj list 
  wc_list <- lapply(curve_nsfit, `[[`, 2) 
  mat_names <- paste0("X", sapply(mats, numeric_to_matname))
  wj_full <- setNames(lapply(seq_along(mats), function(i) {
    do.call(cbind, lapply(wc_list, function(x) x[, i])) 
  }), mat_names)
  
  return(Wj = wj_full)
}

blcc_fe <- function(Xt, Wt, mats, covreg=zero_mean_covreg_em,
                    init = "adaptive", max_iter = 1000, tol = 1e-10, S0 = NULL, B = NULL,
                    S0_shrink_diag = 0, verb = FALSE, term = FALSE) {
  ## baseline cross-curve feature extraction
  BS0_list <- lapply(Wt, function(Wj) {
    covreg(W=Wj, X=Xt, init = init, max_iter = max_iter, tol = tol, S0 = S0, B = B, 
           ridge_S0=1e-8, ridge_X = 1e-10, S0_shrink_diag = S0_shrink_diag,
           check_every=1, use_loglik=TRUE, verbose=FALSE, store_path=TRUE)
           #verb = verb, term = term)
  })
  return(BS0_list) 
}

sigma_jt_optim <- function(BS0j, x) {
  ## compute single-day sigma_{jt}   
  B <- BS0j$B
  S0 <- BS0j$S0
  Bx <- B %*% x
  Sigma <- S0 + tcrossprod(Bx)  
  return(Sigma)
}

jSigma_optim <- function(BS0j, X) {
  # compute full period sigma_{jt}
  N <- nrow(X)
  day_list <- vector("list", N)
  
  jSigma <- lapply(seq_len(N), function(i) {
    sigma_jt_optim(BS0j = BS0j, x = X[i, ])
  })
  return(jSigma)
}

tSigma_optim <- function(BS0, X) {
  # compute full period, all sigma_{jt} 
  N <- nrow(X) 
  sigma_j_full <- lapply(BS0, function(j) jSigma_optim(j, X))
  sigma_t_full <- lapply(seq_len(N), function(t) {
    bdiag(lapply(sigma_j_full, function(j) j[[t]]))
  })
  return(list(tSigma = sigma_t_full, jSigma = sigma_j_full))
}

#############################################
#      Cointegration Spread Smoothers       #
#############################################

fe_kp_quantile <- function(time, knot_count = 5) {
  x <- as.numeric(as.Date(time))
  probs <- seq(0, 1, length.out = knot_count + 2)  
  knots <- quantile(x, probs = probs)[-c(1, knot_count + 2)] 
  
  return(knots)
}

fe_spline_basis <- function(spline_type = c("ns", "bs"), degree = 3) {
  spline_type <- match.arg(spline_type)
  spline_basis <- switch(
    spline_type,
    ns = function(x, knots) splines::ns(x, knots = knots),
    bs = function(x, knots) splines::bs(x, knots = knots, degree = degree),
  )
  return(spline_basis) 
}

fe_xt_spline <- function(time, xt, knot_points, spline_fit = c("ns", "bs"), degree = 3) {
  ## cspread smoothing  
  spline_basis <- fe_spline_basis(spline_fit, degree)
  time <- as.numeric(as.Date(time))
  
  basis_matrix <- spline_basis(time, knots = knot_points)
  
  spline_model <- lm(xt ~ basis_matrix)
  xt_smooth <- predict(spline_model)
  return(xt_smooth) 
}

fe_xt_rm <- function(X, k = 21, passes = 2) {
  h <- (k - 1) %/% 2
  w <- rep(1, k)
  
  apply(X, 2, function(col) {
    y <- col
    n <- length(y)
    for (p in seq_len(passes)) {
      num <- convolve(y, w, type = "open")[(h + 1):(h + n)]
      den <- convolve(rep(1, n), w, type = "open")[(h + 1):(h + n)]
      y <- num / den
    }
    y
  })
}

fe_xt_hp <- function(X, lambda = 1600) {
  if (!is.numeric(lambda) || length(lambda) != 1 || is.na(lambda) || lambda < 0) {
    stop("lambda must be a non-negative numeric scalar")
  }
  
  vector_input <- is.null(dim(X))
  X_mat <- if (vector_input) matrix(X, ncol = 1) else as.matrix(X)
  storage.mode(X_mat) <- "double"
  
  n <- nrow(X_mat)
  if (n < 3 || lambda == 0) {
    xt_smooth <- X_mat
  } else {
    D <- diff(diag(n), differences = 2)
    hp_matrix <- diag(n) + lambda * crossprod(D)
    xt_smooth <- solve(hp_matrix, X_mat)
  }
  
  dimnames(xt_smooth) <- dimnames(X_mat)
  
  if (vector_input) {
    return(as.vector(xt_smooth))
  }
  return(xt_smooth)
}

fe_xt_henderson <- function(X, k = 13) {
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k <= 0 || k %% 2 == 0 || k != as.integer(k)) {
    stop("k must be a positive odd integer")
  }
  
  vector_input <- is.null(dim(X))
  X_mat <- if (vector_input) matrix(X, ncol = 1) else as.matrix(X)
  storage.mode(X_mat) <- "double"
  
  h <- (k - 1) %/% 2
  offsets <- -h:h
  m <- h
  weights <- 315 *
    ((m + 1)^2 - offsets^2) *
    ((m + 2)^2 - offsets^2) *
    ((m + 3)^2 - offsets^2) *
    (3 * (m + 2)^2 - 11 * offsets^2 - 16) /
    (
      8 * (m + 2) *
        ((m + 2)^2 - 1) *
        (4 * (m + 2)^2 - 1) *
        (4 * (m + 2)^2 - 9) *
        (4 * (m + 2)^2 - 25)
    )
  
  n <- nrow(X_mat)
  xt_smooth <- matrix(NA_real_, nrow = n, ncol = ncol(X_mat), dimnames = dimnames(X_mat))
  
  for (j in seq_len(ncol(X_mat))) {
    for (i in seq_len(n)) {
      idx <- max(1, i - h):min(n, i + h)
      weight_idx <- idx - i + h + 1
      local_weights <- weights[weight_idx]
      xt_smooth[i, j] <- sum(X_mat[idx, j] * local_weights) / sum(local_weights)
    }
  }
  
  if (vector_input) {
    return(as.vector(xt_smooth))
  }
  return(xt_smooth)
}
