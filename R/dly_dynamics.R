#############################################
#      State and Transition Estimation      #
#############################################

dly_state_matrix <- function(level_series, slope_series, row_template = NULL) {
  state_mat <- cbind(L = as.numeric(level_series), S = as.numeric(slope_series))
  if (!is.null(row_template)) {
    rownames(state_mat) <- row_template
  }
  state_mat
}

dly_fit_var1 <- function(state_matrix) {
  state_matrix <- as.matrix(state_matrix)

  if (nrow(state_matrix) < 2) {
    stop("DLY global dynamics: at least two observations are required to form lagged pairs; received ",
         nrow(state_matrix), ".", call. = FALSE)
  }

  row_template <- rownames(state_matrix)
  lag_state <- state_matrix[-nrow(state_matrix), , drop = FALSE]
  response <- state_matrix[-1, , drop = FALSE]

  X <- cbind(Intercept = 1, lag_state)
  coef_mat <- apply(response, 2, function(y) {
    stats::lm.fit(X, y)$coefficients
  })

  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, ncol = 1)
    colnames(coef_mat) <- colnames(response)
  }

  if (any(!is.finite(coef_mat))) {
    stop("Non-finite coefficients generated in DLY global dynamics; check for rank-deficient regressors.",
         call. = FALSE)
  }

  intercept <- as.numeric(coef_mat["Intercept", ])
  names(intercept) <- colnames(response)
  phi <- Matrix::t(coef_mat[colnames(lag_state), , drop = FALSE])
  colnames(phi) <- colnames(lag_state)
  rownames(phi) <- colnames(response)

  fitted_core <- X %*% coef_mat
  fitted_state <- state_matrix
  fitted_state[-1, ] <- fitted_core
  residual_state <- state_matrix - fitted_state
  rownames(fitted_state) <- row_template
  rownames(residual_state) <- row_template

  list(
    intercept = intercept,
    phi = phi,
    fitted = fitted_state,
    residual = residual_state
  )
}

dly_fit_rs <- function(residual_matrix, global_matrix) {
  residual_matrix <- as.matrix(residual_matrix)
  global_matrix <- as.matrix(global_matrix)

  if (nrow(residual_matrix) != nrow(global_matrix)) {
    stop("residual_matrix and global_matrix must align by time.")
  }

  if (nrow(residual_matrix) < 2) {
    stop("DLY residual dynamics: at least two observations are required to form lagged pairs; received ",
         nrow(residual_matrix), ".", call. = FALSE)
  }

  row_template <- rownames(residual_matrix)
  residual_lag <- residual_matrix[-nrow(residual_matrix), , drop = FALSE]
  global_lag <- global_matrix[-nrow(global_matrix), , drop = FALSE]
  response <- residual_matrix[-1, , drop = FALSE]

  colnames(residual_lag) <- paste0(colnames(residual_lag), "_lag")
  colnames(global_lag) <- paste0("g", colnames(global_lag), "_lag")
  X_no_intercept <- cbind(residual_lag, global_lag)
  X <- cbind(Intercept = 1, X_no_intercept)

  coef_mat <- apply(response, 2, function(y) {
    stats::lm.fit(X, y)$coefficients
  })

  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, ncol = 1)
    colnames(coef_mat) <- colnames(response)
  }

  if (any(!is.finite(coef_mat))) {
    stop("Non-finite coefficients generated in DLY residual dynamics; check for rank-deficient regressors.",
         call. = FALSE)
  }

  intercept <- as.numeric(coef_mat["Intercept", ])
  names(intercept) <- colnames(response)
  phi <- Matrix::t(coef_mat[paste0(colnames(residual_matrix), "_lag"), , drop = FALSE])
  kappa <- Matrix::t(coef_mat[paste0("g", colnames(global_matrix), "_lag"), , drop = FALSE])
  colnames(phi) <- colnames(residual_matrix)
  rownames(phi) <- colnames(response)
  colnames(kappa) <- colnames(global_matrix)
  rownames(kappa) <- colnames(response)

  fitted_core <- X %*% coef_mat
  fitted_state <- residual_matrix
  fitted_state[-1, ] <- fitted_core
  residual_state <- residual_matrix - fitted_state
  rownames(fitted_state) <- row_template
  rownames(residual_state) <- row_template

  list(
    intercept = intercept,
    phi = phi,
    kappa = kappa,
    fitted = fitted_state,
    residual = residual_state
  )
}

#############################################
#     State Transitions and Propagation     #
#############################################

dly_step_var1 <- function(prev_state, intercept, phi) {
  as.numeric(intercept + phi %*% as.numeric(prev_state))
}

dly_step_rs <- function(prev_residual, prev_global, intercept, phi, kappa) {
  as.numeric(intercept + phi %*% as.numeric(prev_residual) + kappa %*% as.numeric(prev_global))
}

dly_propagate_rs <- function(residual_models, actual_residuals, global_path,
                             curves) {
  stats::setNames(lapply(curves, function(curve) {
    model <- residual_models[[curve]]
    actual_state <- as.matrix(actual_residuals[[curve]])
    propagated_state <- actual_state

    for (t in seq_len(nrow(actual_state))[-1]) {
      propagated_state[t, ] <- dly_step_rs(
        prev_residual = propagated_state[t - 1, ],
        prev_global = global_path[t - 1, ],
        intercept = model$intercept,
        phi = model$phi,
        kappa = model$kappa
      )
    }

    rownames(propagated_state) <- rownames(actual_state)
    colnames(propagated_state) <- colnames(actual_state)
    propagated_state
  }), curves)
}

#############################################
#              State Unpacking              #
#############################################

dly_unpack_country_state <- function(country_state_list, curves) {
  level_factor <- do.call(cbind, lapply(curves, function(curve) country_state_list[[curve]][, "L"]))
  slope_factor <- do.call(cbind, lapply(curves, function(curve) country_state_list[[curve]][, "S"]))
  colnames(level_factor) <- curves
  colnames(slope_factor) <- curves

  row_template <- rownames(country_state_list[[curves[[1]]]])
  rownames(level_factor) <- row_template
  rownames(slope_factor) <- row_template

  list(level = level_factor, slope = slope_factor)
}

#############################################
#           Compatibility Aliases           #
#############################################
dly_fit_var1_block <- dly_fit_var1
dly_fit_country_dynamic_block <- dly_fit_rs
dly_transition_var1 <- dly_step_var1
dly_transition_country_residual <- dly_step_rs
dly_propagate_residual_states <- dly_propagate_rs
