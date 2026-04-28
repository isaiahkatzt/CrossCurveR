#############################################
#             Shock Application             #
#############################################

st_build_shock_matrix <- function(curves, model_factors = c("L", "S", "C"),
                                  shock_curves, shock_factors, magnitude) {
  if (!is.character(curves) || length(curves) == 0 || anyNA(curves) || any(!nzchar(curves))) {
    stop("`curves` must be a non-empty character vector.")
  }

  if (anyDuplicated(curves)) {
    stop("`curves` must contain unique curve names.")
  }

  allowed_factor_sets <- list(
    c("L", "S"),
    c("L", "S", "C")
  )

  valid_model <- any(vapply(
    allowed_factor_sets,
    function(x) identical(model_factors, x),
    logical(1)
  ))

  if (!valid_model) {
    stop("`model_factors` must be exactly c('L', 'S') or c('L', 'S', 'C').")
  }

  if (!is.character(shock_curves) || length(shock_curves) == 0 ||
      anyNA(shock_curves) || any(!nzchar(shock_curves))) {
    stop("`shock_curves` must be a non-empty character vector.")
  }

  if (anyDuplicated(shock_curves)) {
    stop("`shock_curves` must contain unique curve names.")
  }

  if (!all(shock_curves %in% curves)) {
    stop("`shock_curves` must be a subset of `curves`.")
  }

  if (!is.character(shock_factors) || length(shock_factors) == 0 ||
      anyNA(shock_factors) || any(!nzchar(shock_factors))) {
    stop("`shock_factors` must be a non-empty character vector.")
  }

  if (anyDuplicated(shock_factors)) {
    stop("`shock_factors` must contain unique factor names.")
  }

  if (!all(shock_factors %in% model_factors)) {
    stop("`shock_factors` must be a subset of `model_factors`.")
  }

  shock_values <- rep(NA_real_, length(model_factors))
  names(shock_values) <- model_factors

  if (!is.numeric(magnitude) || anyNA(magnitude)) {
    stop("`magnitude` must be numeric and cannot contain NA values.")
  }

  if (length(magnitude) == 1) {
    if (length(shock_factors) != 1) {
      stop("Scalar `magnitude` is only allowed when exactly one factor is shocked.")
    }
    shock_values[shock_factors] <- as.numeric(magnitude)
  } else if (is.null(names(magnitude))) {
    if (length(magnitude) != length(shock_factors)) {
      stop("Unnamed `magnitude` must have the same length as `shock_factors`.")
    }
    shock_values[shock_factors] <- as.numeric(magnitude)
  } else {
    magnitude_names <- names(magnitude)

    if (anyNA(magnitude_names) || any(!nzchar(magnitude_names))) {
      stop("Named `magnitude` entries must all have non-empty names.")
    }

    if (anyDuplicated(magnitude_names)) {
      stop("Named `magnitude` entries must have unique names.")
    }

    if (length(magnitude) != length(shock_factors) ||
        !setequal(magnitude_names, shock_factors)) {
      stop("Named `magnitude` must match `shock_factors` exactly.")
    }

    shock_values[shock_factors] <- as.numeric(magnitude[shock_factors])
  }

  D <- length(curves)
  Beta_dim <- length(model_factors)

  shock_matrix <- matrix(NA_real_, nrow = Beta_dim, ncol = D)
  colnames(shock_matrix) <- curves
  rownames(shock_matrix) <- model_factors
  shock_matrix[, shock_curves] <- shock_values

  return(shock_matrix)
}

st_warn_large_shocks <- function(shock_matrix) {
  shock_thresholds <- c(L = 1, S = 0.5, C = 0.5)
  available_factors <- intersect(rownames(shock_matrix), names(shock_thresholds))

  has_large_shock <- any(vapply(available_factors, function(factor_name) {
    factor_values <- shock_matrix[factor_name, ]
    any(abs(factor_values[!is.na(factor_values)]) > shock_thresholds[[factor_name]])
  }, logical(1)))

  if (has_large_shock) {
    warning("Large shock magnitudes may yield unstable results.", call. = FALSE)
  }
}

st_resolve_shock_rows <- function(row_count, shock_window = NULL,
                                  shock_window_position = c("first", "last")) {
  shock_window_position <- match.arg(shock_window_position)

  if (is.null(shock_window)) {
    return(seq_len(row_count))
  }

  if (!is.numeric(shock_window) || length(shock_window) != 1 || is.na(shock_window) ||
      shock_window <= 0 || shock_window %% 1 != 0) {
    stop("`shock_window` must be NULL or a positive integer.")
  }

  shock_window <- as.integer(shock_window)

  if (shock_window > row_count) {
    stop("`shock_window` cannot exceed the number of available shock rows.")
  }

  if (shock_window_position == "first") {
    return(seq_len(shock_window))
  }

  seq.int(from = row_count - shock_window + 1L, to = row_count)
}

st_additive <- function(curves, curve_dns_factor, shock_curves, shock_factors,
                        shock_magnitude, shock_window = NULL,
                        shock_window_position = c("first", "last")) {
  if (!is.list(curve_dns_factor) || length(curve_dns_factor) == 0) {
    stop("`curve_dns_factor` must be a non-empty list.")
  }

  if (is.null(names(curve_dns_factor)) || anyNA(names(curve_dns_factor)) ||
      any(!nzchar(names(curve_dns_factor)))) {
    stop("`curve_dns_factor` must be a named list.")
  }

  factor_dim <- ncol(as.matrix(curve_dns_factor[[1]]))

  if (factor_dim == 2) {
    model_factors <- c("L", "S")
  } else if (factor_dim == 3) {
    model_factors <- c("L", "S", "C")
  } else {
    stop("`curve_dns_factor` entries must be 2- or 3-dimensional.")
  }

  shock_matrix <- st_build_shock_matrix(
    curves = curves,
    model_factors = model_factors,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    magnitude = shock_magnitude
  )

  st_warn_large_shocks(shock_matrix)
  shock_matrix[is.na(shock_matrix)] <- 0

  shocked_curve_dns_factor <- setNames(lapply(names(curve_dns_factor), function(curve_name) {
    curve_factor_mat <- as.matrix(curve_dns_factor[[curve_name]])

    if (ncol(curve_factor_mat) != factor_dim) {
      stop("All `curve_dns_factor` entries must have the same number of columns.")
    }

    if (!curve_name %in% colnames(shock_matrix)) {
      stop("Each named `curve_dns_factor` entry must correspond to a column in the shock matrix.")
    }

    shock_rows <- st_resolve_shock_rows(
      row_count = nrow(curve_factor_mat),
      shock_window = shock_window,
      shock_window_position = shock_window_position
    )

    shocked_mat <- curve_factor_mat
    shocked_mat[shock_rows, ] <- sweep(
      curve_factor_mat[shock_rows, , drop = FALSE],
      2,
      shock_matrix[, curve_name],
      FUN = "+"
    )

    shocked_mat
  }), names(curve_dns_factor))

  return(list(
    shock_matrix = shock_matrix,
    shocked_factors = shocked_curve_dns_factor
  ))
}

st_multiplicative <- function(curves, curve_dns_factor, shock_curves, shock_factors,
                              shock_magnitude, shock_window = NULL,
                              shock_window_position = c("first", "last")) {
  if (!is.list(curve_dns_factor) || length(curve_dns_factor) == 0) {
    stop("`curve_dns_factor` must be a non-empty list.")
  }

  if (is.null(names(curve_dns_factor)) || anyNA(names(curve_dns_factor)) ||
      any(!nzchar(names(curve_dns_factor)))) {
    stop("`curve_dns_factor` must be a named list.")
  }

  factor_dim <- ncol(as.matrix(curve_dns_factor[[1]]))

  if (factor_dim == 2) {
    model_factors <- c("L", "S")
  } else if (factor_dim == 3) {
    model_factors <- c("L", "S", "C")
  } else {
    stop("`curve_dns_factor` entries must be 2- or 3-dimensional.")
  }

  shock_matrix <- st_build_shock_matrix(
    curves = curves,
    model_factors = model_factors,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    magnitude = shock_magnitude
  )

  st_warn_large_shocks(shock_matrix)
  shock_matrix[is.na(shock_matrix)] <- 0
  shock_matrix <- shock_matrix + 1

  shocked_curve_dns_factor <- setNames(lapply(names(curve_dns_factor), function(curve_name) {
    curve_factor_mat <- as.matrix(curve_dns_factor[[curve_name]])

    if (ncol(curve_factor_mat) != factor_dim) {
      stop("All `curve_dns_factor` entries must have the same number of columns.")
    }

    if (!curve_name %in% colnames(shock_matrix)) {
      stop("Each named `curve_dns_factor` entry must correspond to a column in the shock matrix.")
    }

    shock_rows <- st_resolve_shock_rows(
      row_count = nrow(curve_factor_mat),
      shock_window = shock_window,
      shock_window_position = shock_window_position
    )

    shocked_mat <- curve_factor_mat
    shocked_mat[shock_rows, ] <- sweep(
      curve_factor_mat[shock_rows, , drop = FALSE],
      2,
      shock_matrix[, curve_name],
      FUN = "*"
    )

    shocked_mat
  }), names(curve_dns_factor))

  return(list(
    shock_matrix = shock_matrix,
    shocked_factors = shocked_curve_dns_factor
  ))
}
