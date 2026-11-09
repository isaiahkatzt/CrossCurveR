#############################################
#       Cached Feature Compatibility        #
#############################################

.st_check_cached_coefficients <- function(BS0, X, curves, mats) {
  if (!is.list(BS0) || length(BS0) != length(mats)) {
    stop("Cached BS0 must contain one coefficient block per baseline maturity.", call. = FALSE)
  }
  if (!is.null(names(BS0)) && !identical(
      .cc_tenor_order(names(BS0), mats, "Cached BS0"), seq_along(mats))) {
    stop("Cached BS0 maturity order does not match the baseline design.", call. = FALSE)
  }
  for (fit in BS0) {
    B <- .cc_matrix(fit$B, "Cached B")
    S0 <- .cc_matrix(fit$S0, "Cached S0")
    if (!identical(dim(B), c(length(curves), ncol(X))) ||
        !identical(dim(S0), c(length(curves), length(curves)))) {
      stop("Cached covariance coefficient dimensions do not match the baseline design.", call. = FALSE)
    }
    axes <- list(list(rownames(B), curves), list(colnames(B), colnames(X)),
                 list(rownames(S0), curves), list(colnames(S0), curves))
    if (any(vapply(axes, function(axis) !is.null(axis[[1]]) && !is.null(axis[[2]]) &&
        !identical(axis[[1]], axis[[2]]), logical(1)))) {
      stop("Cached covariance coefficient identities/order do not match the baseline design.",
           call. = FALSE)
    }
    recorded <- fit$input_spec$feature_columns
    if (!is.null(recorded) && !identical(recorded, colnames(X))) {
      stop("Cached covariance coefficients were fitted on a different feature layout.", call. = FALSE)
    }
  }
  invisible(NULL)
}
