#############################################
#       Estimation Input Contracts          #
#############################################

.cc_names <- function(x, label) {
  if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x)) ||
      anyDuplicated(x)) {
    stop(label, " must contain unique, nonempty names.", call. = FALSE)
  }
  invisible(x)
}

.cc_mats <- function(mats) {
  if (!is.numeric(mats) || is.complex(mats) || !length(mats) ||
      any(!is.finite(mats)) || any(mats <= 0) || any(diff(mats) <= 0)) {
    stop("`mats` must be finite, positive, unique maturities in increasing order (months).",
         call. = FALSE)
  }
  invisible(mats)
}

.cc_decay <- function(x, label = "`lambdas`", scalar = FALSE) {
  if (!is.numeric(x) || is.complex(x) || !length(x) ||
      any(!is.finite(x)) || any(x <= 0) || (scalar && length(x) != 1L)) {
    stop(label, " must contain finite positive decay parameters",
         if (scalar) " (one value)." else ".", call. = FALSE)
  }
  invisible(x)
}

# Compare dates by instants, without changing the user's stored index. Numeric
# indices are observation indices, not implicitly Unix times or Excel dates.
.cc_index <- function(x, label = "Observation index") {
  if (inherits(x, "Date")) {
    values <- as.numeric(x) * 86400
    kind <- "instant:"
  } else if (inherits(x, c("POSIXct", "POSIXlt"))) {
    values <- as.numeric(as.POSIXct(x))
    kind <- "instant:"
  } else if (is.character(x)) {
    if (anyNA(x) || any(!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x))) {
      stop(label, " must use Date, POSIXct, ISO dates, or numeric observation indices.",
           call. = FALSE)
    }
    dates <- as.Date(x, format = "%Y-%m-%d")
    if (anyNA(dates) || any(as.character(dates) != x)) {
      stop(label, " contains invalid ISO dates.", call. = FALSE)
    }
    values <- as.numeric(dates) * 86400
    kind <- "instant:"
  } else if (is.numeric(x) && !is.complex(x)) {
    values <- as.numeric(x)
    kind <- "index:"
  } else {
    stop(label, " has an unsupported type.", call. = FALSE)
  }
  if (!length(values) || any(!is.finite(values)) || any(diff(values) <= 0)) {
    stop(label, " must be nonmissing, unique, and strictly increasing.", call. = FALSE)
  }
  paste0(kind, sprintf("%.17g", values))
}

.cc_named <- function(x, keys, label) {
  .cc_names(keys, "Curve names")
  if (!is.list(x) || is.data.frame(x) || length(x) != length(keys)) {
    stop(label, " must be a list with one entry per curve.", call. = FALSE)
  }
  if (is.null(names(x))) {
    # Retain explicitly labelled positional lists, assigning identities once.
    names(x) <- keys
  }
  .cc_names(names(x), paste0(label, " names"))
  if (!setequal(names(x), keys)) {
    stop(label, " names must match `curves` exactly.", call. = FALSE)
  }
  x[keys]
}

.cc_reference <- function(reference, curves) {
  if (!is.character(reference) || length(reference) != 1L || is.na(reference) ||
      !reference %in% curves) {
    stop("`reference` must identify one of the supplied curves.", call. = FALSE)
  }
  invisible(reference)
}

.cc_tenor_order <- function(labels, mats, label) {
  .cc_names(labels, paste0(label, " maturity names"))
  months <- vapply(labels, .matname_to_months, numeric(1), USE.NAMES = FALSE)
  order <- vapply(mats, function(m) {
    hit <- which(abs(months - m) <= 1e-10 * max(1, m))
    if (length(hit) != 1L) {
      stop(label, " maturity columns must match `mats` exactly, without duplicate tenors.",
           call. = FALSE)
    }
    hit
  }, integer(1))
  if (length(labels) != length(mats) || anyDuplicated(order)) {
    stop(label, " maturity columns must match `mats` exactly.", call. = FALSE)
  }
  order
}

.cc_matrix <- function(x, label) {
  if (length(dim(x)) > 2L) {
    stop(label, " must be two dimensional.", call. = FALSE)
  }
  if (is.data.frame(x) && !all(vapply(x, is.numeric, logical(1)))) {
    stop(label, " must contain finite numeric values.", call. = FALSE)
  }
  x <- as.matrix(x)
  if (!is.numeric(x) || is.complex(x) || length(dim(x)) != 2L ||
      !nrow(x) || !ncol(x) || any(!is.finite(x))) {
    stop(label, " must be a nonempty finite numeric matrix.", call. = FALSE)
  }
  if (!is.null(colnames(x))) .cc_names(colnames(x), paste0(label, " column names"))
  x
}

.cc_panel <- function(panel, mats, label = "Yield panel") {
  .cc_mats(mats)
  if (!is.data.frame(panel) || ncol(panel) != length(mats) + 1L) {
    stop(label, " must be a data frame with an index first and one column per maturity.",
         call. = FALSE)
  }
  .cc_names(names(panel), paste0(label, " column names"))
  .cc_index(panel[[1]], paste0(label, " observation index"))
  order <- .cc_tenor_order(names(panel)[-1], mats, label)
  if (!identical(order, seq_along(mats))) panel <- panel[c(1L, order + 1L)]
  .cc_matrix(panel[-1], label)
  panel
}

.cc_rates <- function(x, mats, label = "Yield matrix") {
  .cc_mats(mats)
  x <- .cc_matrix(x, label)
  if (ncol(x) != length(mats)) {
    stop(label, " must have one column per maturity.", call. = FALSE)
  }
  if (!is.null(colnames(x))) {
    order <- .cc_tenor_order(colnames(x), mats, label)
    if (!identical(order, seq_along(mats))) x <- x[, order, drop = FALSE]
  }
  x
}

.cc_yields <- function(yields, mats, curves = names(yields), reference = NULL) {
  .cc_mats(mats)
  yields <- .cc_named(yields, curves, "`yields`")
  if (!is.null(reference)) .cc_reference(reference, curves)
  for (curve in curves) {
    yields[[curve]] <- .cc_panel(yields[[curve]], mats, paste0("Curve `", curve, "`"))
  }
  index <- .cc_index(yields[[1]][[1]])
  for (curve in curves[-1]) {
    if (!identical(index, .cc_index(yields[[curve]][[1]]))) {
      stop("Curve `", curve, "` observation index does not match the other curves. ",
           "Align observations explicitly before fitting.", call. = FALSE)
    }
  }
  yields
}

.cc_positional_rows <- function(rows) {
  is.null(rows) || (!anyNA(rows) && all(grepl("^[0-9]+$", rows)))
}

.cc_features <- function(X, time = NULL, label = "`Xt`") {
  X <- .cc_matrix(X, label)
  if (!is.null(time)) {
    index <- .cc_index(time)
    if (nrow(X) != length(time)) {
      stop(label, " must have one row per observation.", call. = FALSE)
    }
    rows <- rownames(X)
    if (!.cc_positional_rows(rows) &&
        !identical(rows, as.character(time))) {
      if (!identical(.cc_index(rows, paste0(label, " row index")), index)) {
        stop(label, " row index does not match the yield observation index.", call. = FALSE)
      }
    }
  }
  X
}

.cc_input_spec <- function(yields, mats, curves, X = NULL, feature_spec = NULL, decay = NULL) {
  list(
    version = 1L,
    curves = curves,
    time = yields[[1]][[1]],
    mats = mats,
    units = list(maturity = "months", yield = "percent",
                 ns_decay = "inverse months", svensson_decay = "months",
                 covariance = "squared percentage points"),
    rate_definition = "as supplied; no compounding or par-to-zero conversion",
    observation_step = "one successive observation",
    decay = decay,
    features = if (!is.null(X)) list(
      columns = colnames(X), count = ncol(X),
      row_alignment = if (.cc_positional_rows(rownames(X))) "positional" else "indexed",
      preprocessing = feature_spec
    ) else NULL
  )
}

#############################################
#        Precomputed Input Contracts        #
#############################################

.cc_dly_sc <- function(sc, mats, curves) {
  sc <- .cc_named(sc, curves, "Country fits")
  if (is.null(mats)) {
    mats <- sort(vapply(sc[[1]]$yield_names, .matname_to_months,
                        numeric(1), USE.NAMES = FALSE))
  }
  .cc_mats(mats)
  index <- NULL
  for (curve in curves) {
    fit <- sc[[curve]]
    label <- paste0("Country fit `", curve, "`")
    current_index <- .cc_index(fit$time, paste0(label, " observation index"))
    if (!is.null(index) && !identical(index, current_index)) {
      stop(label, " observation index does not match the other country fits.", call. = FALSE)
    }
    index <- current_index
    .cc_decay(fit$lambda, paste0(label, " decay"), scalar = TRUE)
    order <- .cc_tenor_order(fit$yield_names, mats, label)
    expected <- c(length(fit$time), length(mats))
    beta <- .cc_features(fit$betas$NSbetas, fit$time, paste0(label, " factors"))
    if (!identical(dim(beta), c(expected[1], 2L))) {
      stop(label, " factors must have one row per observation and two columns.", call. = FALSE)
    }
    for (part in c("phi", "cross_phi")) {
      mat <- .cc_matrix(fit$phi[[part]], paste0(label, " ", part))
      if (!identical(dim(mat), c(2L, expected[2]))) {
        stop(label, " loadings must have two rows and one column per maturity.", call. = FALSE)
      }
      if (!identical(order, seq_along(mats))) fit$phi[[part]] <- mat[, order, drop = FALSE]
    }
    for (part in c("m_yield", "NSfit", "W")) {
      mat <- if (part == "m_yield") fit$betas[[part]] else fit$nsfit[[part]]
      mat <- .cc_features(mat, fit$time, paste0(label, " ", part))
      if (!identical(dim(mat), expected)) {
        stop(label, " yield and residual matrices must match its dates and maturities.", call. = FALSE)
      }
      if (!identical(order, seq_along(mats))) {
        if (part == "m_yield") fit$betas[[part]] <- mat[, order, drop = FALSE]
        else fit$nsfit[[part]] <- mat[, order, drop = FALSE]
      }
    }
    fit$yield_names <- fit$yield_names[order]
    sc[[curve]] <- fit
  }
  list(sc = sc, mats = mats)
}

.cc_cov_inputs <- function(W, X, S0, B) {
  W <- .cc_matrix(W, "`W`")
  if (ncol(W) < 2L) {
    stop("Covariance regression requires at least two response columns in `W`; univariate responses are not supported.", call. = FALSE)
  }
  X <- .cc_matrix(X, "`X`")
  if (nrow(W) != nrow(X)) stop("W and X must have the same number of rows.", call. = FALSE)
  wr <- rownames(W)
  xr <- rownames(X)
  if (!.cc_positional_rows(wr) && !.cc_positional_rows(xr) && !identical(wr, xr)) {
    stop("W and X row indices must match in the same order.", call. = FALSE)
  }
  align <- function(mat, target, margin, label) {
    supplied <- dimnames(mat)[[margin]]
    if (!is.null(supplied)) .cc_names(supplied, paste0(label, " names"))
    if (!is.null(target) && !is.null(supplied)) {
      if (!setequal(target, supplied)) stop(label, " names do not match the design.", call. = FALSE)
      if (margin == 1L) mat <- mat[match(target, supplied), , drop = FALSE]
      else mat <- mat[, match(target, supplied), drop = FALSE]
    }
    mat
  }
  if (!is.null(S0)) {
    S0 <- .cc_matrix(S0, "`S0`")
    if (!identical(dim(S0), c(ncol(W), ncol(W)))) {
      stop("S0 must have one row and column per response in W.", call. = FALSE)
    }
    S0 <- align(align(S0, colnames(W), 1L, "S0 rows"), colnames(W), 2L, "S0 columns")
    if (max(abs(S0 - t(S0))) > 1e-8 * max(abs(S0))) {
      stop("S0 must be symmetric.", call. = FALSE)
    }
    eigenvalues <- eigen(0.5 * (S0 + t(S0)), symmetric = TRUE, only.values = TRUE)$values
    if (min(eigenvalues) < -1e-8 * max(abs(eigenvalues))) {
      stop("S0 must be positive semidefinite before regularization.", call. = FALSE)
    }
  }
  if (!is.null(B)) {
    B <- .cc_matrix(B, "`B`")
    if (!identical(dim(B), c(ncol(W), ncol(X)))) {
      stop("B must have dimensions ncol(W) x ncol(X).", call. = FALSE)
    }
    B <- align(align(B, colnames(W), 1L, "B rows"), colnames(X), 2L, "B columns")
  }
  list(W = W, X = X, S0 = S0, B = B)
}

.cc_cov_spec <- function(W, X) {
  list(response_columns = colnames(W), feature_columns = colnames(X),
       row_index = rownames(W), feature_row_index = rownames(X),
       response_units = "as supplied", feature_units = "as supplied",
       normalization = "none", intercept = "as supplied")
}
