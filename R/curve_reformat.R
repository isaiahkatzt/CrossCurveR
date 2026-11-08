#############################################
#      Initial Yield Curve Formatting       #
#############################################

full_matrix_yield <- function(yields, curves) {
  yields <- .cc_named(yields, curves, "`yields`")
  mats <- sort(vapply(names(yields[[1]])[-1], .matname_to_months, numeric(1), USE.NAMES = FALSE))
  yields <- .cc_yields(yields, mats, curves)
  stats::setNames(lapply(yields, function(yield) {
    as.matrix(yield[-1])}
    ), curves)
}

vecY <- function(yields, curves, mats, ts=TRUE) {
  .cc_mats(mats)
  ## embed yield matrix Y as MD-dimensional vector
  mat_str <- sapply(mats, numeric_to_matname)
  mds <- as.vector(Matrix::t(outer(curves, mat_str, paste, sep='.')))

  if (!ts) {
    yields <- .cc_named(yields, curves, "`yields`")
    yields <- lapply(yields, .cc_rates, mats = mats)
    if (length(unique(vapply(yields, nrow, integer(1)))) != 1L) {
      stop("Yield matrices must have the same number of observations.", call. = FALSE)
    }
    Y_tilde <- do.call(cbind, yields)
    colnames(Y_tilde) <- mds
    return(list(tY=Y_tilde))
  }
  else {
    yields <- .cc_yields(yields, mats, curves)
    Y_tilde <- do.call(cbind, c(yields[1], lapply(yields[-1], function(y) y[, -1, drop = FALSE])))
    colnames(Y_tilde) <- c("time", mds)
    return(list(tY=Y_tilde[,-1], time=Y_tilde[,1]))
  }
}

#############################################
#        Permutation and Curve Cuts         #
#############################################

PY_full <- function(vecY, P) {
  ## yields x permutation matrix
  vm_Y <- as.matrix(vecY)
  permuted_colnames <- colnames(vecY)[apply(P, 1, which.max)]
  PYt <- vm_Y %*% Matrix::t(P)
  PY <- data.frame(as.matrix(PYt))
  colnames(PY) <- permuted_colnames
  return(list(py = PY, cn = permuted_colnames))
}

cf_curve_cut <- function(vecY, curves, mats) {
  M <- length(mats)
  D <- length(curves)

  cut_points <- rep(seq_len(D), each = M)

  cut_curve <- split(seq_len(ncol(vecY)), cut_points)
  cut_curve <- lapply(cut_curve, function(j) vecY[, j, drop = FALSE])

  names(cut_curve) <- curves
  return(cut_curve)
}

cf_phi_cut <- function(vphi, curves, mats) {
  M <- length(mats)
  D <- length(curves)

  cut_points <- rep(seq_len(D), each = M)

  cut_phi <- lapply(seq_len(D), function(j) {
    lapply(vphi, function(cuts) cuts[cut_points == j, , drop = FALSE] )
  })

  names(cut_phi) <- curves
  return(cut_phi)
}
