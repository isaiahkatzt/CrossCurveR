#### initial yield curve formatting   

full_matrix_yield <- function(yields, curves) {
  setNames(lapply(yields, function(yield) {
    as.matrix(yield[-1])}
    ), curves)
}

numeric_to_matname <- function(numeric_mat) {
  ## convert numeric maturities to character string 
  if (numeric_mat < 12) {
    prefix <- ifelse(numeric_mat < 10, "0", "")
    colname <- paste0(prefix, as.character(numeric_mat), "M")
  } else {
    numeric_mat_month <- numeric_mat / 12
    prefix <- ifelse(numeric_mat_month < 10, "0", "")
    colname <- paste0(prefix, as.character(numeric_mat_month), "Y")
  }
  return(colname)
}

vecY <- function(yields, curves, mats) {
  ## embed yield matrix Y as MD-dimensional vector 
  mat_str <- sapply(mats, numeric_to_matname)
  mds <- as.vector(t(outer(curves, mat_str, paste, sep='.')))
  Y_tilde <- do.call(cbind, c(yields[1], lapply(yields[-1], function(y) y[, -1, drop = FALSE])))
  colnames(Y_tilde) <- c("time", mds) 
    
  return(list(tY = Y_tilde[,-1], time = Y_tilde[,1]))
}

PY_full <- function(vecY, P) {
  vm_Y <- as.matrix(vecY) 
  permuted_colnames <- colnames(vecY)[apply(P, 1, which.max)]
  PYt <- vm_Y %*% t(P)
  PY <- data.frame(as.matrix(PYt))
  colnames(PY) <- permuted_colnames
  return(list(py = PY, cn = permuted_colnames))
}

