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

vecY <- function(yields) {
  ## embed yield matrix Y as MD-dimensional vector 
  Y_tilde <- do.call(cbind, c(yields[1], lapply(yields[-1], function(y) y[, -1, drop = FALSE])))
  return(list(tY = Y_tilde[,-1], time = Y_tilde[,1]))
}

#### STOPPING POINT HERE; NOTE SHIFT IN VECY STRUCTURE 

PY_full <- function(vecY, P) {
  ## full timeframe permuted y 
  time <- vecY[, 1]
  vecY_matrix <- as.matrix(vecY[, -1]) 
  permuted_colnames <- colnames(vecY)[-1][apply(P, 1, which.max)]
  PYt <- vecY_matrix %*% t(P)
  PY <- data.frame(time = time, as.matrix(PYt))
  colnames(PY) <- c("time", permuted_colnames)
  
  return(PY)
}
