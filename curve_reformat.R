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
