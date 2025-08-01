#### initial yield curve formatting   

full_matrix_yield <- function(yields, curves) {
  setNames(lapply(yields, function(yield) {
    as.matrix(yield[-1])}
    ), curves)
}