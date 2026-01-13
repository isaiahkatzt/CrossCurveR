
#############################################
#         Full Cross-Curve Estimation       #
#############################################

blce_rescale_yield <- function(normalized, sqrt_sigma_t, P, ptype = c("tenor", "curve"), time, curves, mats) {
  ptype = match.arg(ptype)
  vec_normalized <- vecY(normalized, curves, mats, ts=FALSE)[['tY']]
  N <- nrow(vec_normalized)
  rescale_yield_matrix <- matrix(NA, nrow = N, ncol = ncol(vec_normalized))
  ## rescaling and curve / tenor permutation 
  
  if (ptype == "tenor") {
    permuted_colnames <- c("time", colnames(vec_normalized)[apply(P, 1, which.max)])
    
    for (dt in seq_len(N)) {
      rescale_yield_matrix[dt, ] <- 
        as.numeric(t(sqrt_sigma_t[[dt]] %*% P %*% as.numeric(vec_normalized[dt, ])))
    }
  }
  
  else if (ptype == "curve") {
    permuted_colnames <- c("time", colnames(vec_normalized))
    
    for (dt in seq_len(N)) {
      rescale_yield_matrix[dt, ] <- 
        as.numeric(t(t(P) %*% sqrt_sigma_t[[dt]] %*% P %*% as.numeric(vec_normalized[dt, ])))
    }
  }
  
  rescale_yield <- data.frame(time, rescale_yield_matrix)
  colnames(rescale_yield) <- permuted_colnames
  return(rescale_yield) 
}

