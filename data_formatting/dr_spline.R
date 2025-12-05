#############################################
#         Spline-Based Bootstrap            #
#############################################

dr_row_spline <- function(row, method=c("fmm", "natural", "hermite", "periodic"), fmat, complete=TRUE) {
  row_nt <- row[-1] 
  method <- match.arg(method)
  if (method == "hermite") {
    method <- "monoH.FC"
  }
  
  if (sum(!is.na(row_nt)) < 2) {
    warning("Insufficient data to construct spline fit for row(s).")
    return(row_nt) 
  }
  
  row_mask <- !is.na(row_nt) 
  rspline <- splinefun(x = fmat[row_mask], y = row_nt[row_mask], method = method) 
  
  int_row <- rspline(fmat) 
  return(int_row) 
}

dr_row_loglin <- function(row, method=NULL, fmat, complete=TRUE) {
  row_nt <- row[-1] 
  
  if (sum(!is.na(row_nt)) < 2) {
    warning("Insufficient data to construct log-linear fit for row(s).")
    return(row_nt) 
  }
  
  row_mask <- !is.na(row_nt) 
  rloglin <- approxfun(x = fmat[row_mask], y = log(as.numeric(row_nt[row_mask]))) 
  
  int_row <- exp(rloglin(fmat))
  
  return(int_row) 
}
