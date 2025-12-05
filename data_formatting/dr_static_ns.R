#############################################
#           Model-Based Bootstrap           #
#############################################

dr_row_model <- function(row, method=c("nelson", "svensson"), fmat, complete=FALSE) {
  row_nt <- row[-1] 
  method <- match.arg(method)
  if (sum(!is.na(row_nt)) < 2) {
    warning("Insufficient data to construct model-based fit for row(s).")
    return(row_nt) 
  }
  
  row_mask <- !is.na(row_nt) 
  
  if (sum(row_mask) == length(row_nt)){
    return(row_nt) 
  }
  
  build_model <- 
    switch(method, 
           nelson = list(
             fit = function(r, m) Nelson.Siegel(r, m), 
             loadings = function(params, m) NS_loadings(params[4], m),
             factors = function(params) params[1:3] 
           ),
           svensson = list(
             fit = function(r, m) Svensson(r, m), 
             loadings = function(params, m) Svensson_loadings(params[5], params[6], m), 
             factors = function(params) params[1:4] 
           )
    )
  
  pfit <- build_model$fit(row_nt[row_mask], fmat[row_mask]) 
  loadings <- build_model$loadings(pfit, fmat)
  factors <- build_model$factors(pfit) 
  int_row <- factors %*% loadings
  
  if (complete) {
    return(int_row)
  } 
  else {
    o_yield <- row_nt[row_mask] 
    int_row[row_mask] <- o_yield
    return(int_row)
  }
}
