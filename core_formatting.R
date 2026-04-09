#############################################
#         Cross-Module Formatting           #
#############################################

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

matname_to_numeric <- function(colname){
  ## convert column (string) maturities to numerical 
  num <- as.numeric(sub("[A-Z]+", "", colname))
  time_unit <- sub("[0-9]+", "", colname)
  
  if (time_unit == "M") { return(num/12) }
  else { return(num) }
}

#############################################
#       Cross-Module Parametric Fits        #
#############################################

NS_loadings_ls <- function(lambda, mats){
  l1 <- rep(1, length(mats)) 
  l2 <- (1 - exp(-mats * lambda)) / (lambda * mats) 
  rbind(l1, l2) 
}

NS_loadings <- function(lambda, mats){
  l1 <- rep(1, length(mats))
  l2 <- (1 - exp(-mats * lambda)) / (lambda * mats) 
  l3 <- (1 - exp(-mats * lambda)) / (lambda * mats) - exp(-lambda * mats)
  rbind(l1, l2, l3) 
}

Svensson_loadings <- function(lambda1, lambda2, mats){
  l1 <- rep(1, length(mats)) 
  l2 <- (1 - exp(-mats / lambda1)) / (mats / lambda1) 
  l3 <- ((1 - exp(-mats / lambda1)) / (mats / lambda1)) - exp(-mats / lambda1) 
  l4 <- ((1 - exp(-mats / lambda2)) / (mats / lambda2)) - exp(-mats / lambda2)
  rbind(l1, l2, l3, l4) 
}
