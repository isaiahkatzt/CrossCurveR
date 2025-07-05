source("packages.R")

#### baseline single-curve dns estimatoin  

l2_loss_lambda <- function(lambda_grid, daily_yield, mats) {
  n <- length(lambda_grid)
  loss_row <- numeric(n)  
  factor_loadings <- lapply(lambda_grid, NS_loadings, mats)
  
  for (i in seq_along(lambda_grid)) {
    X <- t(factor_loadings[[i]])  
    fit <- lm.fit(X, daily_yield)  
    loss_row[i] <- sum(fit$residuals^2)  
  }
  loss_row
}

lambda_grid_search <- function(lambda_grid, yields, mats) {
  ## fast search lambda grid 
  yields <- as.matrix(yields[, -1, drop = FALSE])  # Ensure matrix structure
  N <- nrow(yields)
  n <- length(lambda_grid)
  
  ## vectorized loss computation 
  loss_matrix <- t(vapply(seq_len(N), function(i) l2_loss_lambda(lambda_grid, yields[i, ], mats), numeric(n)))
  
  ## optimal lambda selection 
  optimal_lambda <- lambda_grid[which.min(colSums(loss_matrix))]
  return(optimal_lambda)
}

