#### core Nelson-Siegel parameters  

NS_loadings <- function(lambda, M){
  l1 <- rep(1, length(M))
  l2 <- (1 - exp(-M * lambda)) / (lambda * M) 
  l3 <- (1 - exp(-M * lambda)) / (lambda * M) - exp(-lambda * M)
  rbind(l1, l2, l3) 
}
