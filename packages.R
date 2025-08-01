packages <- c(
  "tidyverse", 
  "dplyr",
  "plotly", 
  "reshape2",
  "ggplot2",
  "Matrix",
  "vars",
  "urca",
  "zoo", 
  "MASS", 
  "splines",
  "tsDyn",
  "parallel",
  "purrr",
  "twosamples",
  "CovTools", 
  "Hotelling", 
  "stats", 
  "lubridate", 
  "tseries", 
  "YieldCurve", 
  "dgof", 
  "matrixStats"
)

install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

invisible(lapply(packages, install_and_load))