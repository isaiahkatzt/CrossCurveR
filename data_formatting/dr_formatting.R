source("packages.R")
source("core_formatting.R")
source("data_formatting/dr_spline.R") 
source("data_formatting/dr_static_ns.R") 

#############################################
#           Input Data Formatting           #
#############################################

dr_normalize_timestamp <- function(ts, tz = "UTC") {
  if (inherits(ts, "POSIXct") || inherits(ts, "POSIXlt")) {
    return(as.POSIXct(ts, tz = tz))
  }
  if (inherits(ts, "Date")) {
    return(as.POSIXct(ts, tz = tz))
  }
  if (is.numeric(ts)) {
    return(as.POSIXct(ts, origin = "1970-01-01", tz = tz))
  }
  if (is.character(ts)) {
    formats <- c(
      "%Y-%m-%d %H:%M:%S",
      "%Y-%m-%d %H:%M",
      "%Y-%m-%d",
      "%d/%m/%Y %H:%M:%S",
      "%d/%m/%Y %H:%M",
      "%d/%m/%Y",
      "%m/%d/%Y %H:%M:%S",
      "%m/%d/%Y %H:%M",
      "%m/%d/%Y",
      "%Y.%m.%d %H:%M:%S",
      "%Y.%m.%d"
    )
    for (fmt in formats) {
      out <- as.POSIXct(ts, format = fmt, tz = tz)
      if (!any(is.na(out))) return(out)
    }
    stop("Timestamp column could not be parsed. Ensure it is in a standard date/time format.")
  }
  stop("Unsupported timestamp type: ", class(ts))
}

.validate_input_data <- function(data= "", mats = c()) {
  
  time <- dr_normalize_timestamp(data[, 1]) 
  
  if (!(length(data[, -1]) == length(mats))){
    stop("input data columns must match input maturities") 
  }
  
  return(time) 
}

dr_yield_format <- function(data = "", mats = c(), start_date = "", end_date = "") {
  ## format yields for later use 
  time <- .validate_input_data(data = data, mats = mats)
  mat_names <- sapply(mats, numeric_to_matname)
  data[, 1] <- time
  
  colnames(data) <- c("time", mat_names)

  data <- data %>%
    dplyr::filter(time >= as.Date(start_date) & time <= as.Date(end_date)) %>% 
    mutate(time = strtrim(time, 10)) %>% 
    mutate(across(everything(),
                  ~ ifelse(is.nan(.), NA, .)))
  return(data)
}

#############################################
#      Single-Day Missing Interpolation     #
#############################################

dr_daily <- function(data, reference=NULL, tkey="time", threshold = 2){
  ## linear interpolation for sub-thresholded vals  
  if (!is.null(reference)) {
    ref_match <- match(reference, as.Date(data[[tkey]]))
    augmented <- matrix(NA, nrow = length(reference), ncol = ncol(data))
    augmented[!is.na(ref_match), ] <- as.matrix(data[ref_match[!is.na(ref_match)], ])
    augmented <- data.frame(augmented)
    colnames(augmented) <- colnames(data)
    augmented[[tkey]] <- reference
  } else {
    augmented <- data
  }
  
  num_cols <- setdiff(colnames(augmented), tkey)
  augmented[num_cols] <- lapply(augmented[num_cols], as.numeric)
  
  yield_matrix <- as.matrix(augmented[num_cols])
  nr <- nrow(yield_matrix)
  nc <- ncol(yield_matrix)
  
  sparse_yields <- rowSums(!is.na(yield_matrix)) < (threshold + 1)
  preceding  <- rbind(matrix(NA, 1, nc), yield_matrix[-nr, ])
  succeeding <- rbind(yield_matrix[-1, ], matrix(NA, 1, nc))
  
  int_mask <- sparse_yields & is.na(yield_matrix) & 
    !is.na(preceding) & !is.na(succeeding)
  
  for (col_id in seq_len(nc)) {
    row_mask <- int_mask[, col_id]
    if (any(row_mask)) {
      yprev <- preceding[row_mask,  col_id]
      ynext <- succeeding[row_mask, col_id]
      yield_matrix[row_mask, col_id] <- (yprev + ynext) / 2
    }
  }
  
  augmented[num_cols] <- yield_matrix
  
  return(augmented)
}

#############################################
#               Curve Builder               #
#############################################

dr_build_curve <- function(data, nmats=NULL, tkey="time", build=c("spline", "log", "parametric"), 
                           par_method = c("nelson", "svensson"), 
                           spline_method = c("fmm", "natural", "hermite", "periodic"), 
                           complete = FALSE,
                           precision = 4) {
  
  mats <- sapply(colnames(data)[-1], matname_to_numeric, USE.NAMES=FALSE) * 12
  dmats <- setdiff(nmats, mats) 
  build = match.arg(build) 
  
  for (m in dmats) {
    mstr <- numeric_to_matname(m) 
    data[[mstr]] <- NA 
  }
  
  fmats <- sort(c(mats, dmats))
  fstr <- sapply(fmats, numeric_to_matname, USE.NAMES=FALSE) 
  ordered_data <- data[c(tkey, fstr)] 
  
  build_setup <- switch(
    build, 
    spline = list(method = match.arg(spline_method), 
                  build_func = dr_row_spline), 
    parametric = list(method = match.arg(par_method), 
                      build_func = dr_row_model), 
    log = list(method=NULL, build_func = dr_row_loglin))
  
  int_data <- apply(ordered_data, 1, build_setup$build_func, build_setup$method, fmats, complete)
  
  
  int_df <- data.frame(ordered_data[[tkey]], t(int_data)) %>% 
    mutate(across(-1, ~ round(as.numeric(.x), digits = precision))) 
  colnames(int_df) <- colnames(ordered_data) 
  
  return(int_df) 
}

