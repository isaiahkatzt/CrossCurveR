source("packages.R")
source("core_formatting.R")

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
  ## reformat yields 
  
  time <- .validate_input_data(data = data, mats = mats)
  mat_names <- sapply(mats, numeric_to_matname)
  data[, 1] <- time
  
  colnames(data) <- c("time", mat_names)

  data <- data %>%
    dplyr::filter(time >= as.Date(start_date) & time <= as.Date(end_date))
  return(data)
}

