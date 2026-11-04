dr_normalize_maturity_names <- function(data, tkey = "time") {
  if (!is.data.frame(data) || ncol(data) < 2L) {
    stop("`data` must be a data frame with a time column and at least one maturity column.",
         call. = FALSE)
  }
  time_index <- which(names(data) == tkey)
  if (length(time_index) != 1L) {
    stop("`tkey` must identify exactly one time column in `data`.", call. = FALSE)
  }

  maturity_indices <- setdiff(seq_len(ncol(data)), time_index)
  labels <- names(data)[maturity_indices]
  months <- vapply(labels, .matname_to_months, numeric(1), USE.NAMES = FALSE)
  canonical <- vapply(months, numeric_to_matname, character(1), USE.NAMES = FALSE)

  if (anyDuplicated(canonical)) {
    stop("Maturity columns identify duplicate tenors after normalization.",
         call. = FALSE)
  }

  # Only the time column is moved; maturity values retain their column order.
  data <- data[c(time_index, maturity_indices)]
  names(data) <- c(tkey, paste0("X", canonical))
  data
}

dr_read_yields <- function(file, tkey = "time") {
  # Preserve original headers so invalid or duplicate names cannot be hidden
  # by read.csv() name repair before validation.
  data <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  data <- dr_normalize_maturity_names(data, tkey = tkey)
  names(data)[1] <- "time"
  data
}
