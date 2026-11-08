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

.matname_to_months <- function(colname) {
  if (!is.character(colname) || length(colname) != 1L || is.na(colname)) {
    stop("Maturity name must be a single non-missing character string.", call. = FALSE)
  }

  label <- toupper(trimws(colname))
  if (!grepl("^X?[0-9]+(\\.[0-9]+)?[MY]$", label)) {
    stop(sprintf(
      'Unrecognized maturity name "%s". Use labels such as "01M", "12M", "01Y", or "X01Y".',
      colname
    ), call. = FALSE)
  }

  label <- sub("^X", "", label)
  amount <- as.numeric(sub("[MY]$", "", label))
  months <- if (endsWith(label, "M")) amount else amount * 12
  if (!is.finite(months) || months <= 0) {
    stop(sprintf('Maturity name "%s" must represent a positive, finite maturity.',
                 colname), call. = FALSE)
  }

  # Undo representation error when a whole month was written in years.
  nearest_month <- round(months)
  if (nearest_month > 0 &&
      abs(months - nearest_month) <= 16 * .Machine$double.eps * max(1, months)) {
    months <- nearest_month
  }
  months
}

matname_to_numeric <- function(colname) {
  # Preserve the public helper's return unit: years.
  .matname_to_months(colname) / 12
}

#############################################
#       Cross-Module Parametric Fits        #
#############################################

NS_loadings_ls <- function(lambda, mats){
  .cc_mats(mats)
  .cc_decay(lambda, "`lambda`", scalar = TRUE)
  l1 <- rep(1, length(mats))
  l2 <- (1 - exp(-mats * lambda)) / (lambda * mats)
  rbind(l1, l2)
}

NS_loadings <- function(lambda, mats){
  .cc_mats(mats)
  .cc_decay(lambda, "`lambda`", scalar = TRUE)
  l1 <- rep(1, length(mats))
  l2 <- (1 - exp(-mats * lambda)) / (lambda * mats)
  l3 <- (1 - exp(-mats * lambda)) / (lambda * mats) - exp(-lambda * mats)
  rbind(l1, l2, l3)
}

Svensson_loadings <- function(lambda1, lambda2, mats){
  .cc_mats(mats)
  .cc_decay(lambda1, "`lambda1`", scalar = TRUE)
  .cc_decay(lambda2, "`lambda2`", scalar = TRUE)
  l1 <- rep(1, length(mats))
  l2 <- (1 - exp(-mats / lambda1)) / (mats / lambda1)
  l3 <- ((1 - exp(-mats / lambda1)) / (mats / lambda1)) - exp(-mats / lambda1)
  l4 <- ((1 - exp(-mats / lambda2)) / (mats / lambda2)) - exp(-mats / lambda2)
  rbind(l1, l2, l3, l4)
}
