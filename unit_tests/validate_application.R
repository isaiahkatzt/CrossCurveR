# Capture or compare a reproducible application subset using the installed package.
# Usage: Rscript --vanilla unit_tests/validate_application.R record|check reference.rds
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[[1]] %in% c("record", "check")) {
  stop("Usage: validate_application.R record|check reference.rds", call. = FALSE)
}
mode <- args[[1]]
reference_file <- args[[2]]
if (mode == "record" && file.exists(reference_file)) {
  stop("Reference already exists; choose a new path to preserve the baseline.", call. = FALSE)
}
reference <- if (mode == "check") readRDS(reference_file) else NULL
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1]])
repo <- dirname(dirname(normalizePath(script)))
setwd(repo)
library(CrossCurveR)

# Execute the notebook's actual computational chunks, without copying its model code.
# Restrict the replay to 2015–2016 and the 2015 global level shock.
lines <- readLines(file.path(repo, "applications", "application_tontine_1yr.Rmd"))
selected <- c("libraries", "configuration", "data", "valuation-windows",
  "tontine-event-path", "baseline-models", "return-factor-helpers",
  "baseline-valuation", "shock-specs", "shock-valuations", "shock-summaries",
  "country-scenario-bond-pl-display")
e <- new.env(parent = globalenv())
executed <- character()
for (first in grep("^```\\{r ", lines)) {
  label <- sub("^```\\{r ([^,}]+).*", "\\1", lines[[first]])
  if (!label %in% selected) next
  last <- which(seq_along(lines) > first & lines == "```")[[1]]
  invisible(eval(parse(text = lines[seq.int(first + 1L, last - 1L)]), e))
  if (label == "configuration") {
    e$purchase_years <- 2015:2016
    e$q <- e$q_full[match(e$purchase_years, 2011:2020)]
  }
  if (label == "shock-specs") e$shock_examples <- e$shock_examples["global_l_neg30pct_2015"]
  executed <- c(executed, label)
}
stopifnot(identical(executed, selected),
  all(is.finite(e$baseline_rung_summary$average_position_end)),
  all(is.finite(e$shock_rung_summary$average_position_end)),
  all(is.finite(e$baseline_return_factors$gross_return_factor)),
  all(is.finite(e$shock_return_factors$gross_return_factor)))
diagnostics <- lapply(e$baseline_curve_models$MCE$fit, function(window) {
  lapply(window$fit$BS0, function(fit) fit[c("iter", "converged")])
})
data_hashes <- tools::md5sum(file.path(repo, "reconstructed_data", unlist(e$curve_files)))
names(data_hashes) <- names(e$curve_files)
results <- list(
  configuration = mget(c("purchase_years", "curves", "mats", "lambdas", "cutoffs",
    "fit_seed", "event_seed", "mce_baseline_controls", "mce_shock_controls",
    "N_members", "initial_position_per_member", "q", "eta_params",
    "withdrawal_fraction", "withdrawal_penalty_rate"), e),
  data_hashes = data_hashes,
  tables = Filter(is.data.frame, as.list(e)),
  baseline_panels = lapply(e$baseline_curve_models, `[[`, "curve_panels"),
  shock_panels = lapply(e$shock_curve_models, function(models) lapply(models, `[[`, "curve_panels")),
  diagnostics = diagnostics
)
artifact <- list(results = results, provenance = list(
  package_path = find.package("CrossCurveR"),
  package_version = as.character(utils::packageVersion("CrossCurveR")),
  notebook_hash = tools::md5sum(file.path(repo, "applications", "application_tontine_1yr.Rmd")),
  session = capture.output(sessionInfo())))
if (mode == "record") {
  saveRDS(artifact, reference_file)
  cat("Application reference saved:", reference_file, "\n")
} else {
  saveRDS(artifact, paste0(reference_file, ".candidate.rds"))
  # Structural refactors must preserve data, configuration, and convergence decisions.
  stopifnot(identical(results$data_hashes, reference$results$data_hashes),
    identical(results$configuration, reference$results$configuration),
    identical(results$diagnostics, reference$results$diagnostics))
  comparison <- all.equal(reference$results, results, tolerance = 1e-10, check.attributes = TRUE)
  if (!isTRUE(comparison)) {
    stop(paste(c("Application regression comparison failed:", comparison), collapse = "\n"), call. = FALSE)
  }
  cat("Application comparison passed (tolerance 1e-10); candidate saved alongside reference.\n")
}
