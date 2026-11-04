# Run against an explicitly installed package; no renv activation or auto-install.
local({
  script <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(script)) {
    file_arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (!length(file_arg)) stop("Run this file with Rscript or source().")
    script <- sub("^--file=", "", file_arg[[1]])
  }
  repo_root <- dirname(dirname(normalizePath(script)))
  testthat::test_dir(
    file.path(repo_root, "tests", "testthat"), package = "CrossCurveR",
    load_package = "installed", reporter = "summary"
  )
  testthat::test_dir(
    file.path(repo_root, "unit_tests"), package = "CrossCurveR",
    load_package = "installed", reporter = "summary"
  )
})
