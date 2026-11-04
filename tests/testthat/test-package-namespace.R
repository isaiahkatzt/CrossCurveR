test_that("the public API resolves to the package namespace", {
  exports <- getNamespaceExports("CrossCurveR")
  expect_true(all(c("dly_fit", "mc_fit_end", "mc_fit_exo", "mc_fit_shock",
                    "dr_build_curve") %in% exports))
  expect_identical(environment(CrossCurveR::dly_fit), asNamespace("CrossCurveR"))
  expect_identical(environment(CrossCurveR::mc_fit_exo), asNamespace("CrossCurveR"))
  expect_false("dynamic_mean_covreg" %in% exports)
})

test_that("an installed package works in a separate base-only R session", {
  # Rscript starts outside the repository, without profiles or attached helpers.
  # Preserve the installed/check library path when launching the child process.
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", "--default-packages=NULL",
      shQuote(test_path("fixtures", "namespace-integration.R"))),
    stdout = TRUE, stderr = TRUE,
    env = paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep)))
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  expect_identical(status, 0L, info = paste(output, collapse = "\n"))
  expect_true(any(grepl("Namespace integration check passed", output, fixed = TRUE)))
})
