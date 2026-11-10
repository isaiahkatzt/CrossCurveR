test_that("the public API resolves to the package namespace", {
  exports <- getNamespaceExports("CrossCurveR")
  expect_true(all(c("dly_fit", "mc_fit_end", "mc_fit_exo", "mc_fit_shock",
                    "dr_build_curve", "sc_fit_shock") %in% exports))
  expect_true(all(c("NS_loadings", "NS_loadings_ls", "Svensson_loadings",
                    "bl_single_curve", "zero_mean_covreg", "zero_mean_covreg_em",
                    "cr_build", "cr_synthetic", "st_factor_profile", "st_mean_profile",
                    "st_dly_shock", "dly_sc_fit_all", "dly_fit_from_sc") %in% exports))
  expect_identical(environment(CrossCurveR::dly_fit), asNamespace("CrossCurveR"))
  expect_identical(environment(CrossCurveR::mc_fit_exo), asNamespace("CrossCurveR"))
})

test_that("estimation and shock plumbing remain internal", {
  helpers <- c(
    "vecY", "PY_full", "cf_curve_cut", "cf_phi_cut", "bln_build_P",
    "sqrt_inv_build", "bln_YB", "bln_phi_hat", "bln_phi_BT", "bln_phi_CT",
    "bln_normalized_fit", "blce_rescale_yield", "bl_phi", "bl_NSbetas",
    "blmc_VAR", "blmc_ECM", "blcc_fe",
    "st_sc_state", "st_prepare_single_curve_state",
    "st_xt_state", "st_build_xt_state",
    "st_shock_state", "st_build_shocked_curve_state",
    "st_project_fit", "st_finalize_shocked_fit",
    "st_align_Xt", "st_align_Xt_to_baseline",
    "st_check_baseline", "st_validate_baseline_state", "st_resolve_config",
    "sc_NSdaily", "sc_NSwindow", "numeric_to_matname", "matname_to_numeric",
    "dr_normalize_timestamp", "lambda_grid_search"
  )
  expect_length(intersect(helpers, getNamespaceExports("CrossCurveR")), 0L)
  ns <- asNamespace("CrossCurveR")
  expect_true(all(vapply(helpers, function(name) {
    is.function(get(name, ns, inherits = FALSE))
  }, logical(1))))
  expect_error(getExportedValue("CrossCurveR", "bl_phi"), "not an exported object")
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
