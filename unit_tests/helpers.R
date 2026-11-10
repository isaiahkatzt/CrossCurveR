library(CrossCurveR)
withr::local_envvar(c(LC_ALL = NA), .local_envir = testthat::teardown_env())

tol <- 1e-6
repo_root <- dirname(normalizePath(testthat::test_path()))
source(file.path(repo_root, "applications", "tontine_setup.R"), local = TRUE)
source(file.path(repo_root, "applications", "tontine_helpers.R"), local = TRUE)
source(file.path(repo_root, "tests", "testthat", "helper-feature-reuse.R"), local = TRUE)

with_temp_bindings <- function(bindings, expr, env = environment(tontine_fit_baseline_curve_models)) {
  expr <- substitute(expr)
  old_values <- list()
  had_value <- logical(length(bindings))
  names(had_value) <- names(bindings)

  for (nm in names(bindings)) {
    had_value[[nm]] <- exists(nm, envir = env, inherits = FALSE)
    if (had_value[[nm]]) {
      old_values[[nm]] <- get(nm, envir = env, inherits = FALSE)
    }
    assign(nm, bindings[[nm]], envir = env)
  }

  on.exit({
    for (nm in names(bindings)) {
      if (had_value[[nm]]) {
        assign(nm, old_values[[nm]], envir = env)
      } else if (exists(nm, envir = env, inherits = FALSE)) {
        rm(list = nm, envir = env)
      }
    }
  }, add = TRUE)

  eval(expr, envir = parent.frame())
}
