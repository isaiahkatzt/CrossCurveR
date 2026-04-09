Sys.unsetenv("LC_ALL")
Sys.setenv(LANGUAGE = "en")

library(testthat)

helper_file <- sys.frame(1)$ofile
if (is.null(helper_file)) {
  helper_file <- if (file.exists("helpers.R")) "helpers.R" else file.path("unit_tests", "helpers.R")
}
repo_root <- dirname(dirname(normalizePath(helper_file)))
old_wd <- setwd(repo_root)
on.exit(setwd(old_wd), add = TRUE)

source("packages.R")

source_project_file <- function(path) {
  sys.source(path, envir = .GlobalEnv)
}

source_project_file(file.path(repo_root, "core_formatting.R"))
source_project_file(file.path(repo_root, "data_formatting/dr_spline.R"))
source_project_file(file.path(repo_root, "data_formatting/dr_static_ns.R"))
source_project_file(file.path(repo_root, "data_formatting/dr_formatting.R"))
source_project_file(file.path(repo_root, "yield_estimation/curve_reformat.R"))
source_project_file(file.path(repo_root, "yield_estimation/matrix_helpers.R"))
source_project_file(file.path(repo_root, "yield_estimation/cr_core_estim.R"))
source_project_file(file.path(repo_root, "yield_estimation/bl_single_curve.R"))
source_project_file(file.path(repo_root, "yield_estimation/bl_feature_extract.R"))
source_project_file(file.path(repo_root, "yield_estimation/bl_normalize.R"))
source_project_file(file.path(repo_root, "yield_estimation/bl_cross_estimation.R"))

tol <- 1e-6
fixture_mats <- c(1, 3, 6, 12, 24, 60)
fixture_curves <- c("usa", "gbr")
fixture_time <- as.Date("2020-01-01") + 0:19

make_ns_curve_df <- function(curve_shift = 0, lambda = 0.35, mats = fixture_mats, time = fixture_time) {
  idx <- seq_along(time) - 1
  beta_mat <- cbind(
    2 + curve_shift + 0.03 * idx,
    -0.8 + 0.01 * idx,
    0.4 - 0.005 * idx
  )

  yield_mat <- beta_mat %*% NS_loadings(lambda, mats)
  curve_df <- data.frame(time = time, as.data.frame(yield_mat))
  colnames(curve_df) <- c("time", sapply(mats, numeric_to_matname))
  curve_df
}

make_pair_yields <- function() {
  list(
    usa = make_ns_curve_df(curve_shift = 0),
    gbr = make_ns_curve_df(curve_shift = 0.2)
  )
}

make_mock_blsc_list <- function(curves = fixture_curves, mats = fixture_mats, time = fixture_time) {
  curve_dfs <- setNames(
    lapply(seq_along(curves), function(i) make_ns_curve_df(curve_shift = 0.1 * (i - 1), mats = mats, time = time)),
    curves
  )

  curve_phi <- setNames(vector("list", length(curves)), curves)
  curve_betas <- setNames(vector("list", length(curves)), curves)
  curve_nsfit <- setNames(vector("list", length(curves)), curves)
  curve_H <- setNames(vector("list", length(curves)), curves)

  for (curve_name in curves) {
    phi_obj <- bl_phi(0.35, mats)
    betas_obj <- bl_NSbetas(curve_dfs[[curve_name]], 0.35, phi_obj$cross_phi, mats)
    nsfit_obj <- bl_NSfit(betas_obj$m_yield, phi_obj$phi, betas_obj$NSbetas, 0.35, mats)

    curve_phi[[curve_name]] <- phi_obj
    curve_betas[[curve_name]] <- betas_obj
    curve_nsfit[[curve_name]] <- nsfit_obj
    curve_H[[curve_name]] <- diag(length(mats))
  }

  list(
    phi = curve_phi,
    betas = curve_betas,
    nsfit = curve_nsfit,
    H = curve_H
  )
}

make_check_phi <- function(n, mats = fixture_mats, lambda = 0.35) {
  phi_matrix <- t(bl_phi(lambda, mats)$phi)
  rep(list(phi_matrix), n)
}

with_temp_bindings <- function(bindings, expr, env = .GlobalEnv) {
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
