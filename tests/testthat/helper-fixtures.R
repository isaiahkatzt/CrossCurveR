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
