#############################################
#          Single-Curve Estimation          #
#############################################

dly_phi <- function(lambda_optim, mats) {
  phi <- NS_loadings_ls(lambda_optim, mats)
  cross_phi <- .bl_loading_inverse(Matrix::tcrossprod(phi),
    sprintf("DLY country loading regression (lambda=%g; mats=[%s] months)",
            lambda_optim, paste(mats, collapse = ", "))) %*% phi
  return(list(phi = phi, cross_phi = cross_phi))
}

dly_sc_fit <- function(lambda_select, yields, mats) {
  yields <- .cc_panel(yields, mats)
  curve_phi <- dly_phi(lambda_select, mats)
  curve_betas <- bl_NSbetas(yields, lambda_select, curve_phi$cross_phi, mats)
  curve_nsfit <- bl_NSfit(
    curve_betas$m_yield,
    curve_phi$phi,
    curve_betas$NSbetas,
    lambda_select,
    mats
  )

  return(list(
    lambda = lambda_select,
    time = yields[[1]],
    yield_names = colnames(yields)[-1],
    phi = curve_phi,
    betas = curve_betas,
    nsfit = curve_nsfit
  ))
}

dly_sc_fit_all <- function(lambdas, yields, mats, curves = names(yields)) {
  if (is.null(curves)) {
    curves <- paste0("curve", seq_along(yields))
  }

  if (length(curves) != length(yields)) {
    stop("curves must have the same length as yields.")
  }
  yields <- .cc_yields(yields, mats, curves)
  .cc_decay(lambdas)

  curve_lambdas <- vapply(yields, function(x) {
    lambda_grid_search(lambdas, x, mats, loss = l2_loss_ls)
  }, numeric(1))

  curve_NSparam <- stats::setNames(lapply(seq_along(yields), function(i) {
    dly_sc_fit(curve_lambdas[[i]], yields[[i]], mats)
  }), curves)

  return(curve_NSparam)
}

#############################################
#         Country Factor Extraction         #
#############################################

dly_get_betas <- function(dly_sclist, curves = names(dly_sclist)) {
  if (is.null(curves)) {
    curves <- names(dly_sclist)
  }

  stats::setNames(lapply(curves, function(curve) {
    dly_sclist[[curve]]$betas$NSbetas
  }), curves)
}

dly_factor_panel <- function(sc_full, curves = names(sc_full), latent = c("L", "S")) {
  latent <- match.arg(latent, c("L", "S"))
  latent_idx <- ifelse(latent == "L", 1, 2)
  curve_beta <- dly_get_betas(sc_full, curves)

  factor_draw <- do.call(cbind, lapply(curve_beta, function(beta_mat) {
    beta_mat[, latent_idx]
  }))

  colnames(factor_draw) <- curves
  return(factor_draw)
}

#############################################
#           Compatibility Aliases           #
#############################################
dly_single_curve <- dly_sc_fit
dly_all_sc <- dly_sc_fit_all
dly_pull_beta <- dly_get_betas
dly_group_factor <- dly_factor_panel
