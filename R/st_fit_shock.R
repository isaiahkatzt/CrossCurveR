#############################################
#         DNS and MCE Shock Fitting         #
#############################################

mc_fit_shock_base <- function(yields, lambdas, cutoffs, reference,
                              ECM_estim = "ML", ECM_type = "eigen", ECM_alpha = 0.1,
                              X_normalize = TRUE, X_trunc = FALSE,
                              CR_algo = zero_mean_covreg_em, CR_init = "adaptive",
                              Xt_smooth = FALSE, smoother = c("ns", "bs", "rm", "hp", "henderson"),
                              knot_count = 5, k_count = 9, k_pass = 3,
                              hp_lambda = 1600, henderson_k = 13,
                              CR_maxiter = 1000, CR_tol = 1e-8,
                              CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = 0,
                              CR_check_every = 1, CR_use_loglik = TRUE, CR_store_path = TRUE,
                              CR_verb = FALSE, CR_term = TRUE,
                              curves, mats) {
  yields <- st_check_yields(yields, curves, reference)
  yields <- .cc_yields(yields, mats, curves, reference)

  single_curve_state <- st_sc_state(
    yields = yields,
    lambdas = lambdas,
    cutoffs = cutoffs,
    curves = curves,
    mats = mats
  )

  baseline_xt_state <- st_xt_state(
    bl_yields = single_curve_state$bl_yields,
    curves = curves,
    reference = reference,
    reuse_ECM = FALSE,
    baseline_ecm = NULL,
    ECM_estim = ECM_estim,
    ECM_type = ECM_type,
    ECM_alpha = ECM_alpha,
    X_normalize = X_normalize,
    X_trunc = X_trunc,
    Xt_smooth = Xt_smooth,
    smoother = smoother,
    knot_count = knot_count,
    k_count = k_count,
    k_pass = k_pass,
    hp_lambda = hp_lambda,
    henderson_k = henderson_k,
    ytime = single_curve_state$time
  )

  baseline_BS0 <- blcc_fe(
    Xt = baseline_xt_state$Xt,
    Wt = single_curve_state$W,
    mats = mats,
    covreg = CR_algo,
    init = CR_init,
    max_iter = CR_maxiter,
    tol = CR_tol,
    S0 = CR_S0init,
    B = CR_Binit,
    S0_shrink_diag = CR_S0_shrink_diag,
    check_every = CR_check_every,
    use_loglik = CR_use_loglik,
    store_path = CR_store_path,
    verb = CR_verb,
    term = CR_term
  )

  baseline_fit_state <- st_project_fit(
    yields = single_curve_state$yields,
    time = single_curve_state$time,
    phi_hat = single_curve_state$phi_hat,
    Xt = baseline_xt_state$Xt,
    BS0 = baseline_BS0,
    curves = curves,
    mats = mats,
    W = single_curve_state$W
  )

  structure(list(
    curves = curves,
    mats = mats,
    reference = reference,
    yields = yields,
    single_curve_state = single_curve_state,
    xt_state = baseline_xt_state,
    BS0 = baseline_BS0,
    baseline_fit = baseline_fit_state,
    input_spec = .cc_input_spec(yields, mats, curves, baseline_xt_state$Xt,
                                baseline_xt_state$feature_spec,
                                decay = single_curve_state$input_spec$decay),
    config = list(
      lambdas = lambdas,
      cutoffs = cutoffs,
      ECM_estim = ECM_estim,
      ECM_type = ECM_type,
      ECM_alpha = ECM_alpha,
      X_normalize = X_normalize,
      X_trunc = X_trunc,
      CR_algo = CR_algo,
      CR_init = CR_init,
      Xt_smooth = Xt_smooth,
      smoother = smoother,
      knot_count = knot_count,
      k_count = k_count,
      k_pass = k_pass,
      hp_lambda = hp_lambda,
      henderson_k = henderson_k,
      CR_maxiter = CR_maxiter,
      CR_tol = CR_tol,
      CR_Binit = CR_Binit,
      CR_S0init = CR_S0init,
      CR_S0_shrink_diag = CR_S0_shrink_diag,
      CR_check_every = CR_check_every,
      CR_use_loglik = CR_use_loglik,
      CR_store_path = CR_store_path,
      CR_verb = CR_verb,
      CR_term = CR_term
    )
  ), class = "mc_shock_baseline")
}

mc_fit_shocked <- function(baseline_state, shocked_factors,
                           reuse_BS0 = TRUE, reuse_ECM = TRUE,
                           ECM_fallback = c("error", "baseline"), controls = list()) {
  ECM_fallback <- match.arg(ECM_fallback)
  baseline_state <- st_check_baseline(baseline_state)
  single_curve_state <- baseline_state$single_curve_state
  config <- baseline_state$config
  control_names <- intersect(names(formals(mc_fit_shock_base)),
                             names(formals(mc_fit_shock_cached)))
  controls <- stats::setNames(lapply(control_names, function(name) {
    st_resolve_config(controls[[name]], config[[name]])
  }), control_names)

  shocked_state <- st_shock_state(
    single_curve_state = single_curve_state,
    shocked_factors = shocked_factors
  )

  shocked_xt_state <- st_xt_state(
    bl_yields = shocked_state$bl_yields,
    curves = baseline_state$curves,
    reference = baseline_state$reference,
    reuse_ECM = reuse_ECM,
    baseline_ecm = baseline_state$xt_state$cc_ci_features$cc_ecm,
    ECM_estim = controls$ECM_estim,
    ECM_type = controls$ECM_type,
    ECM_alpha = controls$ECM_alpha,
    ECM_fallback = ECM_fallback,
    X_normalize = controls$X_normalize,
    X_trunc = controls$X_trunc,
    Xt_smooth = controls$Xt_smooth,
    smoother = controls$smoother,
    knot_count = controls$knot_count,
    k_count = controls$k_count,
    k_pass = controls$k_pass,
    hp_lambda = controls$hp_lambda,
    henderson_k = controls$henderson_k,
    ytime = single_curve_state$time
  )

  if (reuse_BS0) {
    shocked_Xt <- st_align_Xt(shocked_xt_state$Xt, baseline_state$xt_state$Xt)
    shocked_BS0 <- baseline_state$BS0
  } else {
    shocked_Xt <- shocked_xt_state$Xt
    shocked_BS0 <- blcc_fe(
      Xt = shocked_Xt,
      Wt = shocked_state$W,
      mats = baseline_state$mats,
      covreg = controls$CR_algo,
      init = controls$CR_init,
      max_iter = controls$CR_maxiter,
      tol = controls$CR_tol,
      S0 = controls$CR_S0init,
      B = controls$CR_Binit,
      S0_shrink_diag = controls$CR_S0_shrink_diag,
      check_every = controls$CR_check_every,
      use_loglik = controls$CR_use_loglik,
      store_path = controls$CR_store_path,
      verb = controls$CR_verb,
      term = controls$CR_term
    )
  }

  shocked_fit_state <- st_project_fit(
    yields = shocked_state$yields,
    time = single_curve_state$time,
    phi_hat = shocked_state$phi_hat,
    Xt = shocked_Xt,
    BS0 = shocked_BS0,
    curves = baseline_state$curves,
    mats = baseline_state$mats,
    W = shocked_state$W
  )

  return(list(
    tenor = baseline_state$baseline_fit$tenor,
    curve = baseline_state$baseline_fit$curve,
    sigma_JT = shocked_fit_state$sigma_JT,
    Xt = shocked_fit_state$Xt,
    W = shocked_fit_state$W,
    ns_factor = shocked_state$ns_factor,
    shocked_factors = shocked_factors,
    shocked_yields = shocked_state$yields,
    shocked_tenor = shocked_fit_state$tenor,
    shocked_curve = shocked_fit_state$curve,
    shocked_single_curve_state = shocked_state,
    baseline_state = baseline_state,
    BS0 = shocked_BS0,
    xt_state = shocked_xt_state,
    input_spec = .cc_input_spec(shocked_state$yields, baseline_state$mats,
      baseline_state$curves, shocked_Xt, shocked_xt_state$feature_spec)
  ))
}

mc_fit_shock_cached <- function(baseline_state,
                                shock_curves, shock_factors, shock_magnitude,
                                shock_type = c("additive", "multiplicative"),
                                shock_window = NULL,
                                shock_window_position = c("first", "last"),
                                reuse_BS0 = TRUE, reuse_ECM = TRUE,
                                ECM_estim = NULL, ECM_type = NULL, ECM_alpha = NULL,
                                ECM_fallback = c("error", "baseline"),
                                X_normalize = NULL, X_trunc = NULL,
                                CR_algo = NULL, CR_init = NULL,
                                Xt_smooth = NULL, smoother = NULL,
                                knot_count = NULL, k_count = NULL, k_pass = NULL,
                                hp_lambda = NULL, henderson_k = NULL,
                                CR_maxiter = NULL, CR_tol = NULL,
                                CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = NULL,
                                CR_check_every = NULL, CR_use_loglik = NULL, CR_store_path = NULL,
                                CR_verb = NULL, CR_term = NULL) {
  shock_type <- match.arg(shock_type)
  shock_window_position <- match.arg(shock_window_position)
  ECM_fallback <- match.arg(ECM_fallback)
  baseline_state <- st_check_baseline(baseline_state)
  single_curve_state <- baseline_state$single_curve_state
  control_names <- intersect(names(formals(mc_fit_shock_base)),
                             names(formals(mc_fit_shock_cached)))
  controls <- mget(control_names, envir = environment(), inherits = FALSE)

  shock_fit <- st_factor_shock(
    curves = baseline_state$curves,
    curve_dns_factor = single_curve_state$ns_factor,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window,
    shock_window_position = shock_window_position
  )

  fit <- mc_fit_shocked(
    baseline_state, shock_fit$shocked_factors, reuse_BS0, reuse_ECM, ECM_fallback, controls
  )
  append(fit, list(shock_matrix = shock_fit$shock_matrix), after = match("ns_factor", names(fit)))
}

mc_fit_shock_refit <- function(baseline_state,
                               shock_curves, shock_factors, shock_magnitude,
                               shock_type = c("additive", "multiplicative"),
                               shock_window = NULL,
                               shock_window_position = c("first", "last"),
                               ECM_estim = NULL, ECM_type = NULL, ECM_alpha = NULL,
                               ECM_fallback = c("error", "baseline"),
                               X_normalize = NULL, X_trunc = NULL,
                               CR_algo = NULL, CR_init = NULL,
                               Xt_smooth = NULL, smoother = NULL,
                               knot_count = NULL, k_count = NULL, k_pass = NULL,
                               hp_lambda = NULL, henderson_k = NULL,
                               CR_maxiter = NULL, CR_tol = NULL,
                               CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = NULL,
                               CR_check_every = NULL, CR_use_loglik = NULL, CR_store_path = NULL,
                               CR_verb = NULL, CR_term = NULL) {
  shock_type <- match.arg(shock_type)
  shock_window_position <- match.arg(shock_window_position)
  ECM_fallback <- match.arg(ECM_fallback)

  mc_fit_shock_cached(
    baseline_state = baseline_state,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window,
    shock_window_position = shock_window_position,
    reuse_BS0 = FALSE,
    reuse_ECM = FALSE,
    ECM_estim = ECM_estim,
    ECM_type = ECM_type,
    ECM_alpha = ECM_alpha,
    ECM_fallback = ECM_fallback,
    X_normalize = X_normalize,
    X_trunc = X_trunc,
    CR_algo = CR_algo,
    CR_init = CR_init,
    Xt_smooth = Xt_smooth,
    smoother = smoother,
    knot_count = knot_count,
    k_count = k_count,
    k_pass = k_pass,
    hp_lambda = hp_lambda,
    henderson_k = henderson_k,
    CR_maxiter = CR_maxiter,
    CR_tol = CR_tol,
    CR_Binit = CR_Binit,
    CR_S0init = CR_S0init,
    CR_S0_shrink_diag = CR_S0_shrink_diag,
    CR_check_every = CR_check_every,
    CR_use_loglik = CR_use_loglik,
    CR_store_path = CR_store_path,
    CR_verb = CR_verb,
    CR_term = CR_term
  )
}

mc_fit_shock <- function(yields = NULL, lambdas = NULL, cutoffs = NULL, reference = NULL,
                         shock_curves, shock_factors, shock_magnitude,
                         shock_type = c("additive", "multiplicative"),
                         shock_window = NULL,
                         shock_window_position = c("first", "last"),
                         reuse_BS0 = TRUE, reuse_ECM = TRUE,
                         ECM_estim = "ML", ECM_type = "eigen", ECM_alpha = 0.1,
                         ECM_fallback = c("error", "baseline"),
                         X_normalize = TRUE, X_trunc = FALSE,
                         CR_algo = zero_mean_covreg_em, CR_init = "adaptive",
                         Xt_smooth = FALSE, smoother = c("ns", "bs", "rm", "hp", "henderson"),
                         knot_count = 5, k_count = 9, k_pass = 3,
                         hp_lambda = 1600, henderson_k = 13,
                         CR_maxiter = 1000, CR_tol = 1e-8,
                         CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = 0,
                         CR_check_every = 1, CR_use_loglik = TRUE, CR_store_path = TRUE,
                         CR_verb = FALSE, CR_term = TRUE,
                         curves = NULL, mats = NULL, baseline_state = NULL) {
  control_names <- intersect(names(formals(mc_fit_shock_base)),
                             names(formals(mc_fit_shock_cached)))
  supplied_controls <- intersect(names(match.call())[-1], control_names)
  ECM_fallback <- match.arg(ECM_fallback)

  if (is.null(baseline_state)) {
    if (is.null(yields) || is.null(lambdas) || is.null(cutoffs) ||
        is.null(reference) || is.null(curves) || is.null(mats)) {
      stop(
        "When `baseline_state` is not supplied, `yields`, `lambdas`, `cutoffs`, ",
        "`reference`, `curves`, and `mats` must be provided."
      )
    }

    baseline_state <- mc_fit_shock_base(
      yields = yields,
      lambdas = lambdas,
      cutoffs = cutoffs,
      reference = reference,
      ECM_estim = ECM_estim,
      ECM_type = ECM_type,
      ECM_alpha = ECM_alpha,
      X_normalize = X_normalize,
      X_trunc = X_trunc,
      CR_algo = CR_algo,
      CR_init = CR_init,
      Xt_smooth = Xt_smooth,
      smoother = smoother,
      knot_count = knot_count,
      k_count = k_count,
      k_pass = k_pass,
      hp_lambda = hp_lambda,
      henderson_k = henderson_k,
      CR_maxiter = CR_maxiter,
      CR_tol = CR_tol,
      CR_Binit = CR_Binit,
      CR_S0init = CR_S0init,
      CR_S0_shrink_diag = CR_S0_shrink_diag,
      CR_check_every = CR_check_every,
      CR_use_loglik = CR_use_loglik,
      CR_store_path = CR_store_path,
      CR_verb = CR_verb,
      CR_term = CR_term,
      curves = curves,
      mats = mats
    )
  } else {
    baseline_state <- st_check_baseline(
      baseline_state = baseline_state,
      curves = curves,
      mats = mats,
      reference = reference
    )
    for (name in c("lambdas", "cutoffs")) {
      value <- get(name, inherits = FALSE)
      if (!is.null(value) && !isTRUE(all.equal(value, baseline_state$config[[name]],
                                              tolerance = 0, check.attributes = FALSE))) {
        stop("Supplied `", name, "` does not match the cached baseline. Build a new baseline to change it.",
             call. = FALSE)
      }
    }
    if (!is.null(yields)) {
      yields <- .cc_yields(yields, baseline_state$mats, baseline_state$curves)
      saved <- baseline_state$single_curve_state$yields
      same <- identical(.cc_index(yields[[1]][[1]]), .cc_index(saved[[1]][[1]])) &&
        all(vapply(baseline_state$curves, function(curve) isTRUE(all.equal(
          as.matrix(yields[[curve]][-1]), as.matrix(saved[[curve]][-1]),
          tolerance = 0, check.attributes = FALSE)), logical(1)))
      if (!same) stop("Supplied yields do not match the cached baseline. Build a new baseline.", call. = FALSE)
    }
  }

  args <- list(
    baseline_state = baseline_state,
    shock_curves = shock_curves,
    shock_factors = shock_factors,
    shock_magnitude = shock_magnitude,
    shock_type = shock_type,
    shock_window = shock_window,
    shock_window_position = shock_window_position,
    reuse_BS0 = reuse_BS0,
    reuse_ECM = reuse_ECM,
    ECM_fallback = ECM_fallback
  )
  # Omitted controls inherit the cache, including when its settings differ
  # from this convenience wrapper's defaults. Fresh fits inherit their new cache.
  args <- c(args, mget(supplied_controls, envir = environment(), inherits = FALSE))
  do.call(mc_fit_shock_cached, args)
}

#############################################
#             DNS Profile Fitting           #
#############################################

sc_fit_shock <- function(yields = NULL, lambdas = NULL, cutoffs = NULL,
                         curves = NULL, mats = NULL, baseline_state = NULL,
                         shock_curves = NULL, shock_factors = "L",
                         shock_profile = NULL,
                         shock_type = c("additive", "multiplicative")) {
  shock_type <- match.arg(shock_type)
  if (is.null(baseline_state)) {
    baseline_state <- st_sc_state(yields, lambdas, cutoffs, curves, mats)
  }
  shocked_factors <- baseline_state$ns_factor
  if (!is.null(shock_profile)) {
    magnitudes <- if (is.data.frame(shock_profile)) shock_profile$shock_magnitude else shock_profile
    shocked_factors <- st_factor_profile(
      curves = baseline_state$curves, curve_dns_factor = baseline_state$ns_factor,
      shock_curves = shock_curves, shock_factors = shock_factors,
      shock_magnitudes = magnitudes, shock_type = shock_type
    )$shocked_factors
  }
  shocked_state <- st_shock_state(baseline_state, shocked_factors)
  list(
    yields = shocked_state$yields,
    ns_factor = shocked_state$ns_factor,
    shock_profile = shock_profile,
    shocked_factors = shocked_factors,
    baseline_state = baseline_state
  )
}

#############################################
#           Compatibility Aliases           #
#############################################

mc_fit_shock_baseline <- mc_fit_shock_base
mc_fit_shock_from_baseline <- mc_fit_shock_cached
mc_fit_shock_refit_from_baseline <- mc_fit_shock_refit
