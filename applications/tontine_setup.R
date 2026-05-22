#############################################
#        General Tontine Setup Layer        #
#############################################

tontine_validate_positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0) {
    stop(sprintf("`%s` must be a positive numeric scalar.", name))
  }
  as.numeric(x)
}

tontine_validate_nonnegative_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 0) {
    stop(sprintf("`%s` must be a non-negative numeric scalar.", name))
  }
  as.numeric(x)
}

tontine_validate_probability_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 0 || x > 1) {
    stop(sprintf("`%s` must be a numeric scalar between 0 and 1.", name))
  }
  as.numeric(x)
}

tontine_validate_probability_schedule <- function(x, n, name) {
  if (!is.numeric(x) || length(x) != n || anyNA(x) || any(x < 0 | x > 1)) {
    stop(sprintf("`%s` must be a numeric vector of length %s with values between 0 and 1.", name, n))
  }
  as.numeric(x)
}

tontine_create_member_pool <- function(N_members, initial_position = 0) {
  if (!is.numeric(N_members) || length(N_members) != 1 || is.na(N_members) ||
      N_members <= 0 || N_members != as.integer(N_members)) {
    stop("`N_members` must be a positive integer scalar.")
  }

  initial_position <- tontine_validate_nonnegative_scalar(initial_position, "initial_position")

  data.frame(
    member_id = seq_len(as.integer(N_members)),
    alive = TRUE,
    current_position = rep(initial_position, as.integer(N_members)),
    cumulative_contributions = 0,
    cumulative_withdrawals = 0,
    cumulative_withdrawal_penalties = 0,
    cumulative_payouts = 0,
    death_rung_id = NA_integer_,
    stringsAsFactors = FALSE
  )
}

tontine_create_rung_schedule <- function(purchase_years, bond_maturity_years = 1,
                                         contribution_per_member,
                                         death_rates = NULL,
                                         withdrawal_rates = NULL) {
  if (!is.numeric(purchase_years) || length(purchase_years) == 0 ||
      anyNA(purchase_years) || any(purchase_years != as.integer(purchase_years))) {
    stop("`purchase_years` must be a non-empty integer year vector.")
  }

  purchase_years <- as.integer(purchase_years)
  bond_maturity_years <- tontine_validate_positive_scalar(bond_maturity_years, "bond_maturity_years")
  if (bond_maturity_years != as.integer(bond_maturity_years)) {
    stop("`bond_maturity_years` must be an integer scalar.")
  }

  n <- length(purchase_years)
  contribution_per_member <- if (length(contribution_per_member) == 1) {
    rep(tontine_validate_nonnegative_scalar(contribution_per_member, "contribution_per_member"), n)
  } else {
    if (!is.numeric(contribution_per_member) || length(contribution_per_member) != n ||
        anyNA(contribution_per_member) || any(contribution_per_member < 0)) {
      stop("`contribution_per_member` must be a scalar or non-negative vector aligned with `purchase_years`.")
    }
    as.numeric(contribution_per_member)
  }

  if (is.null(death_rates)) {
    death_rates <- rep(NA_real_, n)
  } else {
    death_rates <- tontine_validate_probability_schedule(death_rates, n, "death_rates")
  }

  if (is.null(withdrawal_rates)) {
    withdrawal_rates <- rep(NA_real_, n)
  } else {
    withdrawal_rates <- tontine_validate_probability_schedule(withdrawal_rates, n, "withdrawal_rates")
  }

  data.frame(
    rung_id = seq_len(n),
    purchase_year = purchase_years,
    maturity_year = purchase_years + as.integer(bond_maturity_years),
    bond_maturity_years = as.integer(bond_maturity_years),
    contribution_per_member = contribution_per_member,
    death_rate = death_rates,
    withdrawal_rate = withdrawal_rates,
    stringsAsFactors = FALSE
  )
}

tontine_create_event_schedule <- function(purchase_years, death_rates,
                                          integrated_eta = NULL,
                                          withdrawal_intensity = NULL,
                                          integrated_lambda = NULL) {
  if (!is.numeric(purchase_years) || length(purchase_years) == 0 ||
      anyNA(purchase_years) || any(purchase_years != as.integer(purchase_years))) {
    stop("`purchase_years` must be a non-empty integer year vector.")
  }

  purchase_years <- as.integer(purchase_years)
  n <- length(purchase_years)
  death_rates <- tontine_validate_probability_schedule(death_rates, n, "death_rates")

  if (is.null(integrated_eta)) {
    integrated_eta <- integrated_lambda
  }

  if (is.null(integrated_eta)) {
    if (is.null(withdrawal_intensity)) {
      stop("Supply `integrated_eta` for the period withdrawal intensity.")
    }
    integrated_eta <- withdrawal_intensity
  }

  if (!is.numeric(integrated_eta) || length(integrated_eta) != n ||
      anyNA(integrated_eta) || any(integrated_eta < 0)) {
    stop("`integrated_eta` must be a non-negative numeric vector aligned with `purchase_years`.")
  }

  data.frame(
    rung_id = seq_len(n),
    purchase_year = purchase_years,
    death_rate = death_rates,
    integrated_eta = as.numeric(integrated_eta),
    integrated_lambda = as.numeric(integrated_eta),
    stringsAsFactors = FALSE
  )
}

tontine_simulate_cir_integrated_intensity <- function(n_rungs, eta0 = NULL, a = NULL,
                                                     b = NULL, sigma, dt = 1,
                                                     substeps = 12,
                                                     seed = NULL,
                                                     warn_feller = TRUE,
                                                     lambda0 = NULL,
                                                     kappa = NULL,
                                                     theta = NULL) {
  if (!is.numeric(n_rungs) || length(n_rungs) != 1 || is.na(n_rungs) ||
      n_rungs <= 0 || n_rungs != as.integer(n_rungs)) {
    stop("`n_rungs` must be a positive integer scalar.")
  }

  if (is.null(eta0)) {
    eta0 <- lambda0
  }
  if (is.null(a)) {
    a <- kappa
  }
  if (is.null(b)) {
    b <- theta
  }

  eta0 <- tontine_validate_nonnegative_scalar(eta0, "eta0")
  a <- tontine_validate_positive_scalar(a, "a")
  b <- tontine_validate_nonnegative_scalar(b, "b")
  sigma <- tontine_validate_nonnegative_scalar(sigma, "sigma")
  dt <- tontine_validate_positive_scalar(dt, "dt")

  if (!is.numeric(substeps) || length(substeps) != 1 || is.na(substeps) ||
      substeps <= 0 || substeps != as.integer(substeps)) {
    stop("`substeps` must be a positive integer scalar.")
  }

  feller_satisfied <- 2 * a * b >= sigma^2
  if (isTRUE(warn_feller) && !feller_satisfied) {
    warning(
      "CIR Feller condition 2 * a * b >= sigma^2 is not satisfied.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_rungs <- as.integer(n_rungs)
  substeps <- as.integer(substeps)
  h <- dt / substeps
  eta_current <- eta0
  out <- vector("list", n_rungs)

  for (rung_idx in seq_len(n_rungs)) {
    eta_start <- eta_current
    integrated_eta <- 0

    for (step_idx in seq_len(substeps)) {
      eta_before <- abs(eta_current)
      shock <- stats::rnorm(1)
      eta_after <- eta_before +
        a * (b - eta_before) * h +
        sigma * sqrt(eta_before) * sqrt(h) * shock
      eta_after <- abs(eta_after)
      integrated_eta <- integrated_eta + 0.5 * (eta_before + eta_after) * h
      eta_current <- eta_after
    }

    out[[rung_idx]] <- data.frame(
      rung_id = rung_idx,
      eta_start = eta_start,
      eta_end = eta_current,
      integrated_eta = integrated_eta,
      lambda_start = eta_start,
      lambda_end = eta_current,
      integrated_lambda = integrated_eta,
      feller_satisfied = feller_satisfied,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

tontine_withdrawal_settings <- function(withdrawal_fraction = 0.05,
                                        withdrawal_penalty_rate = 0.005) {
  withdrawal_fraction <- tontine_validate_probability_scalar(
    withdrawal_fraction,
    "withdrawal_fraction"
  )
  withdrawal_penalty_rate <- tontine_validate_probability_scalar(
    withdrawal_penalty_rate,
    "withdrawal_penalty_rate"
  )

  if (withdrawal_fraction + withdrawal_penalty_rate > 1) {
    stop("`withdrawal_fraction + withdrawal_penalty_rate` cannot exceed 1.")
  }

  list(
    withdrawal_fraction = withdrawal_fraction,
    withdrawal_penalty_rate = withdrawal_penalty_rate
  )
}

tontine_payout_settings <- function(allocation = c("position_weighted", "equal_survivor"),
                                    penalty_timing = c("next_rung", "current_rung"),
                                    payout_reduces_position = FALSE,
                                    position_update = c("unchanged", "reduce_by_payout", "replace_with_payout")) {
  allocation <- match.arg(allocation)
  penalty_timing <- match.arg(penalty_timing)
  position_update <- match.arg(position_update)

  if (!is.logical(payout_reduces_position) || length(payout_reduces_position) != 1 ||
      is.na(payout_reduces_position)) {
    stop("`payout_reduces_position` must be TRUE or FALSE.")
  }

  if (isTRUE(payout_reduces_position)) {
    position_update <- "reduce_by_payout"
  }

  list(
    allocation = allocation,
    penalty_timing = penalty_timing,
    payout_reduces_position = identical(position_update, "reduce_by_payout"),
    position_update = position_update
  )
}

tontine_initialize_state <- function(member_pool,
                                     withdrawal_settings = tontine_withdrawal_settings(),
                                     payout_settings = tontine_payout_settings()) {
  required_cols <- c(
    "member_id", "alive", "current_position", "cumulative_contributions",
    "cumulative_withdrawals", "cumulative_withdrawal_penalties",
    "cumulative_payouts", "death_rung_id"
  )

  if (!is.data.frame(member_pool) || !all(required_cols %in% colnames(member_pool))) {
    stop("`member_pool` must be created by `tontine_create_member_pool()`.")
  }

  list(
    members = member_pool,
    withdrawal_settings = withdrawal_settings,
    payout_settings = payout_settings,
    penalty_carryforward = 0,
    rung_history = data.frame(),
    member_events = data.frame()
  )
}

tontine_alive_member_ids <- function(state) {
  state$members$member_id[state$members$alive]
}

tontine_validate_member_ids <- function(member_ids, eligible_ids, name) {
  member_ids <- as.integer(member_ids)
  if (length(member_ids) == 0) {
    return(integer())
  }
  if (anyNA(member_ids) || any(!member_ids %in% eligible_ids)) {
    stop(sprintf("`%s` contains members that are not eligible for this event.", name))
  }
  unique(member_ids)
}

tontine_sample_death_ids <- function(state, death_rate) {
  death_rate <- tontine_validate_probability_scalar(death_rate, "death_rate")
  eligible_ids <- tontine_alive_member_ids(state)

  if (length(eligible_ids) == 0 || death_rate == 0) {
    return(integer())
  }

  death_indicator <- stats::runif(length(eligible_ids)) < death_rate
  as.integer(eligible_ids[death_indicator])
}

tontine_sample_withdrawal_ids <- function(state, withdrawal_intensity,
                                          eligible_member_ids = NULL) {
  withdrawal_intensity <- tontine_validate_nonnegative_scalar(
    withdrawal_intensity,
    "withdrawal_intensity"
  )

  if (is.null(eligible_member_ids)) {
    eligible_member_ids <- tontine_alive_member_ids(state)
  } else {
    eligible_member_ids <- tontine_validate_member_ids(
      eligible_member_ids,
      tontine_alive_member_ids(state),
      "eligible_member_ids"
    )
  }

  if (length(eligible_member_ids) == 0 || withdrawal_intensity == 0) {
    return(integer())
  }

  withdrawal_count <- stats::rpois(1, lambda = withdrawal_intensity)
  withdrawal_count <- min(withdrawal_count, length(eligible_member_ids))

  if (withdrawal_count == 0) {
    return(integer())
  }

  as.integer(sample(eligible_member_ids, size = withdrawal_count, replace = FALSE))
}

tontine_sample_withdrawal_ids_from_integrated_intensity <- function(state,
                                                                    integrated_lambda,
                                                                    eligible_member_ids = NULL,
                                                                    warn_on_cap = TRUE) {
  integrated_lambda <- tontine_validate_nonnegative_scalar(
    integrated_lambda,
    "integrated_lambda"
  )

  if (is.null(eligible_member_ids)) {
    eligible_member_ids <- tontine_alive_member_ids(state)
  } else {
    eligible_member_ids <- tontine_validate_member_ids(
      eligible_member_ids,
      tontine_alive_member_ids(state),
      "eligible_member_ids"
    )
  }

  survivor_count <- length(eligible_member_ids)
  poisson_mean <- survivor_count * integrated_lambda
  raw_withdrawal_count <- if (survivor_count == 0 || poisson_mean == 0) {
    0L
  } else {
    stats::rpois(1, lambda = poisson_mean)
  }
  withdrawal_count_capped <- raw_withdrawal_count > survivor_count
  withdrawal_count <- min(raw_withdrawal_count, survivor_count)

  if (withdrawal_count_capped && isTRUE(warn_on_cap)) {
    warning(
      "Poisson withdrawal count exceeded surviving member count; capped withdrawals at surviving count.",
      call. = FALSE
    )
  }

  withdrawal_member_ids <- if (withdrawal_count == 0) {
    integer()
  } else {
    as.integer(sample(eligible_member_ids, size = withdrawal_count, replace = FALSE))
  }

  list(
    withdrawal_member_ids = withdrawal_member_ids,
    integrated_lambda = integrated_lambda,
    poisson_mean = poisson_mean,
    raw_withdrawal_count = raw_withdrawal_count,
    withdrawal_count = withdrawal_count,
    withdrawal_count_capped = withdrawal_count_capped
  )
}

tontine_sample_rung_event_ids <- function(state, death_rate,
                                          integrated_lambda = NULL,
                                          withdrawal_intensity = NULL,
                                          warn_on_cap = TRUE) {
  death_member_ids <- tontine_sample_death_ids(
    state = state,
    death_rate = death_rate
  )

  survivors_after_death <- setdiff(tontine_alive_member_ids(state), death_member_ids)
  if (is.null(integrated_lambda)) {
    if (is.null(withdrawal_intensity)) {
      stop("Supply `integrated_lambda` for the period withdrawal intensity.")
    }
    integrated_lambda <- withdrawal_intensity
  }

  withdrawal_result <- tontine_sample_withdrawal_ids_from_integrated_intensity(
    state = state,
    integrated_lambda = integrated_lambda,
    eligible_member_ids = survivors_after_death,
    warn_on_cap = warn_on_cap
  )

  list(
    death_member_ids = death_member_ids,
    withdrawal_member_ids = withdrawal_result$withdrawal_member_ids,
    death_count = length(death_member_ids),
    survivors_after_death = length(survivors_after_death),
    integrated_lambda = withdrawal_result$integrated_lambda,
    withdrawal_poisson_mean = withdrawal_result$poisson_mean,
    raw_withdrawal_count = withdrawal_result$raw_withdrawal_count,
    withdrawal_count = withdrawal_result$withdrawal_count,
    withdrawal_count_capped = withdrawal_result$withdrawal_count_capped
  )
}

tontine_apply_contributions <- function(state, rung_id, contribution_per_member) {
  contribution_per_member <- tontine_validate_nonnegative_scalar(
    contribution_per_member,
    "contribution_per_member"
  )
  alive_idx <- which(state$members$alive)
  total_contribution <- contribution_per_member * length(alive_idx)

  state$members$current_position[alive_idx] <-
    state$members$current_position[alive_idx] + contribution_per_member
  state$members$cumulative_contributions[alive_idx] <-
    state$members$cumulative_contributions[alive_idx] + contribution_per_member

  if (contribution_per_member == 0 || length(alive_idx) == 0) {
    return(list(
      state = state,
      event = data.frame(),
      total_contribution = total_contribution
    ))
  }

  list(
    state = state,
    event = data.frame(
      rung_id = rung_id,
      event_type = "contribution",
      member_id = state$members$member_id[alive_idx],
      amount = contribution_per_member,
      penalty = 0,
      position_before = state$members$current_position[alive_idx] - contribution_per_member,
      position_after = state$members$current_position[alive_idx],
      stringsAsFactors = FALSE
    ),
    total_contribution = total_contribution
  )
}

tontine_apply_deaths <- function(state, rung_id, death_member_ids = integer()) {
  eligible_ids <- tontine_alive_member_ids(state)
  death_member_ids <- tontine_validate_member_ids(death_member_ids, eligible_ids, "death_member_ids")
  if (length(death_member_ids) == 0) {
    return(list(state = state, event = data.frame(), forfeited_position = 0))
  }

  death_idx <- match(death_member_ids, state$members$member_id)
  position_before <- state$members$current_position[death_idx]
  forfeited_position <- sum(position_before)

  state$members$alive[death_idx] <- FALSE
  state$members$current_position[death_idx] <- 0
  state$members$death_rung_id[death_idx] <- as.integer(rung_id)

  list(
    state = state,
    event = data.frame(
      rung_id = rung_id,
      event_type = "death",
      member_id = death_member_ids,
      amount = 0,
      penalty = 0,
      position_before = position_before,
      position_after = 0,
      stringsAsFactors = FALSE
    ),
    forfeited_position = forfeited_position
  )
}

tontine_apply_withdrawals <- function(state, rung_id, withdrawal_member_ids = integer()) {
  eligible_ids <- tontine_alive_member_ids(state)
  withdrawal_member_ids <- tontine_validate_member_ids(
    withdrawal_member_ids,
    eligible_ids,
    "withdrawal_member_ids"
  )
  if (length(withdrawal_member_ids) == 0) {
    return(list(
      state = state,
      event = data.frame(),
      total_withdrawal = 0,
      total_penalty = 0
    ))
  }

  withdrawal_idx <- match(withdrawal_member_ids, state$members$member_id)
  position_before <- state$members$current_position[withdrawal_idx]
  withdrawal_amount <-
    state$withdrawal_settings$withdrawal_fraction * position_before
  withdrawal_penalty <-
    state$withdrawal_settings$withdrawal_penalty_rate * position_before
  position_after <- position_before - withdrawal_amount - withdrawal_penalty

  state$members$current_position[withdrawal_idx] <- position_after
  state$members$cumulative_withdrawals[withdrawal_idx] <-
    state$members$cumulative_withdrawals[withdrawal_idx] + withdrawal_amount
  state$members$cumulative_withdrawal_penalties[withdrawal_idx] <-
    state$members$cumulative_withdrawal_penalties[withdrawal_idx] + withdrawal_penalty

  list(
    state = state,
    event = data.frame(
      rung_id = rung_id,
      event_type = "withdrawal",
      member_id = withdrawal_member_ids,
      amount = withdrawal_amount,
      penalty = withdrawal_penalty,
      position_before = position_before,
      position_after = position_after,
      stringsAsFactors = FALSE
    ),
    total_withdrawal = sum(withdrawal_amount),
    total_penalty = sum(withdrawal_penalty)
  )
}

tontine_allocate_rung_payout <- function(state, rung_id, maturity_proceeds,
                                         current_penalty_pool = 0) {
  maturity_proceeds <- tontine_validate_nonnegative_scalar(maturity_proceeds, "maturity_proceeds")
  current_penalty_pool <- tontine_validate_nonnegative_scalar(current_penalty_pool, "current_penalty_pool")
  payout_pool <- maturity_proceeds + current_penalty_pool
  alive_idx <- which(state$members$alive)

  if (length(alive_idx) == 0 || payout_pool == 0) {
    return(list(state = state, event = data.frame(), payout_pool = payout_pool))
  }

  if (state$payout_settings$allocation == "equal_survivor") {
    payout_weights <- rep(1 / length(alive_idx), length(alive_idx))
  } else {
    alive_position <- state$members$current_position[alive_idx]
    total_alive_position <- sum(alive_position)
    if (total_alive_position <= 0) {
      payout_weights <- rep(1 / length(alive_idx), length(alive_idx))
    } else {
      payout_weights <- alive_position / total_alive_position
    }
  }

  payout_amount <- payout_pool * payout_weights
  position_before <- state$members$current_position[alive_idx]

  state$members$cumulative_payouts[alive_idx] <-
    state$members$cumulative_payouts[alive_idx] + payout_amount

  if (identical(state$payout_settings$position_update, "reduce_by_payout")) {
    state$members$current_position[alive_idx] <-
      pmax(state$members$current_position[alive_idx] - payout_amount, 0)
  } else if (identical(state$payout_settings$position_update, "replace_with_payout")) {
    state$members$current_position[alive_idx] <- payout_amount
  }

  list(
    state = state,
    event = data.frame(
      rung_id = rung_id,
      event_type = "payout",
      member_id = state$members$member_id[alive_idx],
      amount = payout_amount,
      penalty = 0,
      position_before = position_before,
      position_after = state$members$current_position[alive_idx],
      stringsAsFactors = FALSE
    ),
    payout_pool = payout_pool
  )
}

tontine_process_rung <- function(state, rung,
                                 death_member_ids = integer(),
                                 withdrawal_member_ids = integer(),
                                 maturity_proceeds,
                                 withdrawal_diagnostics = NULL) {
  if (!is.data.frame(rung) || nrow(rung) != 1) {
    stop("`rung` must be a one-row data frame from `tontine_create_rung_schedule()`.")
  }

  rung_id <- as.integer(rung$rung_id[[1]])
  active_start_count <- sum(state$members$alive)
  active_position_start <- sum(state$members$current_position[state$members$alive])
  penalty_available_at_start <- state$penalty_carryforward

  contribution_result <- tontine_apply_contributions(
    state = state,
    rung_id = rung_id,
    contribution_per_member = rung$contribution_per_member[[1]]
  )
  state <- contribution_result$state

  death_result <- tontine_apply_deaths(
    state = state,
    rung_id = rung_id,
    death_member_ids = death_member_ids
  )
  state <- death_result$state

  withdrawal_result <- tontine_apply_withdrawals(
    state = state,
    rung_id = rung_id,
    withdrawal_member_ids = withdrawal_member_ids
  )
  state <- withdrawal_result$state

  current_penalty_pool <- if (state$payout_settings$penalty_timing == "current_rung") {
    penalty_available_at_start + withdrawal_result$total_penalty
  } else {
    penalty_available_at_start
  }

  payout_result <- tontine_allocate_rung_payout(
    state = state,
    rung_id = rung_id,
    maturity_proceeds = maturity_proceeds,
    current_penalty_pool = current_penalty_pool
  )
  state <- payout_result$state

  state$penalty_carryforward <- if (state$payout_settings$penalty_timing == "current_rung") {
    0
  } else {
    withdrawal_result$total_penalty
  }

  active_end_count <- sum(state$members$alive)
  active_position_end <- sum(state$members$current_position[state$members$alive])
  if (is.null(withdrawal_diagnostics)) {
    withdrawal_diagnostics <- list(
      integrated_lambda = NA_real_,
      withdrawal_poisson_mean = NA_real_,
      raw_withdrawal_count = length(withdrawal_member_ids),
      withdrawal_count_capped = FALSE
    )
  }

  rung_summary <- data.frame(
    rung_id = rung_id,
    purchase_year = rung$purchase_year[[1]],
    maturity_year = rung$maturity_year[[1]],
    active_start_count = active_start_count,
    active_position_start = active_position_start,
    contribution_per_member = rung$contribution_per_member[[1]],
    total_contribution = contribution_result$total_contribution,
    death_count = length(death_member_ids),
    death_forfeited_position = death_result$forfeited_position,
    active_after_deaths_count = active_start_count - length(death_member_ids),
    integrated_lambda = withdrawal_diagnostics$integrated_lambda,
    integrated_eta = withdrawal_diagnostics$integrated_lambda,
    withdrawal_poisson_mean = withdrawal_diagnostics$withdrawal_poisson_mean,
    raw_withdrawal_count = withdrawal_diagnostics$raw_withdrawal_count,
    withdrawal_count = length(withdrawal_member_ids),
    withdrawal_count_capped = withdrawal_diagnostics$withdrawal_count_capped,
    total_withdrawal = withdrawal_result$total_withdrawal,
    withdrawal_penalty = withdrawal_result$total_penalty,
    penalty_available_at_start = penalty_available_at_start,
    penalty_paid_with_rung = current_penalty_pool,
    maturity_proceeds = maturity_proceeds,
    payout_pool = payout_result$payout_pool,
    active_end_count = active_end_count,
    active_position_end = active_position_end,
    penalty_carryforward = state$penalty_carryforward,
    stringsAsFactors = FALSE
  )

  event_frames <- list(
    contribution_result$event,
    death_result$event,
    withdrawal_result$event,
    payout_result$event
  )
  event_frames <- event_frames[vapply(event_frames, nrow, integer(1)) > 0]
  rung_events <- if (length(event_frames) == 0) {
    data.frame()
  } else {
    do.call(rbind, event_frames)
  }

  state$rung_history <- rbind(state$rung_history, rung_summary)
  state$member_events <- rbind(state$member_events, rung_events)
  state
}

tontine_process_stochastic_rung <- function(state, rung, death_rate,
                                            integrated_lambda,
                                            maturity_proceeds,
                                            warn_on_cap = TRUE) {
  event_ids <- tontine_sample_rung_event_ids(
    state = state,
    death_rate = death_rate,
    integrated_lambda = integrated_lambda,
    warn_on_cap = warn_on_cap
  )

  tontine_process_rung(
    state = state,
    rung = rung,
    death_member_ids = event_ids$death_member_ids,
    withdrawal_member_ids = event_ids$withdrawal_member_ids,
    maturity_proceeds = maturity_proceeds,
    withdrawal_diagnostics = event_ids
  )
}

tontine_simulate_event_path <- function(member_pool, event_schedule, seed = NULL,
                                        warn_on_cap = TRUE) {
  required_cols <- c("rung_id", "purchase_year", "death_rate", "integrated_lambda")
  if (!is.data.frame(event_schedule) || !all(required_cols %in% colnames(event_schedule))) {
    stop("`event_schedule` must be created by `tontine_create_event_schedule()`.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  state <- tontine_initialize_state(member_pool)
  event_rows <- vector("list", nrow(event_schedule))

  for (i in seq_len(nrow(event_schedule))) {
    event_ids <- tontine_sample_rung_event_ids(
      state = state,
      death_rate = event_schedule$death_rate[[i]],
      integrated_lambda = event_schedule$integrated_lambda[[i]],
      warn_on_cap = warn_on_cap
    )

    death_result <- tontine_apply_deaths(
      state = state,
      rung_id = event_schedule$rung_id[[i]],
      death_member_ids = event_ids$death_member_ids
    )
    state <- death_result$state

    event_rows[[i]] <- data.frame(
      rung_id = event_schedule$rung_id[[i]],
      purchase_year = event_schedule$purchase_year[[i]],
      death_rate = event_schedule$death_rate[[i]],
      integrated_eta = event_ids$integrated_lambda,
      integrated_lambda = event_ids$integrated_lambda,
      withdrawal_poisson_mean = event_ids$withdrawal_poisson_mean,
      raw_withdrawal_count = event_ids$raw_withdrawal_count,
      withdrawal_count = event_ids$withdrawal_count,
      withdrawal_count_capped = event_ids$withdrawal_count_capped,
      death_count = event_ids$death_count,
      survivors_after_death = event_ids$survivors_after_death,
      death_member_ids = I(list(event_ids$death_member_ids)),
      withdrawal_member_ids = I(list(event_ids$withdrawal_member_ids)),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, event_rows)
}

tontine_validate_return_factors <- function(gross_return_factors) {
  required_cols <- c("rung_id", "gross_return_factor")
  if (!is.data.frame(gross_return_factors) || !all(required_cols %in% colnames(gross_return_factors))) {
    stop("`gross_return_factors` must contain `rung_id` and `gross_return_factor`.")
  }

  if (anyNA(gross_return_factors$gross_return_factor) ||
      any(gross_return_factors$gross_return_factor < 0)) {
    stop("`gross_return_factor` values must be non-negative and non-missing.")
  }

  gross_return_factors
}

tontine_replay_1yr_reinvestment_path <- function(member_pool, rung_schedule,
                                                 event_path,
                                                 gross_return_factors,
                                                 withdrawal_settings = tontine_withdrawal_settings(),
                                                 payout_settings = tontine_payout_settings(
                                                   allocation = "position_weighted",
                                                   penalty_timing = "next_rung",
                                                   position_update = "replace_with_payout"
                                                 ),
                                                 model = NA_character_,
                                                 scenario = NA_character_,
                                                 scenario_type = NA_character_,
                                                 shock_label = NA_character_) {
  if (!is.data.frame(rung_schedule) || !"rung_id" %in% colnames(rung_schedule)) {
    stop("`rung_schedule` must contain `rung_id`.")
  }

  required_event_cols <- c(
    "rung_id", "death_member_ids", "withdrawal_member_ids",
    "integrated_lambda", "withdrawal_poisson_mean",
    "raw_withdrawal_count", "withdrawal_count_capped"
  )
  if (!is.data.frame(event_path) || !all(required_event_cols %in% colnames(event_path))) {
    stop("`event_path` must be created by `tontine_simulate_event_path()`.")
  }

  gross_return_factors <- tontine_validate_return_factors(gross_return_factors)

  state <- tontine_initialize_state(
    member_pool = member_pool,
    withdrawal_settings = withdrawal_settings,
    payout_settings = payout_settings
  )
  replay_rows <- vector("list", nrow(rung_schedule))

  for (i in seq_len(nrow(rung_schedule))) {
    rung <- rung_schedule[i, , drop = FALSE]
    rung_id <- rung$rung_id[[1]]
    event_row <- event_path[event_path$rung_id == rung_id, , drop = FALSE]
    return_row <- gross_return_factors[gross_return_factors$rung_id == rung_id, , drop = FALSE]

    if (nrow(event_row) != 1) {
      stop("`event_path` must contain exactly one row per rung.")
    }
    if (nrow(return_row) != 1) {
      stop("`gross_return_factors` must contain exactly one row per rung.")
    }

    active_start_count <- sum(state$members$alive)
    active_position_start <- sum(state$members$current_position[state$members$alive])
    penalty_available_at_start <- state$penalty_carryforward

    contribution_result <- tontine_apply_contributions(
      state = state,
      rung_id = rung_id,
      contribution_per_member = rung$contribution_per_member[[1]]
    )
    state <- contribution_result$state

    death_member_ids <- event_row$death_member_ids[[1]]
    death_result <- tontine_apply_deaths(
      state = state,
      rung_id = rung_id,
      death_member_ids = death_member_ids
    )
    state <- death_result$state

    withdrawal_member_ids <- event_row$withdrawal_member_ids[[1]]
    withdrawal_result <- tontine_apply_withdrawals(
      state = state,
      rung_id = rung_id,
      withdrawal_member_ids = withdrawal_member_ids
    )
    state <- withdrawal_result$state

    survivor_position_after_withdrawals <- sum(state$members$current_position[state$members$alive])
    investable_position <- survivor_position_after_withdrawals + death_result$forfeited_position
    gross_return_factor <- return_row$gross_return_factor[[1]]
    maturity_proceeds <- investable_position * gross_return_factor

    current_penalty_pool <- if (state$payout_settings$penalty_timing == "current_rung") {
      penalty_available_at_start + withdrawal_result$total_penalty
    } else {
      penalty_available_at_start
    }

    payout_result <- tontine_allocate_rung_payout(
      state = state,
      rung_id = rung_id,
      maturity_proceeds = maturity_proceeds,
      current_penalty_pool = current_penalty_pool
    )
    state <- payout_result$state

    state$penalty_carryforward <- if (state$payout_settings$penalty_timing == "current_rung") {
      0
    } else {
      withdrawal_result$total_penalty
    }

    active_end_count <- sum(state$members$alive)
    active_position_end <- sum(state$members$current_position[state$members$alive])
    replay_rows[[i]] <- data.frame(
      model = model,
      scenario = scenario,
      scenario_type = scenario_type,
      shock_label = shock_label,
      rung_id = rung_id,
      purchase_year = rung$purchase_year[[1]],
      maturity_year = rung$maturity_year[[1]],
      active_start_count = active_start_count,
      active_position_start = active_position_start,
      contribution_per_member = rung$contribution_per_member[[1]],
      total_contribution = contribution_result$total_contribution,
      death_count = length(death_member_ids),
      death_forfeited_position = death_result$forfeited_position,
      active_after_deaths_count = active_start_count - length(death_member_ids),
      integrated_lambda = event_row$integrated_lambda[[1]],
      integrated_eta = event_row$integrated_lambda[[1]],
      withdrawal_poisson_mean = event_row$withdrawal_poisson_mean[[1]],
      raw_withdrawal_count = event_row$raw_withdrawal_count[[1]],
      withdrawal_count = length(withdrawal_member_ids),
      withdrawal_count_capped = event_row$withdrawal_count_capped[[1]],
      total_withdrawal = withdrawal_result$total_withdrawal,
      withdrawal_penalty = withdrawal_result$total_penalty,
      survivor_position_after_withdrawals = survivor_position_after_withdrawals,
      investable_position = investable_position,
      gross_return_factor = gross_return_factor,
      penalty_available_at_start = penalty_available_at_start,
      penalty_paid_with_rung = current_penalty_pool,
      maturity_proceeds = maturity_proceeds,
      payout_pool = payout_result$payout_pool,
      active_end_count = active_end_count,
      active_position_end = active_position_end,
      average_position_end = ifelse(active_end_count > 0, active_position_end / active_end_count, NA_real_),
      penalty_carryforward = state$penalty_carryforward,
      stringsAsFactors = FALSE
    )

    event_frames <- list(
      contribution_result$event,
      death_result$event,
      withdrawal_result$event,
      payout_result$event
    )
    event_frames <- event_frames[vapply(event_frames, nrow, integer(1)) > 0]
    rung_events <- if (length(event_frames) > 0) {
      do.call(rbind, event_frames)
    } else {
      data.frame()
    }

    if (nrow(rung_events) > 0) {
      rung_events$model <- model
      rung_events$scenario <- scenario
      rung_events$scenario_type <- scenario_type
      rung_events$shock_label <- shock_label
    }

    state$member_events <- rbind(state$member_events, rung_events)
  }

  rung_history <- do.call(rbind, replay_rows)
  rownames(rung_history) <- NULL
  rownames(state$member_events) <- NULL
  state$rung_history <- rung_history

  state
}
