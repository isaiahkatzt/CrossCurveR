#############################################
#          Tontine Application Helpers      #
#############################################

tontine_project_files <- function() {
  c(
    "core_formatting.R",
    "yield_estimation/curve_reformat.R",
    "yield_estimation/matrix_helpers.R",
    "yield_estimation/cr_core_estim.R",
    "yield_estimation/cr_eval.R",
    "yield_estimation/bl_single_curve.R",
    "yield_estimation/dly_mcm.R",
    "yield_estimation/bl_feature_extract.R",
    "yield_estimation/bl_normalize.R",
    "yield_estimation/bl_cross_estimation.R",
    "stress_testing/st_build_shock.R",
    "stress_testing/st_fit_shock.R"
  )
}

tontine_source_project <- function(root = here::here(), envir = .GlobalEnv,
                                   files = tontine_project_files()) {
  invisible(lapply(files, function(path) {
    sys.source(file.path(root, path), envir = envir)
  }))
}

tontine_mat_names <- function(mats) {
  paste0("X", unname(vapply(mats, numeric_to_matname, character(1))))
}

tontine_read_yields <- function(curve_files, data_dir = here::here("reconstructed_data"),
                                start_date = NULL, end_date = NULL) {
  if (!is.list(curve_files) || is.null(names(curve_files))) {
    stop("`curve_files` must be a named list.")
  }

  yields <- lapply(curve_files, function(file_name) {
    utils::read.csv(file.path(data_dir, file_name), stringsAsFactors = FALSE)
  })

  setNames(lapply(yields, function(curve_df) {
    curve_df$time <- as.Date(curve_df$time)
    curve_df <- curve_df[order(curve_df$time), , drop = FALSE]

    if (!is.null(start_date)) {
      curve_df <- curve_df[curve_df$time >= as.Date(start_date), , drop = FALSE]
    }
    if (!is.null(end_date)) {
      curve_df <- curve_df[curve_df$time <= as.Date(end_date), , drop = FALSE]
    }

    curve_df
  }), names(curve_files))
}

tontine_format_shock_label <- function(shock_spec) {
  if (!is.null(shock_spec$label)) {
    return(shock_spec$label)
  }

  paste(
    paste(shock_spec$shock_curves, collapse = "+"),
    paste(shock_spec$shock_factors, collapse = "+"),
    scales::number(shock_spec$shock_magnitude, accuracy = 0.001, trim = TRUE),
    sep = " "
  )
}

tontine_build_shock_specs <- function(shock_templates, shock_structures,
                                      shared_controls = list()) {
  specs <- list()

  for (template_name in names(shock_templates)) {
    for (structure_name in names(shock_structures)) {
      spec_name <- paste(template_name, structure_name, sep = "__")
      specs[[spec_name]] <- utils::modifyList(
        utils::modifyList(shock_templates[[template_name]], shock_structures[[structure_name]]),
        shared_controls
      )
    }
  }

  specs
}

#############################################
#        Tontine Bond Ladder Instrument     #
#############################################

tontine_first_observed_date <- function(curve_df, year_value) {
  curve_time <- as.Date(curve_df$time)
  year_rows <- curve_time[format(curve_time, "%Y") == as.character(year_value)]

  if (length(year_rows) == 0) {
    stop(sprintf("No reference-curve observations found for year %s.", year_value))
  }

  min(year_rows)
}

tontine_last_observed_date <- function(curve_df, year_value) {
  curve_time <- as.Date(curve_df$time)
  year_rows <- curve_time[format(curve_time, "%Y") == as.character(year_value)]

  if (length(year_rows) == 0) {
    stop(sprintf("No reference-curve observations found for year %s.", year_value))
  }

  max(year_rows)
}

tontine_validate_probability_vector <- function(q, purchase_years) {
  if (!is.numeric(q) || anyNA(q) || length(q) != length(purchase_years)) {
    stop("`q` must be a numeric vector with one value per purchase year.")
  }

  if (any(q < 0 | q > 1)) {
    stop("`q` entries must be probabilities between 0 and 1.")
  }

  as.numeric(q)
}

tontine_validate_bond_maturity_years <- function(bond_maturity_years) {
  if (!is.numeric(bond_maturity_years) || length(bond_maturity_years) != 1 ||
      is.na(bond_maturity_years) || bond_maturity_years <= 0 ||
      bond_maturity_years != as.integer(bond_maturity_years)) {
    stop("`bond_maturity_years` must be a positive integer scalar.")
  }

  as.integer(bond_maturity_years)
}

tontine_bond_tenor_col <- function(bond_maturity_years) {
  paste0("X", numeric_to_matname(tontine_validate_bond_maturity_years(bond_maturity_years) * 12L))
}

tontine_build_valuation_windows <- function(yields, reference_curve, purchase_years,
                                            lookback_years = 2,
                                            min_observations = 1,
                                            date_position = c("first", "last"),
                                            valuation_date_position = NULL) {
  date_position <- match.arg(date_position)
  if (is.null(valuation_date_position)) {
    valuation_date_position <- date_position
  }
  valuation_date_position <- match.arg(valuation_date_position, c("first", "last"))

  if (!is.list(yields) || is.null(names(yields)) || !reference_curve %in% names(yields)) {
    stop("`yields` must be a named list containing `reference_curve`.")
  }

  if (!is.numeric(purchase_years) || length(purchase_years) == 0 ||
      anyNA(purchase_years) || any(purchase_years != as.integer(purchase_years))) {
    stop("`purchase_years` must be a non-empty integer year vector.")
  }

  if (!is.numeric(lookback_years) || length(lookback_years) != 1 ||
      is.na(lookback_years) || lookback_years <= 0 ||
      lookback_years != as.integer(lookback_years)) {
    stop("`lookback_years` must be a positive integer scalar.")
  }

  if (!is.numeric(min_observations) || length(min_observations) != 1 ||
      is.na(min_observations) || min_observations < 1) {
    stop("`min_observations` must be a positive numeric scalar.")
  }

  purchase_years <- as.integer(purchase_years)
  lookback_years <- as.integer(lookback_years)
  reference_panel <- yields[[reference_curve]]
  reference_time <- as.Date(reference_panel$time)
  resolve_purchase_date <- if (date_position == "first") {
    tontine_first_observed_date
  } else {
    tontine_last_observed_date
  }
  resolve_valuation_date <- if (valuation_date_position == "first") {
    tontine_first_observed_date
  } else {
    tontine_last_observed_date
  }

  windows <- do.call(rbind, lapply(seq_along(purchase_years), function(i) {
    purchase_date <- resolve_purchase_date(reference_panel, purchase_years[[i]])
    valuation_date <- resolve_valuation_date(reference_panel, purchase_years[[i]])
    maturity_date <- resolve_purchase_date(reference_panel, purchase_years[[i]] + 1L)
    estimation_start <- seq(
      from = as.Date(valuation_date),
      by = sprintf("-%d years", lookback_years),
      length.out = 2
    )[[2]]
    estimation_end <- as.Date(valuation_date)
    n_reference_observations <- sum(reference_time >= estimation_start & reference_time <= estimation_end)

    if (n_reference_observations < min_observations) {
      stop(sprintf(
        "Valuation window for purchase year %s has fewer than %s reference observations.",
        purchase_years[[i]], min_observations
      ))
    }

    data.frame(
      rung_id = i,
      purchase_year = purchase_years[[i]],
      purchase_date = as.Date(purchase_date),
      valuation_date = as.Date(valuation_date),
      maturity_date = as.Date(maturity_date),
      lookback_years = lookback_years,
      estimation_start = as.Date(estimation_start),
      estimation_end = as.Date(estimation_end),
      n_reference_observations = n_reference_observations,
      stringsAsFactors = FALSE
    )
  }))

  rownames(windows) <- NULL
  windows
}

tontine_slice_yields_for_window <- function(yields, valuation_window,
                                            curves = names(yields)) {
  if (nrow(valuation_window) != 1) {
    stop("`valuation_window` must contain exactly one row.")
  }

  estimation_start <- as.Date(valuation_window$estimation_start[[1]])
  estimation_end <- as.Date(valuation_window$estimation_end[[1]])

  setNames(lapply(curves, function(curve_name) {
    curve_df <- yields[[curve_name]]
    if (is.null(curve_df)) {
      stop(sprintf("`yields` is missing curve `%s`.", curve_name))
    }

    curve_time <- as.Date(curve_df$time)
    out <- curve_df[curve_time >= estimation_start & curve_time <= estimation_end, , drop = FALSE]
    out[order(as.Date(out$time)), , drop = FALSE]
  }), curves)
}

tontine_simulate_population <- function(N_members, q, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  alive <- rep(TRUE, N_members)
  expected_p_k <- cumprod(1 - q)
  population_path <- vector("list", length(q))

  for (i in seq_along(q)) {
    population_start <- sum(alive)
    alive_idx <- which(alive)
    deaths <- stats::runif(population_start) < q[[i]]
    alive[alive_idx[deaths]] <- FALSE
    population_end <- sum(alive)

    population_path[[i]] <- data.frame(
      rung_id = i,
      q_k = q[[i]],
      expected_p_k = expected_p_k[[i]],
      population_start = population_start,
      deaths = sum(deaths),
      population_end = population_end,
      realized_p_k = population_end / N_members,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, population_path)
}

tontine_build_bond_ladder_schedule <- function(yields, reference_curve, purchase_years,
                                               curves, total_target_par, q,
                                               population_simulation = NULL,
                                               target_rung_par_per_member = NULL,
                                               bond_maturity_years = 1) {
  if (!is.list(yields) || is.null(names(yields)) || !reference_curve %in% names(yields)) {
    stop("`yields` must be a named list containing `reference_curve`.")
  }

  if (!is.character(curves) || length(curves) == 0 || anyNA(curves) || any(!nzchar(curves))) {
    stop("`curves` must be a non-empty character vector.")
  }

  if (anyDuplicated(curves)) {
    stop("`curves` must contain unique curve names.")
  }

  if (!is.numeric(purchase_years) || length(purchase_years) == 0 ||
      anyNA(purchase_years) || any(purchase_years != as.integer(purchase_years))) {
    stop("`purchase_years` must be a non-empty integer year vector.")
  }

  if (!is.numeric(total_target_par) || length(total_target_par) != 1 ||
      is.na(total_target_par) || total_target_par <= 0) {
    stop("`total_target_par` must be a positive numeric scalar.")
  }

  purchase_years <- as.integer(purchase_years)
  bond_maturity_years <- tontine_validate_bond_maturity_years(bond_maturity_years)
  q <- tontine_validate_probability_vector(q, purchase_years)
  p_k <- cumprod(1 - q)

  if (!is.null(population_simulation)) {
    required_population_cols <- c(
      "rung_id", "q_k", "expected_p_k", "population_start",
      "deaths", "population_end", "realized_p_k"
    )

    if (!all(required_population_cols %in% colnames(population_simulation))) {
      stop("`population_simulation` must be created by `tontine_simulate_population()`.")
    }

    if (nrow(population_simulation) != length(purchase_years)) {
      stop("`population_simulation` must contain one row per purchase year.")
    }

    if (!all(population_simulation$rung_id == seq_along(purchase_years))) {
      stop("`population_simulation$rung_id` must align with `purchase_years`.")
    }

    if (any(abs(population_simulation$q_k - q) > sqrt(.Machine$double.eps))) {
      stop("`population_simulation$q_k` must match `q`.")
    }

    if (is.null(target_rung_par_per_member) ||
        !is.numeric(target_rung_par_per_member) ||
        length(target_rung_par_per_member) != 1 ||
        is.na(target_rung_par_per_member) ||
        target_rung_par_per_member <= 0) {
      stop("`target_rung_par_per_member` must be a positive scalar when `population_simulation` is supplied.")
    }
  }

  reference_panel <- yields[[reference_curve]]
  target_par <- if (is.null(population_simulation)) {
    rep(as.numeric(total_target_par), length(purchase_years))
  } else {
    population_simulation$population_start * as.numeric(target_rung_par_per_member)
  }
  next_rung_required_par <- if (is.null(population_simulation)) {
    rep(NA_real_, length(purchase_years))
  } else {
    population_simulation$population_end * as.numeric(target_rung_par_per_member)
  }
  mortality_credit_par <- if (is.null(population_simulation)) {
    rep(NA_real_, length(purchase_years))
  } else {
    population_simulation$deaths * as.numeric(target_rung_par_per_member)
  }
  per_survivor_maturity_value <- if (is.null(population_simulation)) {
    rep(NA_real_, length(purchase_years))
  } else {
    ifelse(
      population_simulation$population_end > 0,
      target_par / population_simulation$population_end,
      NA_real_
    )
  }
  mortality_credit_per_survivor <- if (is.null(population_simulation)) {
    rep(NA_real_, length(purchase_years))
  } else {
    per_survivor_maturity_value - as.numeric(target_rung_par_per_member)
  }

  rung_schedule <- data.frame(
    rung_id = seq_along(purchase_years),
    purchase_year = purchase_years,
    purchase_date = as.Date(vapply(purchase_years, function(year_value) {
      as.character(tontine_first_observed_date(reference_panel, year_value))
    }, character(1))),
    maturity_date = as.Date(vapply(purchase_years + bond_maturity_years, function(year_value) {
      as.character(tontine_first_observed_date(reference_panel, year_value))
    }, character(1))),
    bond_maturity_years = bond_maturity_years,
    target_par = target_par,
    q_k = q,
    p_k = p_k,
    target_rung_par_per_member = if (is.null(target_rung_par_per_member)) {
      NA_real_
    } else {
      as.numeric(target_rung_par_per_member)
    },
    population_start = if (is.null(population_simulation)) NA_integer_ else population_simulation$population_start,
    deaths = if (is.null(population_simulation)) NA_integer_ else population_simulation$deaths,
    population_end = if (is.null(population_simulation)) NA_integer_ else population_simulation$population_end,
    realized_p_k = if (is.null(population_simulation)) NA_real_ else population_simulation$realized_p_k,
    next_rung_required_par = next_rung_required_par,
    mortality_credit_par = mortality_credit_par,
    per_survivor_maturity_value = per_survivor_maturity_value,
    mortality_credit_per_survivor = mortality_credit_per_survivor,
    stringsAsFactors = FALSE
  )

  schedule <- do.call(rbind, lapply(seq_len(nrow(rung_schedule)), function(i) {
    data.frame(
      rung_id = rep(rung_schedule$rung_id[[i]], length(curves)),
      purchase_year = rep(rung_schedule$purchase_year[[i]], length(curves)),
      purchase_date = rep(rung_schedule$purchase_date[[i]], length(curves)),
      maturity_date = rep(rung_schedule$maturity_date[[i]], length(curves)),
      bond_maturity_years = rep(rung_schedule$bond_maturity_years[[i]], length(curves)),
      country = curves,
      target_par = rep(rung_schedule$target_par[[i]], length(curves)),
      country_par = rung_schedule$target_par[[i]] / length(curves),
      q_k = rep(rung_schedule$q_k[[i]], length(curves)),
      p_k = rep(rung_schedule$p_k[[i]], length(curves)),
      target_rung_par_per_member = rep(rung_schedule$target_rung_par_per_member[[i]], length(curves)),
      population_start = rep(rung_schedule$population_start[[i]], length(curves)),
      deaths = rep(rung_schedule$deaths[[i]], length(curves)),
      population_end = rep(rung_schedule$population_end[[i]], length(curves)),
      realized_p_k = rep(rung_schedule$realized_p_k[[i]], length(curves)),
      next_rung_required_par = rep(rung_schedule$next_rung_required_par[[i]], length(curves)),
      mortality_credit_par = rung_schedule$mortality_credit_par[[i]] / length(curves),
      per_survivor_maturity_value = rep(rung_schedule$per_survivor_maturity_value[[i]], length(curves)),
      mortality_credit_per_survivor = rep(rung_schedule$mortality_credit_per_survivor[[i]], length(curves)),
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(population_simulation)) {
    schedule$survival_adjusted_par <- schedule$country_par * schedule$p_k
    schedule$purchase_par <- schedule$survival_adjusted_par
  } else {
    schedule$survival_adjusted_par <- schedule$country_par
    schedule$purchase_par <- schedule$country_par
  }
  rownames(schedule) <- NULL
  schedule[, c(
    "rung_id", "purchase_year", "purchase_date", "maturity_date",
    "bond_maturity_years", "country", "target_par", "country_par", "q_k", "p_k",
    "target_rung_par_per_member", "population_start", "deaths",
    "population_end", "realized_p_k", "next_rung_required_par",
    "mortality_credit_par", "per_survivor_maturity_value",
    "mortality_credit_per_survivor", "survival_adjusted_par", "purchase_par"
  )]
}

tontine_extract_curve_rate <- function(curve_panel, target_date, tenor_col) {
  if (!tenor_col %in% colnames(curve_panel)) {
    stop(sprintf("Curve panel is missing required tenor column `%s`.", tenor_col))
  }

  panel_time <- as.Date(curve_panel$time)
  row_idx <- which(panel_time == as.Date(target_date))

  if (length(row_idx) != 1) {
    stop("Curve panel must contain exactly one row for each ladder purchase date.")
  }

  as.numeric(curve_panel[row_idx, tenor_col, drop = TRUE])
}

interpolate_curve_rate <- function(curve_panel, target_date, target_months,
                                   mats, mat_names, rate_scale = c("percent", "decimal")) {
  rate_scale <- match.arg(rate_scale)
  if (!is.data.frame(curve_panel) || !"time" %in% colnames(curve_panel)) {
    stop("`curve_panel` must be a data frame with a `time` column.")
  }
  if (!is.numeric(target_months) || length(target_months) != 1 ||
      is.na(target_months) || target_months <= 0) {
    stop("`target_months` must be a positive numeric scalar.")
  }
  if (!is.numeric(mats) || length(mats) == 0 || anyNA(mats) || any(mats <= 0)) {
    stop("`mats` must be a positive numeric maturity vector in months.")
  }
  if (!is.character(mat_names) || length(mat_names) != length(mats) ||
      anyNA(mat_names) || any(!nzchar(mat_names))) {
    stop("`mat_names` must be a character vector aligned with `mats`.")
  }
  if (!all(mat_names %in% colnames(curve_panel))) {
    stop("`curve_panel` is missing one or more maturity columns from `mat_names`.")
  }

  panel_time <- as.Date(curve_panel$time)
  row_idx <- which(panel_time == as.Date(target_date))
  if (length(row_idx) != 1) {
    stop("Curve panel must contain exactly one row for `target_date`.")
  }

  node_years <- as.numeric(mats) / 12
  target_years <- as.numeric(target_months) / 12
  zero_rates <- as.numeric(curve_panel[row_idx, mat_names, drop = TRUE]) / 100
  log_discount_nodes <- -zero_rates * node_years
  log_discount <- stats::approx(
    x = node_years,
    y = log_discount_nodes,
    xout = target_years,
    rule = 2,
    ties = "ordered"
  )$y

  interpolated_decimal_rate <- as.numeric(-log_discount / target_years)
  if (identical(rate_scale, "percent")) {
    return(100 * interpolated_decimal_rate)
  }

  interpolated_decimal_rate
}

tontine_mat_months_from_names <- function(mat_names) {
  suffix <- sub("^X", "", mat_names)
  out <- ifelse(
    grepl("M$", suffix),
    as.numeric(sub("M$", "", suffix)),
    as.numeric(sub("Y$", "", suffix)) * 12
  )

  if (anyNA(out)) {
    stop("`mat_names` must use maturity labels such as `X06M`, `X01Y`, or `X05Y`.")
  }

  out
}

tontine_interpolate_zero_rate <- function(curve_row, maturity_years, mats, mat_names) {
  if (!all(mat_names %in% colnames(curve_row))) {
    stop("Curve row is missing one or more maturity columns.")
  }

  maturity_years <- as.numeric(maturity_years)
  if (length(maturity_years) != 1 || is.na(maturity_years) || maturity_years <= 0) {
    stop("`maturity_years` must be a positive scalar.")
  }

  node_years <- as.numeric(mats) / 12
  zero_rates <- as.numeric(curve_row[1, mat_names, drop = TRUE]) / 100
  log_discount_nodes <- -zero_rates * node_years
  log_discount <- stats::approx(
    x = node_years,
    y = log_discount_nodes,
    xout = maturity_years,
    rule = 2,
    ties = "ordered"
  )$y

  as.numeric(-log_discount / maturity_years)
}

tontine_price_bond_ladder <- function(curve_panels, ladder_schedule, mat_names,
                                      model_name, scenario_name, scenario_type,
                                      shock_label = NA_character_) {
  if (!is.list(curve_panels) || is.null(names(curve_panels))) {
    stop("`curve_panels` must be a named list of curve data frames.")
  }

  required_cols <- c(
    "rung_id", "purchase_year", "purchase_date", "maturity_date",
    "country", "target_par", "country_par", "q_k", "p_k",
    "survival_adjusted_par"
  )

  if (!all(required_cols %in% colnames(ladder_schedule))) {
    stop("`ladder_schedule` must be created by `tontine_build_bond_ladder_schedule()`.")
  }

  missing_curves <- setdiff(unique(ladder_schedule$country), names(curve_panels))
  if (length(missing_curves) > 0) {
    stop("`curve_panels` is missing one or more ladder countries.")
  }

  out <- ladder_schedule
  out$model <- model_name
  out$scenario <- scenario_name
  out$scenario_type <- scenario_type
  out$shock_label <- shock_label
  if (!"bond_maturity_years" %in% colnames(out)) {
    out$bond_maturity_years <- 1L
  }
  out$tenor_col <- vapply(out$bond_maturity_years, tontine_bond_tenor_col, character(1))

  missing_tenors <- setdiff(unique(out$tenor_col), mat_names)
  if (length(missing_tenors) > 0) {
    stop(sprintf(
      "`mat_names` must include required ladder tenor column(s): %s.",
      paste(missing_tenors, collapse = ", ")
    ))
  }

  out$zero_rate <- vapply(seq_len(nrow(out)), function(i) {
    rate_pct <- tontine_extract_curve_rate(
      curve_panel = curve_panels[[out$country[[i]]]],
      target_date = out$purchase_date[[i]],
      tenor_col = out$tenor_col[[i]]
    )
    rate_pct / 100
  }, numeric(1))
  out$zero_rate_1y <- out$zero_rate
  out$priced_par <- if ("purchase_par" %in% colnames(out)) {
    out$purchase_par
  } else {
    out$survival_adjusted_par
  }
  out$discount_factor <- exp(-out$zero_rate * out$bond_maturity_years)
  out$purchase_cost <- out$priced_par * out$discount_factor
  out$maturity_proceeds <- out$priced_par

  optional_population_cols <- intersect(
    c(
      "target_rung_par_per_member", "population_start", "deaths",
      "population_end", "realized_p_k", "next_rung_required_par",
      "mortality_credit_par", "per_survivor_maturity_value",
      "mortality_credit_per_survivor", "purchase_par"
    ),
    colnames(out)
  )

  out[, c(
    "model", "scenario", "scenario_type", "shock_label",
    "rung_id", "purchase_year", "purchase_date", "maturity_date",
    "bond_maturity_years", "tenor_col", "country", "target_par", "country_par", "q_k", "p_k",
    optional_population_cols,
    "survival_adjusted_par", "priced_par", "zero_rate", "zero_rate_1y", "discount_factor",
    "purchase_cost", "maturity_proceeds"
  )]
}

tontine_mark_to_market_ladder <- function(country_holdings, curve_panels, mats,
                                          mat_names, valuation_date,
                                          scenario_name, scenario_type,
                                          shock_label = NA_character_) {
  valuation_date <- as.Date(valuation_date)
  required_cols <- c(
    "model", "rung_id", "purchase_year", "purchase_date", "maturity_date",
    "country", "priced_par", "purchase_cost"
  )

  if (!all(required_cols %in% colnames(country_holdings))) {
    stop("`country_holdings` is missing required mark-to-market columns.")
  }

  outstanding <- country_holdings[
    as.Date(country_holdings$purchase_date) < valuation_date &
      as.Date(country_holdings$maturity_date) > valuation_date,
    ,
    drop = FALSE
  ]

  if (nrow(outstanding) == 0) {
    return(data.frame())
  }

  outstanding$scenario <- scenario_name
  outstanding$scenario_type <- scenario_type
  outstanding$shock_label <- shock_label
  outstanding$valuation_date <- valuation_date
  outstanding$remaining_maturity_years <-
    as.numeric(as.Date(outstanding$maturity_date) - valuation_date) / 365.25
  outstanding$mtm_zero_rate <- vapply(seq_len(nrow(outstanding)), function(i) {
    curve_panel <- curve_panels[[outstanding$country[[i]]]]
    if (is.null(curve_panel)) {
      stop(sprintf("`curve_panels` is missing curve `%s`.", outstanding$country[[i]]))
    }
    row_idx <- which(as.Date(curve_panel$time) == valuation_date)
    if (length(row_idx) != 1) {
      stop("Curve panel must contain exactly one row for the mark-to-market valuation date.")
    }

    tontine_interpolate_zero_rate(
      curve_row = curve_panel[row_idx, , drop = FALSE],
      maturity_years = outstanding$remaining_maturity_years[[i]],
      mats = mats,
      mat_names = mat_names
    )
  }, numeric(1))
  outstanding$mtm_discount_factor <- exp(
    -outstanding$mtm_zero_rate * outstanding$remaining_maturity_years
  )
  outstanding$mtm_value <- outstanding$priced_par * outstanding$mtm_discount_factor
  outstanding$mtm_profit_loss_since_purchase <- outstanding$mtm_value - outstanding$purchase_cost

  optional_cols <- intersect(
    c("bond_maturity_years", "tenor_col", "target_par", "country_par"),
    colnames(outstanding)
  )

  outstanding[, c(
    "model", "scenario", "scenario_type", "shock_label",
    "valuation_date", "rung_id", "purchase_year", "purchase_date", "maturity_date",
    optional_cols, "country", "priced_par", "purchase_cost",
    "remaining_maturity_years", "mtm_zero_rate", "mtm_discount_factor",
    "mtm_value", "mtm_profit_loss_since_purchase"
  )]
}

tontine_compare_mtm_shock <- function(baseline_holdings, baseline_curve_models,
                                      shocked_curve_models, mats, mat_names,
                                      valuation_date, shock_label) {
  common_models <- intersect(names(shocked_curve_models), names(baseline_curve_models))

  do.call(rbind, lapply(common_models, function(model_name) {
    model_holdings <- baseline_holdings[
      baseline_holdings$model == baseline_curve_models[[model_name]]$model,
      ,
      drop = FALSE
    ]

    baseline_mtm <- tontine_mark_to_market_ladder(
      country_holdings = model_holdings,
      curve_panels = baseline_curve_models[[model_name]]$curve_panels,
      mats = mats,
      mat_names = mat_names,
      valuation_date = valuation_date,
      scenario_name = paste("Baseline MTM", baseline_curve_models[[model_name]]$model),
      scenario_type = "baseline_mtm",
      shock_label = NA_character_
    )
    shocked_mtm <- tontine_mark_to_market_ladder(
      country_holdings = model_holdings,
      curve_panels = shocked_curve_models[[model_name]]$curve_panels,
      mats = mats,
      mat_names = mat_names,
      valuation_date = valuation_date,
      scenario_name = paste(shock_label, shocked_curve_models[[model_name]]$model),
      scenario_type = "shock_mtm",
      shock_label = shock_label
    )

    key_cols <- c("model", "rung_id", "purchase_year", "country")
    if (nrow(baseline_mtm) == 0 || nrow(shocked_mtm) == 0) {
      return(data.frame(
        model = character(),
        rung_id = integer(),
        purchase_year = integer(),
        country = character(),
        mtm_zero_rate_baseline = numeric(),
        mtm_value_baseline = numeric(),
        mtm_zero_rate_shocked = numeric(),
        mtm_value_shocked = numeric(),
        shock_label = character(),
        valuation_date = as.Date(character()),
        mtm_value_delta = numeric(),
        mtm_value_delta_pct = numeric(),
        mtm_zero_rate_delta = numeric(),
        stringsAsFactors = FALSE
      ))
    }

    out <- merge(
      baseline_mtm[, c(key_cols, "mtm_zero_rate", "mtm_value"), drop = FALSE],
      shocked_mtm[, c(key_cols, "mtm_zero_rate", "mtm_value"), drop = FALSE],
      by = key_cols,
      suffixes = c("_baseline", "_shocked")
    )
    out$shock_label <- shock_label
    out$valuation_date <- valuation_date
    out$mtm_value_delta <- out$mtm_value_shocked - out$mtm_value_baseline
    out$mtm_value_delta_pct <- out$mtm_value_delta / out$mtm_value_baseline
    out$mtm_zero_rate_delta <- out$mtm_zero_rate_shocked - out$mtm_zero_rate_baseline
    out
  }))
}

tontine_resolve_initial_cash <- function(baseline_country_holdings) {
  if (!"rung_id" %in% colnames(baseline_country_holdings) ||
      !"purchase_cost" %in% colnames(baseline_country_holdings)) {
    stop("`baseline_country_holdings` must contain `rung_id` and `purchase_cost`.")
  }

  first_rung <- min(baseline_country_holdings$rung_id)
  sum(baseline_country_holdings$purchase_cost[baseline_country_holdings$rung_id == first_rung])
}

tontine_apply_ladder_cash_account <- function(country_holdings, initial_cash) {
  if (!is.numeric(initial_cash) || length(initial_cash) != 1 || is.na(initial_cash)) {
    stop("`initial_cash` must be a numeric scalar.")
  }

  required_cols <- c(
    "model", "scenario", "scenario_type", "shock_label",
    "rung_id", "purchase_year", "purchase_date", "maturity_date",
    "purchase_cost", "maturity_proceeds"
  )

  if (!all(required_cols %in% colnames(country_holdings))) {
    stop("`country_holdings` is missing required ladder cash-account columns.")
  }

  rung_keys <- unique(country_holdings[, c(
    "model", "scenario", "scenario_type", "shock_label",
    "rung_id", "purchase_year", "purchase_date", "maturity_date"
  )])
  rung_keys <- rung_keys[order(rung_keys$rung_id), , drop = FALSE]

  rung_purchase_cost <- stats::aggregate(
    purchase_cost ~ rung_id,
    data = country_holdings,
    FUN = sum
  )
  rung_maturity_proceeds <- stats::aggregate(
    maturity_proceeds ~ rung_id,
    data = country_holdings,
    FUN = sum
  )

  rung_account <- merge(rung_keys, rung_purchase_cost, by = "rung_id", all.x = TRUE)
  rung_account <- merge(rung_account, rung_maturity_proceeds, by = "rung_id", all.x = TRUE)
  if ("mortality_credit_par" %in% colnames(country_holdings) &&
      any(!is.na(country_holdings$mortality_credit_par))) {
    rung_mortality_credit <- stats::aggregate(
      mortality_credit_par ~ rung_id,
      data = country_holdings,
      FUN = sum
    )
    rung_account <- merge(rung_account, rung_mortality_credit, by = "rung_id", all.x = TRUE)
  } else {
    rung_account$mortality_credit_par <- NA_real_
  }
  rung_account <- rung_account[order(rung_account$rung_id), , drop = FALSE]

  n <- nrow(rung_account)
  cash_start <- numeric(n)
  maturity_proceeds_from_prior <- numeric(n)
  available_cash <- numeric(n)
  cash_end <- numeric(n)
  purchase_dates <- as.Date(rung_account$purchase_date)
  maturity_dates <- as.Date(rung_account$maturity_date)

  for (i in seq_len(n)) {
    cash_start[[i]] <- if (i == 1) initial_cash else cash_end[[i - 1]]
    maturity_proceeds_from_prior[[i]] <- if (i == 1) {
      0
    } else {
      sum(
        rung_account$maturity_proceeds[
          maturity_dates > purchase_dates[[i - 1]] &
            maturity_dates <= purchase_dates[[i]]
        ],
        na.rm = TRUE
      )
    }
    available_cash[[i]] <- cash_start[[i]] + maturity_proceeds_from_prior[[i]]
    cash_end[[i]] <- available_cash[[i]] - rung_account$purchase_cost[[i]]
  }

  rung_account$initial_cash <- initial_cash
  rung_account$cash_start <- cash_start
  rung_account$maturity_proceeds_from_prior <- maturity_proceeds_from_prior
  rung_account$available_cash <- available_cash
  rung_account$rung_profit_loss <- rung_account$maturity_proceeds - rung_account$purchase_cost
  rung_account$cumulative_profit_loss <- cumsum(rung_account$rung_profit_loss)
  rung_account$cash_end <- cash_end
  rung_account$reserve_draw <- pmax(-cash_end, 0)
  rung_account$cumulative_reserve_draw <- cummax(rung_account$reserve_draw)

  rung_account[, c(
    "model", "scenario", "scenario_type", "shock_label",
    "rung_id", "purchase_year", "purchase_date", "maturity_date",
    "initial_cash", "cash_start", "maturity_proceeds_from_prior",
    "available_cash", "purchase_cost", "maturity_proceeds", "mortality_credit_par",
    "rung_profit_loss", "cumulative_profit_loss", "cash_end",
    "reserve_draw", "cumulative_reserve_draw"
  )]
}

#############################################
#        Baseline Curve Model Validation    #
#############################################

tontine_curve_fit_to_panel <- function(time_index, fitted_matrix, mat_names) {
  fitted_matrix <- as.matrix(fitted_matrix)

  if (ncol(fitted_matrix) != length(mat_names)) {
    stop("`fitted_matrix` column count must match `mat_names`.")
  }

  out <- data.frame(
    time = as.Date(time_index),
    as.data.frame(fitted_matrix, check.names = FALSE),
    check.names = FALSE
  )
  colnames(out) <- c("time", mat_names)
  out
}

tontine_curve_yield_list_to_panels <- function(curve_yields, curves, mat_names) {
  setNames(lapply(curves, function(curve_name) {
    if (is.null(curve_yields[[curve_name]])) {
      stop(sprintf("Missing curve yield panel for `%s`.", curve_name))
    }

    curve_df <- curve_yields[[curve_name]]
    curve_df$time <- as.Date(curve_df$time)
    curve_df <- curve_df[order(curve_df$time), , drop = FALSE]

    if (ncol(curve_df) != length(mat_names) + 1) {
      stop("Curve yield panels must contain one time column plus one column per maturity.")
    }

    colnames(curve_df) <- c("time", mat_names)
    curve_df
  }), curves)
}

tontine_split_curve_panel <- function(panel_df, curves, mat_names) {
  tenor_suffix <- sub("^X", "", mat_names)

  setNames(lapply(curves, function(curve_name) {
    curve_cols <- paste0(curve_name, ".", tenor_suffix)

    if (!all(curve_cols %in% colnames(panel_df))) {
      stop(sprintf("Combined curve panel is missing one or more columns for `%s`.", curve_name))
    }

    curve_df <- data.frame(
      time = as.Date(panel_df$time),
      panel_df[, curve_cols, drop = FALSE],
      check.names = FALSE
    )
    colnames(curve_df) <- c("time", mat_names)
    curve_df
  }), curves)
}

tontine_observed_curve_panels <- function(yields, curves, mat_names) {
  tontine_curve_yield_list_to_panels(yields, curves, mat_names)
}

tontine_fit_dns_baseline_curves <- function(yields, curves, mats, mat_names,
                                            lambdas, cutoffs) {
  dns_fit <- st_prepare_single_curve_state(
    yields = yields,
    lambdas = lambdas,
    cutoffs = cutoffs,
    curves = curves,
    mats = mats
  )

  curve_panels <- setNames(lapply(curves, function(curve_name) {
    tontine_curve_fit_to_panel(
      time_index = dns_fit$time,
      fitted_matrix = dns_fit$bl_yields$nsfit[[curve_name]]$NSfit,
      mat_names = mat_names
    )
  }), curves)

  list(
    model = "DNS",
    fit = dns_fit,
    curve_panels = curve_panels
  )
}

tontine_fit_dly_baseline_curves <- function(yields, curves, mats, mat_names,
                                            lambdas, reference_curve,
                                            component = "global_plus_idiosyncratic") {
  dly_fit_obj <- dly_fit(
    yields = yields,
    mats = mats,
    lambdas = lambdas,
    curves = curves,
    reference = reference_curve
  )

  if (!component %in% names(dly_fit_obj$yields)) {
    stop(sprintf("DLY fit does not contain yield component `%s`.", component))
  }

  list(
    model = "DLY",
    fit = dly_fit_obj,
    component = component,
    curve_panels = tontine_curve_yield_list_to_panels(
      curve_yields = dly_fit_obj$yields[[component]][curves],
      curves = curves,
      mat_names = mat_names
    )
  )
}

tontine_fit_mce_baseline_curves <- function(yields, curves, mats, mat_names,
                                            lambdas, cutoffs, reference_curve,
                                            controls = list()) {
  get_control <- function(name, default) {
    if (is.null(controls[[name]])) default else controls[[name]]
  }

  mce_fit_obj <- mc_fit_shock_baseline(
    yields = yields,
    lambdas = lambdas,
    cutoffs = cutoffs,
    reference = reference_curve,
    curves = curves,
    mats = mats,
    ECM_estim = get_control("ECM_estim", "ML"),
    ECM_type = get_control("ECM_type", "eigen"),
    ECM_alpha = get_control("ECM_alpha", 0.1),
    X_normalize = get_control("X_normalize", TRUE),
    X_trunc = get_control("X_trunc", FALSE),
    Xt_smooth = get_control("Xt_smooth", FALSE),
    smoother = get_control("smoother", NULL),
    hp_lambda = get_control("hp_lambda", 1600),
    CR_maxiter = get_control("CR_maxiter", 1000),
    CR_tol = get_control("CR_tol", 1e-8),
    CR_S0_shrink_diag = get_control("CR_S0_shrink_diag", 0),
    CR_check_every = get_control("CR_check_every", 1),
    CR_use_loglik = get_control("CR_use_loglik", TRUE),
    CR_store_path = get_control("CR_store_path", TRUE)
  )

  list(
    model = "MCE",
    fit = mce_fit_obj,
    curve_panels = tontine_split_curve_panel(
      panel_df = mce_fit_obj$baseline_fit$curve,
      curves = curves,
      mat_names = mat_names
    )
  )
}

tontine_fit_baseline_curve_models <- function(yields, curves, mats, mat_names,
                                              lambdas, cutoffs, reference_curve,
                                              mce_controls = list(),
                                              include_models = c(
                                                "Observed", "DNS",
                                                "DLY", "MCE"
                                              )) {
  model_results <- list()

  if ("Observed" %in% include_models) {
    model_results$Observed <- list(
      model = "Observed",
      fit = NULL,
      curve_panels = tontine_observed_curve_panels(yields, curves, mat_names)
    )
  }

  if ("DNS" %in% include_models) {
    model_results$DNS <- tontine_fit_dns_baseline_curves(
      yields = yields,
      curves = curves,
      mats = mats,
      mat_names = mat_names,
      lambdas = lambdas,
      cutoffs = cutoffs
    )
  }

  if ("DLY" %in% include_models) {
    model_results$DLY <- tontine_fit_dly_baseline_curves(
      yields = yields,
      curves = curves,
      mats = mats,
      mat_names = mat_names,
      lambdas = lambdas,
      reference_curve = reference_curve
    )
  }

  if ("MCE" %in% include_models) {
    model_results$MCE <- tontine_fit_mce_baseline_curves(
      yields = yields,
      curves = curves,
      mats = mats,
      mat_names = mat_names,
      lambdas = lambdas,
      cutoffs = cutoffs,
      reference_curve = reference_curve,
      controls = mce_controls
    )
  }

  model_results
}

tontine_extract_purchase_curve_rows <- function(curve_panels, curves, mat_names,
                                                purchase_date) {
  setNames(lapply(curves, function(curve_name) {
    curve_panel <- curve_panels[[curve_name]]
    if (is.null(curve_panel)) {
      stop(sprintf("Curve panels are missing curve `%s`.", curve_name))
    }

    row_idx <- which(as.Date(curve_panel$time) == as.Date(purchase_date))
    if (length(row_idx) != 1) {
      stop("Curve panel must contain exactly one row for each window purchase date.")
    }

    out <- curve_panel[row_idx, c("time", mat_names), drop = FALSE]
    out$time <- as.Date(out$time)
    out
  }), curves)
}

tontine_extract_window_curve_rows <- function(curve_panels, curves, mat_names,
                                              valuation_window) {
  if (nrow(valuation_window) != 1) {
    stop("`valuation_window` must contain exactly one row.")
  }

  row_dates <- unique(as.Date(c(
    valuation_window$purchase_date[[1]],
    if ("valuation_date" %in% colnames(valuation_window)) {
      valuation_window$valuation_date[[1]]
    } else {
      valuation_window$purchase_date[[1]]
    }
  )))

  setNames(lapply(curves, function(curve_name) {
    curve_panel <- curve_panels[[curve_name]]
    if (is.null(curve_panel)) {
      stop(sprintf("Curve panels are missing curve `%s`.", curve_name))
    }

    out <- do.call(rbind, lapply(row_dates, function(row_date) {
      row_idx <- which(as.Date(curve_panel$time) == row_date)
      if (length(row_idx) != 1) {
        stop("Curve panel must contain exactly one row for each window report date.")
      }

      curve_row <- curve_panel[row_idx, c("time", mat_names), drop = FALSE]
      curve_row$time <- as.Date(curve_row$time)
      curve_row
    }))
    rownames(out) <- NULL
    out[order(as.Date(out$time)), , drop = FALSE]
  }), curves)
}

tontine_bind_purchase_curve_rows <- function(purchase_rows, curves) {
  setNames(lapply(curves, function(curve_name) {
    out <- do.call(rbind, lapply(purchase_rows, `[[`, curve_name))
    out <- out[!duplicated(as.Date(out$time)), , drop = FALSE]
    rownames(out) <- NULL
    out[order(as.Date(out$time)), , drop = FALSE]
  }), curves)
}

tontine_fit_rolling_baseline_curve_models <- function(yields, valuation_windows,
                                                      curves, mats, mat_names,
                                                      lambdas, cutoffs,
                                                      reference_curve,
                                                      mce_controls = list(),
                                                      include_models = c(
                                                        "Observed", "DNS",
                                                        "DLY", "MCE"
                                                      ),
                                                      fit_seed = NULL) {
  if (!all(c("rung_id", "purchase_date", "estimation_start", "estimation_end") %in%
           colnames(valuation_windows))) {
    stop("`valuation_windows` must be created by `tontine_build_valuation_windows()`.")
  }

  fit_one_window <- function(model_name, window_yields) {
    if (model_name == "Observed") {
      return(list(
        model = "Observed",
        fit = NULL,
        curve_panels = tontine_observed_curve_panels(window_yields, curves, mat_names)
      ))
    }

    if (model_name == "DNS") {
      return(tontine_fit_dns_baseline_curves(
        yields = window_yields,
        curves = curves,
        mats = mats,
        mat_names = mat_names,
        lambdas = lambdas,
        cutoffs = cutoffs
      ))
    }

    if (model_name == "DLY") {
      return(tontine_fit_dly_baseline_curves(
        yields = window_yields,
        curves = curves,
        mats = mats,
        mat_names = mat_names,
        lambdas = lambdas,
        reference_curve = reference_curve
      ))
    }

    if (model_name == "MCE") {
      return(tontine_fit_mce_baseline_curves(
        yields = window_yields,
        curves = curves,
        mats = mats,
        mat_names = mat_names,
        lambdas = lambdas,
        cutoffs = cutoffs,
        reference_curve = reference_curve,
        controls = mce_controls
      ))
    }

    stop(sprintf("Unsupported baseline model `%s`.", model_name))
  }

  model_results <- list()

  for (model_name in include_models) {
    window_fits <- vector("list", nrow(valuation_windows))
    purchase_rows <- vector("list", nrow(valuation_windows))

    for (i in seq_len(nrow(valuation_windows))) {
      valuation_window <- valuation_windows[i, , drop = FALSE]
      window_yields <- tontine_slice_yields_for_window(
        yields = yields,
        valuation_window = valuation_window,
        curves = curves
      )
      if (!is.null(fit_seed)) {
        set.seed(as.integer(fit_seed) + i)
      }
      window_fit <- fit_one_window(model_name, window_yields)
      window_fits[[i]] <- window_fit
      purchase_rows[[i]] <- tontine_extract_window_curve_rows(
        curve_panels = window_fit$curve_panels,
        curves = curves,
        mat_names = mat_names,
        valuation_window = valuation_window
      )
    }

    model_results[[model_name]] <- list(
      model = window_fits[[1]]$model,
      fit = window_fits,
      valuation_windows = valuation_windows,
      curve_panels = tontine_bind_purchase_curve_rows(purchase_rows, curves)
    )
  }

  model_results
}

tontine_is_decay_shock <- function(shock_spec) {
  identical(shock_spec$shock_profile, "linear_decay") ||
    identical(shock_spec$shock_profile, "linear_increase") ||
    (!is.null(shock_spec$shock_start_year) && !is.null(shock_spec$shock_end_year))
}

tontine_resolve_shock_anchor_date <- function(valuation_windows, year_value, anchor) {
  anchor <- match.arg(
    anchor,
    c("purchase_date", "valuation_date", "maturity_date", "estimation_start", "estimation_end")
  )
  window <- valuation_windows[valuation_windows$purchase_year == as.integer(year_value), , drop = FALSE]

  if (nrow(window) != 1) {
    stop(sprintf("Could not resolve a unique valuation window for shock year %s.", year_value))
  }

  as.Date(window[[anchor]][[1]])
}

tontine_resolve_decay_shock_dates <- function(valuation_windows, shock_spec) {
  start_anchor <- if (is.null(shock_spec$shock_start_anchor)) "purchase_date" else shock_spec$shock_start_anchor
  end_anchor <- if (is.null(shock_spec$shock_end_anchor)) "purchase_date" else shock_spec$shock_end_anchor

  start_date <- tontine_resolve_shock_anchor_date(
    valuation_windows = valuation_windows,
    year_value = shock_spec$shock_start_year,
    anchor = start_anchor
  )
  end_date <- tontine_resolve_shock_anchor_date(
    valuation_windows = valuation_windows,
    year_value = shock_spec$shock_end_year,
    anchor = end_anchor
  )

  if (end_date <= start_date) {
    stop("Decay shock end date must be after its start date.")
  }

  list(start_date = start_date, end_date = end_date)
}

tontine_build_calendar_decay_profile <- function(dates, shock_start_date, shock_end_date,
                                                 shock_start_magnitude = 1,
                                                 shock_end_magnitude = 0,
                                                 profile = c("linear_decay", "linear_increase")) {
  profile <- match.arg(profile)
  dates <- as.Date(dates)
  shock_start_date <- as.Date(shock_start_date)
  shock_end_date <- as.Date(shock_end_date)

  if (!is.numeric(shock_start_magnitude) || length(shock_start_magnitude) != 1 ||
      is.na(shock_start_magnitude) ||
      !is.numeric(shock_end_magnitude) || length(shock_end_magnitude) != 1 ||
      is.na(shock_end_magnitude)) {
    stop("Decay shock magnitudes must be numeric scalars.")
  }

  magnitude <- rep(0, length(dates))
  in_decay <- dates >= shock_start_date & dates <= shock_end_date

  if (any(in_decay)) {
    elapsed <- as.numeric(dates[in_decay] - shock_start_date)
    duration <- as.numeric(shock_end_date - shock_start_date)
    magnitude[in_decay] <- shock_start_magnitude +
      (shock_end_magnitude - shock_start_magnitude) * (elapsed / duration)
  }

  data.frame(
    time = dates,
    shock_magnitude = magnitude,
    stringsAsFactors = FALSE
  )
}

tontine_build_full_window_decay_profile <- function(dates, valuation_date,
                                                    shock_start_date, shock_end_date,
                                                    shock_start_magnitude = 1,
                                                    shock_end_magnitude = 0,
                                                    profile = c("linear_decay", "linear_increase")) {
  profile <- match.arg(profile)
  dates <- as.Date(dates)
  valuation_date <- as.Date(valuation_date)

  valuation_profile <- tontine_build_calendar_decay_profile(
    dates = valuation_date,
    shock_start_date = shock_start_date,
    shock_end_date = shock_end_date,
    shock_start_magnitude = shock_start_magnitude,
    shock_end_magnitude = shock_end_magnitude,
    profile = profile
  )

  data.frame(
    time = dates,
    shock_magnitude = rep(valuation_profile$shock_magnitude[[1]], length(dates)),
    stringsAsFactors = FALSE
  )
}

tontine_profile_has_shock <- function(shock_profile) {
  if (is.null(shock_profile)) {
    return(FALSE)
  }

  any(abs(as.numeric(shock_profile$shock_magnitude)) > sqrt(.Machine$double.eps))
}

tontine_uses_average_level_fraction_shock <- function(shock_spec) {
  identical(shock_spec$shock_basis, "average_level_fraction")
}

tontine_resolve_average_level_shock_profile <- function(curve_factor, shock_profile,
                                                        shock_factor = "L",
                                                        curves = names(curve_factor)) {
  if (!is.data.frame(shock_profile) || !"shock_magnitude" %in% colnames(shock_profile)) {
    stop("`shock_profile` must contain `shock_magnitude`.")
  }

  out <- shock_profile
  out$shock_fraction <- as.numeric(out$shock_magnitude)
  out$shock_magnitude <- st_average_level_shock_profile(
    curve_factor = curve_factor,
    shock_fraction_profile = out$shock_fraction,
    shock_factor = shock_factor,
    curves = curves
  )
  out
}

tontine_mce_shock_from_baseline_profile <- function(baseline_state, shock_spec,
                                                    shock_profile,
                                                    shock_type = c("additive", "multiplicative"),
                                                    controls = list()) {
  shock_type <- match.arg(shock_type)
  baseline_state <- st_validate_baseline_state(baseline_state)
  config <- baseline_state$config
  single_curve_state <- baseline_state$single_curve_state

  get_control <- function(name, default) {
    if (is.null(controls[[name]])) default else controls[[name]]
  }

  reuse_BS0 <- get_control("reuse_BS0", TRUE)
  reuse_ECM <- get_control("reuse_ECM", TRUE)
  ECM_fallback <- match.arg(get_control("ECM_fallback", "baseline"), c("error", "baseline"))

  if (tontine_uses_average_level_fraction_shock(shock_spec)) {
    shock_profile <- tontine_resolve_average_level_shock_profile(
      curve_factor = single_curve_state$ns_factor,
      shock_profile = shock_profile,
      shock_factor = shock_spec$shock_factors,
      curves = baseline_state$curves
    )
  }

  shock_fit <- st_apply_factor_shock_profile(
    curves = baseline_state$curves,
    curve_dns_factor = single_curve_state$ns_factor,
    shock_curves = shock_spec$shock_curves,
    shock_factors = shock_spec$shock_factors,
    shock_magnitudes = shock_profile$shock_magnitude,
    shock_type = shock_type
  )

  shocked_state <- st_build_shocked_curve_state(
    single_curve_state = single_curve_state,
    shocked_factors = shock_fit$shocked_factors
  )

  shocked_xt_state <- st_build_xt_state(
    bl_yields = shocked_state$bl_yields,
    curves = baseline_state$curves,
    reference = baseline_state$reference,
    reuse_ECM = reuse_ECM,
    baseline_ecm = baseline_state$xt_state$cc_ci_features$cc_ecm,
    ECM_estim = st_resolve_config(get_control("ECM_estim", NULL), config$ECM_estim),
    ECM_type = st_resolve_config(get_control("ECM_type", NULL), config$ECM_type),
    ECM_alpha = st_resolve_config(get_control("ECM_alpha", NULL), config$ECM_alpha),
    ECM_fallback = ECM_fallback,
    X_normalize = st_resolve_config(get_control("X_normalize", NULL), config$X_normalize),
    X_trunc = st_resolve_config(get_control("X_trunc", NULL), config$X_trunc),
    Xt_smooth = st_resolve_config(get_control("Xt_smooth", NULL), config$Xt_smooth),
    smoother = st_resolve_config(get_control("smoother", NULL), config$smoother),
    knot_count = st_resolve_config(get_control("knot_count", NULL), config$knot_count),
    k_count = st_resolve_config(get_control("k_count", NULL), config$k_count),
    k_pass = st_resolve_config(get_control("k_pass", NULL), config$k_pass),
    hp_lambda = st_resolve_config(get_control("hp_lambda", NULL), config$hp_lambda),
    henderson_k = st_resolve_config(get_control("henderson_k", NULL), config$henderson_k),
    ytime = single_curve_state$time
  )

  if (reuse_BS0) {
    shocked_Xt <- st_align_Xt_to_baseline(shocked_xt_state$Xt, baseline_state$xt_state$Xt)
    shocked_BS0 <- baseline_state$BS0
  } else {
    shocked_Xt <- shocked_xt_state$Xt
    shocked_BS0 <- blcc_fe(
      Xt = shocked_Xt,
      Wt = shocked_state$W,
      mats = baseline_state$mats,
      covreg = st_resolve_config(get_control("CR_algo", NULL), config$CR_algo),
      init = st_resolve_config(get_control("CR_init", NULL), config$CR_init),
      max_iter = st_resolve_config(get_control("CR_maxiter", NULL), config$CR_maxiter),
      tol = st_resolve_config(get_control("CR_tol", NULL), config$CR_tol),
      S0 = st_resolve_config(get_control("CR_S0init", NULL), config$CR_S0init),
      B = st_resolve_config(get_control("CR_Binit", NULL), config$CR_Binit),
      S0_shrink_diag = st_resolve_config(get_control("CR_S0_shrink_diag", NULL), config$CR_S0_shrink_diag),
      check_every = st_resolve_config(get_control("CR_check_every", NULL), config$CR_check_every),
      use_loglik = st_resolve_config(get_control("CR_use_loglik", NULL), config$CR_use_loglik),
      store_path = st_resolve_config(get_control("CR_store_path", NULL), config$CR_store_path),
      verb = st_resolve_config(get_control("CR_verb", NULL), config$CR_verb),
      term = st_resolve_config(get_control("CR_term", NULL), config$CR_term)
    )
  }

  shocked_fit_state <- st_finalize_shocked_fit(
    yields = shocked_state$yields,
    time = single_curve_state$time,
    phi_hat = shocked_state$phi_hat,
    Xt = shocked_Xt,
    BS0 = shocked_BS0,
    curves = baseline_state$curves,
    mats = baseline_state$mats,
    W = shocked_state$W
  )

  list(
    tenor = baseline_state$baseline_fit$tenor,
    curve = baseline_state$baseline_fit$curve,
    sigma_JT = shocked_fit_state$sigma_JT,
    Xt = shocked_fit_state$Xt,
    W = shocked_fit_state$W,
    ns_factor = shocked_state$ns_factor,
    shock_profile = shock_profile,
    shocked_factors = shock_fit$shocked_factors,
    shocked_yields = shocked_state$yields,
    shocked_tenor = shocked_fit_state$tenor,
    shocked_curve = shocked_fit_state$curve,
    shocked_single_curve_state = shocked_state,
    baseline_state = baseline_state
  )
}

tontine_fit_rolling_shocked_curve_models <- function(yields, valuation_windows,
                                                     curves, mats, mat_names,
                                                     lambdas, cutoffs,
                                                     reference_curve,
                                                     shock_spec,
                                                     mce_controls = list(),
                                                     mce_shock_controls = list(),
                                                     include_models = c(
                                                       "DNS", "DLY",
                                                       "MCE"
                                                     ),
                                                     fit_seed = NULL) {
  is_decay_shock <- tontine_is_decay_shock(shock_spec)
  if (!is_decay_shock && is.null(shock_spec$shock_year)) {
    stop("`shock_spec` must contain `shock_year`.")
  }

  if (!is.null(shock_spec$shock_window)) {
    stop("Application shocks now use row-wise profiles; `shock_window` is no longer supported here.")
  }

  shock_type <- if (is.null(shock_spec$shock_type)) "additive" else shock_spec$shock_type
  shock_label <- if (is.null(shock_spec$label)) {
    tontine_format_shock_label(shock_spec)
  } else {
    shock_spec$label
  }
  decay_dates <- if (is_decay_shock) {
    tontine_resolve_decay_shock_dates(valuation_windows, shock_spec)
  } else {
    NULL
  }

  build_window_shock_profile <- function(window_time, valuation_window) {
    if (!is_decay_shock) {
      return(data.frame(
        time = as.Date(window_time),
        shock_magnitude = rep(shock_spec$shock_magnitude, length(window_time)),
        stringsAsFactors = FALSE
      ))
    }

    start_magnitude <- if (is.null(shock_spec$shock_start_magnitude)) {
      shock_spec$shock_magnitude
    } else {
      shock_spec$shock_start_magnitude
    }
    end_magnitude <- if (is.null(shock_spec$shock_end_magnitude)) 0 else shock_spec$shock_end_magnitude
    shock_profile_name <- if (is.null(shock_spec$shock_profile)) "linear_decay" else shock_spec$shock_profile
    shock_effect_timing <- if (is.null(shock_spec$shock_effect_timing)) {
      "calendar"
    } else {
      shock_spec$shock_effect_timing
    }

    if (!shock_effect_timing %in% c("calendar", "full_window")) {
      stop("`shock_effect_timing` must be either `calendar` or `full_window`.")
    }

    if (shock_effect_timing == "full_window") {
      valuation_date <- if ("valuation_date" %in% colnames(valuation_window)) {
        valuation_window$valuation_date[[1]]
      } else {
        valuation_window$purchase_date[[1]]
      }

      return(tontine_build_full_window_decay_profile(
        dates = window_time,
        valuation_date = valuation_date,
        shock_start_date = decay_dates$start_date,
        shock_end_date = decay_dates$end_date,
        shock_start_magnitude = start_magnitude,
        shock_end_magnitude = end_magnitude,
        profile = shock_profile_name
      ))
    }

    tontine_build_calendar_decay_profile(
      dates = window_time,
      shock_start_date = decay_dates$start_date,
      shock_end_date = decay_dates$end_date,
      shock_start_magnitude = start_magnitude,
      shock_end_magnitude = end_magnitude,
      profile = shock_profile_name
    )
  }

  fit_one_window <- function(model_name, window_yields, valuation_window, window_index) {
    if (!is.null(fit_seed)) {
      set.seed(as.integer(fit_seed) + window_index)
    }

    shock_this_window <- if (is_decay_shock) {
      valuation_window$estimation_end[[1]] >= decay_dates$start_date &&
        valuation_window$estimation_start[[1]] <= decay_dates$end_date
    } else {
      valuation_window$purchase_year[[1]] == shock_spec$shock_year
    }
    if (model_name == "Observed") {
      stop("Observed curves are not shocked in the tontine shock applications.")
    }

    if (model_name == "DNS") {
      dns_fit <- tontine_fit_dns_baseline_curves(
        yields = window_yields,
        curves = curves,
        mats = mats,
        mat_names = mat_names,
        lambdas = lambdas,
        cutoffs = cutoffs
      )
      if (!shock_this_window) {
        return(dns_fit)
      }

      shock_profile <- build_window_shock_profile(dns_fit$fit$time, valuation_window)
      if (tontine_uses_average_level_fraction_shock(shock_spec)) {
        shock_profile <- tontine_resolve_average_level_shock_profile(
          curve_factor = dns_fit$fit$ns_factor,
          shock_profile = shock_profile,
          shock_factor = shock_spec$shock_factors,
          curves = dns_fit$fit$curves
        )
      }
      if (!tontine_profile_has_shock(shock_profile)) {
        return(dns_fit)
      }
      dns_profile_shock <- st_apply_factor_shock_profile(
        curves = dns_fit$fit$curves,
        curve_dns_factor = dns_fit$fit$ns_factor,
        shock_curves = shock_spec$shock_curves,
        shock_factors = shock_spec$shock_factors,
        shock_magnitudes = shock_profile$shock_magnitude,
        shock_type = shock_type
      )
      shocked_state <- st_build_shocked_curve_state(
        single_curve_state = dns_fit$fit,
        shocked_factors = dns_profile_shock$shocked_factors
      )
      dns_shock <- list(
        yields = shocked_state$yields,
        ns_factor = shocked_state$ns_factor,
        shock_profile = shock_profile,
        shocked_factors = dns_profile_shock$shocked_factors,
        baseline_state = dns_fit$fit
      )

      dns_fit$curve_panels <- tontine_curve_yield_list_to_panels(
        curve_yields = dns_shock$yields,
        curves = curves,
        mat_names = mat_names
      )
      dns_fit$shock_fit <- dns_shock
      return(dns_fit)
    }

    if (model_name == "DLY") {
      dly_fit_obj <- tontine_fit_dly_baseline_curves(
        yields = window_yields,
        curves = curves,
        mats = mats,
        mat_names = mat_names,
        lambdas = lambdas,
        reference_curve = reference_curve
      )
      if (!shock_this_window) {
        return(dly_fit_obj)
      }

      shock_profile <- build_window_shock_profile(
        dly_fit_obj$fit$single_curve[[curves[[1]]]]$time,
        valuation_window
      )
      if (tontine_uses_average_level_fraction_shock(shock_spec)) {
        shock_profile <- tontine_resolve_average_level_shock_profile(
          curve_factor = dly_fit_obj$fit$country_betas$raw,
          shock_profile = shock_profile,
          shock_factor = shock_spec$shock_factors,
          curves = dly_fit_obj$fit$curves
        )
      }
      if (!tontine_profile_has_shock(shock_profile)) {
        return(dly_fit_obj)
      }
      dly_shock <- if (tontine_uses_average_level_fraction_shock(shock_spec)) {
        st_dly_global_shock_from_fit_profile(
          dly_fit_obj = dly_fit_obj$fit,
          shock_factors = shock_spec$shock_factors,
          shock_profile = shock_profile,
          shock_type = shock_type
        )
      } else {
        st_dly_shock_from_fit_profile(
          dly_fit_obj = dly_fit_obj$fit,
          shock_curves = shock_spec$shock_curves,
          shock_factors = shock_spec$shock_factors,
          shock_profile = shock_profile,
          shock_type = shock_type
        )
      }

      dly_fit_obj$curve_panels <- tontine_curve_yield_list_to_panels(
        curve_yields = dly_shock$fit$yields[[dly_fit_obj$component]][curves],
        curves = curves,
        mat_names = mat_names
      )
      dly_fit_obj$shock_fit <- dly_shock
      return(dly_fit_obj)
    }

    if (model_name == "MCE") {
      mce_fit_obj <- tontine_fit_mce_baseline_curves(
        yields = window_yields,
        curves = curves,
        mats = mats,
        mat_names = mat_names,
        lambdas = lambdas,
        cutoffs = cutoffs,
        reference_curve = reference_curve,
        controls = mce_controls
      )
      if (!shock_this_window) {
        return(mce_fit_obj)
      }

      shock_profile <- build_window_shock_profile(
        mce_fit_obj$fit$single_curve_state$time,
        valuation_window
      )
      if (!tontine_profile_has_shock(shock_profile)) {
        return(mce_fit_obj)
      }
      mce_shock <- tontine_mce_shock_from_baseline_profile(
        baseline_state = mce_fit_obj$fit,
        shock_spec = shock_spec,
        shock_profile = shock_profile,
        shock_type = shock_type,
        controls = mce_shock_controls
      )

      mce_fit_obj$curve_panels <- tontine_split_curve_panel(
        panel_df = mce_shock$shocked_curve,
        curves = curves,
        mat_names = mat_names
      )
      mce_fit_obj$shock_fit <- mce_shock
      return(mce_fit_obj)
    }

    stop(sprintf("Unsupported baseline model `%s`.", model_name))
  }

  model_results <- list()

  for (model_name in include_models) {
    window_fits <- vector("list", nrow(valuation_windows))
    purchase_rows <- vector("list", nrow(valuation_windows))

    for (i in seq_len(nrow(valuation_windows))) {
      valuation_window <- valuation_windows[i, , drop = FALSE]
      window_yields <- tontine_slice_yields_for_window(
        yields = yields,
        valuation_window = valuation_window,
        curves = curves
      )
      window_fit <- fit_one_window(model_name, window_yields, valuation_window, i)
      window_fits[[i]] <- window_fit
      purchase_rows[[i]] <- tontine_extract_window_curve_rows(
        curve_panels = window_fit$curve_panels,
        curves = curves,
        mat_names = mat_names,
        valuation_window = valuation_window
      )
    }

    model_results[[model_name]] <- list(
      model = window_fits[[1]]$model,
      fit = window_fits,
      valuation_windows = valuation_windows,
      shock_label = shock_label,
      shock_spec = shock_spec,
      curve_panels = tontine_bind_purchase_curve_rows(purchase_rows, curves)
    )
  }

  model_results
}

tontine_price_baseline_models <- function(baseline_curve_models, ladder_schedule,
                                          mat_names, scenario_type = "baseline",
                                          scenario_name_prefix = "Baseline",
                                          shock_label = NA_character_,
                                          initial_cash = NULL) {
  country_holdings <- do.call(rbind, lapply(names(baseline_curve_models), function(model_name) {
    model_obj <- baseline_curve_models[[model_name]]

    tontine_price_bond_ladder(
      curve_panels = model_obj$curve_panels,
      ladder_schedule = ladder_schedule,
      mat_names = mat_names,
      model_name = model_obj$model,
      scenario_name = paste(scenario_name_prefix, model_obj$model),
      scenario_type = scenario_type,
      shock_label = shock_label
    )
  }))
  rownames(country_holdings) <- NULL

  if (is.null(initial_cash)) {
    initial_cash <- stats::setNames(
      vapply(split(country_holdings, country_holdings$model), tontine_resolve_initial_cash, numeric(1)),
      names(split(country_holdings, country_holdings$model))
    )
  } else {
    priced_models <- unique(country_holdings$model)
    missing_initial_cash <- setdiff(priced_models, names(initial_cash))
    if (length(missing_initial_cash) > 0) {
      stop("`initial_cash` is missing one or more priced models.")
    }
    initial_cash <- initial_cash[priced_models]
  }

  cash_accounts <- do.call(rbind, lapply(names(initial_cash), function(model_name) {
    tontine_apply_ladder_cash_account(
      country_holdings = country_holdings[country_holdings$model == model_name, , drop = FALSE],
      initial_cash = initial_cash[[model_name]]
    )
  }))
  rownames(cash_accounts) <- NULL

  list(
    country_holdings = country_holdings,
    cash_accounts = cash_accounts,
    initial_cash = initial_cash
  )
}

tontine_summarize_baseline_ladder <- function(baseline_valuation) {
  country_holdings <- baseline_valuation$country_holdings
  cash_accounts <- baseline_valuation$cash_accounts

  total_purchase <- stats::aggregate(
    purchase_cost ~ model,
    data = country_holdings,
    FUN = sum
  )
  names(total_purchase)[names(total_purchase) == "purchase_cost"] <- "total_purchase_cost"

  avg_rate <- stats::aggregate(
    zero_rate ~ model,
    data = country_holdings,
    FUN = mean
  )
  names(avg_rate)[names(avg_rate) == "zero_rate"] <- "average_zero_rate"
  avg_rate$average_1y_zero_rate <- avg_rate$average_zero_rate

  ending_cash <- stats::aggregate(
    cash_end ~ model,
    data = cash_accounts,
    FUN = function(x) x[[length(x)]]
  )
  names(ending_cash)[names(ending_cash) == "cash_end"] <- "ending_cash_balance"

  min_cash <- stats::aggregate(
    cash_end ~ model,
    data = cash_accounts,
    FUN = min
  )
  names(min_cash)[names(min_cash) == "cash_end"] <- "minimum_cash_balance"

  max_reserve <- stats::aggregate(
    cash_end ~ model,
    data = cash_accounts,
    FUN = function(x) max(pmax(-x, 0))
  )
  names(max_reserve)[names(max_reserve) == "cash_end"] <- "maximum_reserve_draw"

  total_profit_loss <- stats::aggregate(
    rung_profit_loss ~ model,
    data = cash_accounts,
    FUN = sum
  )
  names(total_profit_loss)[names(total_profit_loss) == "rung_profit_loss"] <- "total_profit_loss"

  if ("mortality_credit_par" %in% colnames(cash_accounts) &&
      any(!is.na(cash_accounts$mortality_credit_par))) {
    total_mortality_credit <- stats::aggregate(
      mortality_credit_par ~ model,
      data = cash_accounts,
      FUN = sum,
      na.rm = TRUE
    )
    names(total_mortality_credit)[names(total_mortality_credit) == "mortality_credit_par"] <- "total_mortality_credit"
    total_mortality_credit$distributable_tontine_gain <- total_mortality_credit$total_mortality_credit
  } else {
    total_mortality_credit <- data.frame(
      model = unique(cash_accounts$model),
      total_mortality_credit = NA_real_,
      distributable_tontine_gain = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  initial_cash <- data.frame(
    model = names(baseline_valuation$initial_cash),
    initial_cash = as.numeric(baseline_valuation$initial_cash),
    stringsAsFactors = FALSE
  )

  Reduce(function(x, y) merge(x, y, by = "model", all = TRUE), list(
    initial_cash,
    total_purchase,
    ending_cash,
    min_cash,
    max_reserve,
    total_profit_loss,
    total_mortality_credit,
    avg_rate
  ))
}

tontine_add_calendar_years <- function(date_value, years) {
  date_value <- as.Date(date_value)
  date_parts <- as.POSIXlt(date_value)
  date_parts$year <- date_parts$year + as.integer(years)
  as.Date(date_parts)
}

tontine_resolve_maturity_date <- function(reference_panel, purchase_date, maturity_year) {
  panel_time <- as.Date(reference_panel$time)
  year_idx <- which(format(panel_time, "%Y") == as.character(as.integer(maturity_year)))

  if (length(year_idx) > 0) {
    return(min(panel_time[year_idx]))
  }

  tontine_add_calendar_years(purchase_date, maturity_year - as.integer(format(purchase_date, "%Y")))
}

tontine_build_rolling_reinvestment_ladder <- function(curve_panels, yields, reference_curve,
                                                      purchase_years, curves, initial_cash,
                                                      q, population_simulation = NULL,
                                                      bond_maturity_years = 1, mat_names,
                                                      model_name, scenario_name,
                                                      scenario_type = "baseline",
                                                      shock_label = NA_character_) {
  if (!is.list(curve_panels) || is.null(names(curve_panels))) {
    stop("`curve_panels` must be a named list of curve data frames.")
  }

  if (!is.list(yields) || is.null(names(yields)) || !reference_curve %in% names(yields)) {
    stop("`yields` must be a named list containing `reference_curve`.")
  }

  if (!is.character(curves) || length(curves) == 0 || anyNA(curves) || any(!nzchar(curves))) {
    stop("`curves` must be a non-empty character vector.")
  }

  missing_curves <- setdiff(curves, names(curve_panels))
  if (length(missing_curves) > 0) {
    stop("`curve_panels` is missing one or more ladder countries.")
  }

  if (!is.numeric(initial_cash) || length(initial_cash) != 1 ||
      is.na(initial_cash) || initial_cash <= 0) {
    stop("`initial_cash` must be a positive numeric scalar.")
  }

  purchase_years <- as.integer(purchase_years)
  bond_maturity_years <- tontine_validate_bond_maturity_years(bond_maturity_years)
  q <- tontine_validate_probability_vector(q, purchase_years)
  p_k <- cumprod(1 - q)
  tenor_col <- tontine_bond_tenor_col(bond_maturity_years)

  if (!tenor_col %in% mat_names) {
    stop(sprintf("`mat_names` must include required ladder tenor column `%s`.", tenor_col))
  }

  if (!is.null(population_simulation)) {
    required_population_cols <- c(
      "rung_id", "q_k", "expected_p_k", "population_start",
      "deaths", "population_end", "realized_p_k"
    )

    if (!all(required_population_cols %in% colnames(population_simulation))) {
      stop("`population_simulation` must be created by `tontine_simulate_population()`.")
    }

    if (nrow(population_simulation) != length(purchase_years)) {
      stop("`population_simulation` must contain one row per purchase year.")
    }

    if (!all(population_simulation$rung_id == seq_along(purchase_years))) {
      stop("`population_simulation$rung_id` must align with `purchase_years`.")
    }

    if (any(abs(population_simulation$q_k - q) > sqrt(.Machine$double.eps))) {
      stop("`population_simulation$q_k` must match `q`.")
    }
  }

  reference_panel <- yields[[reference_curve]]
  purchase_dates <- as.Date(vapply(purchase_years, function(year_value) {
    as.character(tontine_first_observed_date(reference_panel, year_value))
  }, character(1)))
  maturity_dates <- as.Date(vapply(seq_along(purchase_years), function(i) {
    as.character(tontine_resolve_maturity_date(
      reference_panel = reference_panel,
      purchase_date = purchase_dates[[i]],
      maturity_year = purchase_years[[i]] + bond_maturity_years
    ))
  }, character(1)))

  country_holdings <- vector("list", length(purchase_years))
  cash_accounts <- vector("list", length(purchase_years))
  initial_rung_count <- min(bond_maturity_years, length(purchase_years))
  annual_initial_cash <- initial_cash / initial_rung_count
  reserve_cash_start <- initial_cash
  carry_cash_start <- 0
  cumulative_profit_loss <- 0

  for (i in seq_along(purchase_years)) {
    purchase_date <- purchase_dates[[i]]
    maturity_date <- maturity_dates[[i]]
    initial_cash_draw <- if (i <= initial_rung_count) {
      min(annual_initial_cash, reserve_cash_start)
    } else {
      0
    }
    maturity_proceeds_from_prior <- if (i == 1) {
      0
    } else {
      prior_cash_accounts <- do.call(rbind, cash_accounts[seq_len(i - 1)])
      previous_purchase_date <- purchase_dates[[i - 1]]
      sum(
        prior_cash_accounts$maturity_proceeds[
          as.Date(prior_cash_accounts$maturity_date) > previous_purchase_date &
            as.Date(prior_cash_accounts$maturity_date) <= purchase_date
        ],
        na.rm = TRUE
      )
    }
    available_cash <- carry_cash_start + initial_cash_draw + maturity_proceeds_from_prior

    zero_rates <- vapply(curves, function(country_name) {
      tontine_extract_curve_rate(
        curve_panel = curve_panels[[country_name]],
        target_date = purchase_date,
        tenor_col = tenor_col
      ) / 100
    }, numeric(1))
    discount_factors <- exp(-zero_rates * bond_maturity_years)
    total_purchased_par <- available_cash / mean(discount_factors)
    country_par <- total_purchased_par / length(curves)
    purchase_costs <- country_par * discount_factors
    rung_purchase_cost <- sum(purchase_costs)
    rung_maturity_proceeds <- total_purchased_par
    cash_end <- available_cash - rung_purchase_cost
    reserve_cash_end <- reserve_cash_start - initial_cash_draw
    rung_profit_loss <- rung_maturity_proceeds - rung_purchase_cost
    cumulative_profit_loss <- cumulative_profit_loss + rung_profit_loss

    population_start <- if (is.null(population_simulation)) {
      NA_integer_
    } else {
      population_simulation$population_start[[i]]
    }
    deaths <- if (is.null(population_simulation)) {
      NA_integer_
    } else {
      population_simulation$deaths[[i]]
    }
    population_end <- if (is.null(population_simulation)) {
      NA_integer_
    } else {
      population_simulation$population_end[[i]]
    }
    realized_p_k <- if (is.null(population_simulation)) {
      NA_real_
    } else {
      population_simulation$realized_p_k[[i]]
    }
    purchased_par_per_start_member <- if (!is.na(population_start) && population_start > 0) {
      total_purchased_par / population_start
    } else {
      NA_real_
    }
    per_survivor_maturity_value <- if (!is.na(population_end) && population_end > 0) {
      total_purchased_par / population_end
    } else {
      NA_real_
    }
    mortality_credit_par <- if (!is.na(deaths) && !is.na(purchased_par_per_start_member)) {
      deaths * purchased_par_per_start_member
    } else {
      NA_real_
    }
    mortality_credit_per_survivor <- per_survivor_maturity_value - purchased_par_per_start_member

    country_holdings[[i]] <- data.frame(
      model = model_name,
      scenario = scenario_name,
      scenario_type = scenario_type,
      shock_label = shock_label,
      rung_id = i,
      purchase_year = purchase_years[[i]],
      purchase_date = purchase_date,
      maturity_date = maturity_date,
      bond_maturity_years = bond_maturity_years,
      tenor_col = tenor_col,
      country = curves,
      target_par = total_purchased_par,
      country_par = country_par,
      q_k = q[[i]],
      p_k = p_k[[i]],
      target_rung_par_per_member = purchased_par_per_start_member,
      population_start = population_start,
      deaths = deaths,
      population_end = population_end,
      realized_p_k = realized_p_k,
      next_rung_required_par = NA_real_,
      mortality_credit_par = mortality_credit_par / length(curves),
      per_survivor_maturity_value = per_survivor_maturity_value,
      mortality_credit_per_survivor = mortality_credit_per_survivor,
      survival_adjusted_par = country_par,
      purchase_par = country_par,
      priced_par = country_par,
      zero_rate = as.numeric(zero_rates[curves]),
      zero_rate_1y = as.numeric(zero_rates[curves]),
      discount_factor = as.numeric(discount_factors[curves]),
      purchase_cost = as.numeric(purchase_costs[curves]),
      maturity_proceeds = country_par,
      stringsAsFactors = FALSE
    )

    cash_accounts[[i]] <- data.frame(
      model = model_name,
      scenario = scenario_name,
      scenario_type = scenario_type,
      shock_label = shock_label,
      rung_id = i,
      purchase_year = purchase_years[[i]],
      purchase_date = purchase_date,
      maturity_date = maturity_date,
      bond_maturity_years = bond_maturity_years,
      initial_cash = initial_cash,
      initial_rung_count = initial_rung_count,
      reserve_cash_start = reserve_cash_start,
      initial_cash_draw = initial_cash_draw,
      reserve_cash_end = reserve_cash_end,
      cash_start = carry_cash_start,
      maturity_proceeds_from_prior = maturity_proceeds_from_prior,
      available_cash = available_cash,
      purchase_cost = rung_purchase_cost,
      maturity_proceeds = rung_maturity_proceeds,
      purchased_par = total_purchased_par,
      target_rung_par_per_member = purchased_par_per_start_member,
      population_start = population_start,
      deaths = deaths,
      population_end = population_end,
      mortality_credit_par = mortality_credit_par,
      rung_profit_loss = rung_profit_loss,
      cumulative_profit_loss = cumulative_profit_loss,
      cash_end = cash_end,
      reserve_draw = pmax(-cash_end, 0),
      stringsAsFactors = FALSE
    )

    reserve_cash_start <- reserve_cash_end
    carry_cash_start <- cash_end
  }

  country_holdings <- do.call(rbind, country_holdings)
  cash_accounts <- do.call(rbind, cash_accounts)
  rownames(country_holdings) <- NULL
  rownames(cash_accounts) <- NULL

  list(
    country_holdings = country_holdings,
    cash_accounts = cash_accounts,
    initial_cash = stats::setNames(initial_cash, model_name)
  )
}

tontine_price_rolling_reinvestment_models <- function(baseline_curve_models, yields,
                                                      reference_curve, purchase_years,
                                                      curves, initial_cash, q,
                                                      population_simulation = NULL,
                                                      bond_maturity_years = 1,
                                                      mat_names,
                                                      scenario_type = "baseline",
                                                      scenario_name_prefix = "Baseline",
                                                      shock_label = NA_character_) {
  model_valuations <- lapply(names(baseline_curve_models), function(model_name) {
    model_obj <- baseline_curve_models[[model_name]]
    tontine_build_rolling_reinvestment_ladder(
      curve_panels = model_obj$curve_panels,
      yields = yields,
      reference_curve = reference_curve,
      purchase_years = purchase_years,
      curves = curves,
      initial_cash = initial_cash,
      q = q,
      population_simulation = population_simulation,
      bond_maturity_years = bond_maturity_years,
      mat_names = mat_names,
      model_name = model_obj$model,
      scenario_name = paste(scenario_name_prefix, model_obj$model),
      scenario_type = scenario_type,
      shock_label = shock_label
    )
  })

  country_holdings <- do.call(rbind, lapply(model_valuations, `[[`, "country_holdings"))
  cash_accounts <- do.call(rbind, lapply(model_valuations, `[[`, "cash_accounts"))
  initial_cash_by_model <- stats::setNames(
    rep(initial_cash, length(model_valuations)),
    vapply(model_valuations, function(valuation) {
      unique(valuation$cash_accounts$model)
    }, character(1))
  )
  rownames(country_holdings) <- NULL
  rownames(cash_accounts) <- NULL

  list(
    country_holdings = country_holdings,
    cash_accounts = cash_accounts,
    initial_cash = initial_cash_by_model
  )
}

tontine_build_annual_contribution_ladder <- function(curve_panels, yields, reference_curve,
                                                     purchase_years, curves,
                                                     annual_contribution_per_member,
                                                     q, population_simulation = NULL,
                                                     bond_maturity_years = 1,
                                                     mat_names, model_name,
                                                     scenario_name,
                                                     scenario_type = "baseline",
                                                     shock_label = NA_character_) {
  if (!is.list(curve_panels) || is.null(names(curve_panels))) {
    stop("`curve_panels` must be a named list of curve data frames.")
  }

  if (!is.list(yields) || is.null(names(yields)) || !reference_curve %in% names(yields)) {
    stop("`yields` must be a named list containing `reference_curve`.")
  }

  if (!is.character(curves) || length(curves) == 0 || anyNA(curves) || any(!nzchar(curves))) {
    stop("`curves` must be a non-empty character vector.")
  }

  missing_curves <- setdiff(curves, names(curve_panels))
  if (length(missing_curves) > 0) {
    stop("`curve_panels` is missing one or more ladder countries.")
  }

  if (!is.numeric(annual_contribution_per_member) ||
      length(annual_contribution_per_member) != 1 ||
      is.na(annual_contribution_per_member) ||
      annual_contribution_per_member <= 0) {
    stop("`annual_contribution_per_member` must be a positive numeric scalar.")
  }

  purchase_years <- as.integer(purchase_years)
  bond_maturity_years <- tontine_validate_bond_maturity_years(bond_maturity_years)
  q <- tontine_validate_probability_vector(q, purchase_years)
  p_k <- cumprod(1 - q)
  tenor_col <- tontine_bond_tenor_col(bond_maturity_years)

  if (!tenor_col %in% mat_names) {
    stop(sprintf("`mat_names` must include required ladder tenor column `%s`.", tenor_col))
  }

  if (!is.null(population_simulation)) {
    required_population_cols <- c(
      "rung_id", "q_k", "expected_p_k", "population_start",
      "deaths", "population_end", "realized_p_k"
    )

    if (!all(required_population_cols %in% colnames(population_simulation))) {
      stop("`population_simulation` must be created by `tontine_simulate_population()`.")
    }

    if (nrow(population_simulation) != length(purchase_years)) {
      stop("`population_simulation` must contain one row per purchase year.")
    }

    if (!all(population_simulation$rung_id == seq_along(purchase_years))) {
      stop("`population_simulation$rung_id` must align with `purchase_years`.")
    }

    if (any(abs(population_simulation$q_k - q) > sqrt(.Machine$double.eps))) {
      stop("`population_simulation$q_k` must match `q`.")
    }
  }

  reference_panel <- yields[[reference_curve]]
  purchase_dates <- as.Date(vapply(purchase_years, function(year_value) {
    as.character(tontine_first_observed_date(reference_panel, year_value))
  }, character(1)))
  maturity_dates <- as.Date(vapply(seq_along(purchase_years), function(i) {
    as.character(tontine_resolve_maturity_date(
      reference_panel = reference_panel,
      purchase_date = purchase_dates[[i]],
      maturity_year = purchase_years[[i]] + bond_maturity_years
    ))
  }, character(1)))

  country_holdings <- vector("list", length(purchase_years))
  payout_accounts <- vector("list", length(purchase_years))

  for (i in seq_along(purchase_years)) {
    purchase_date <- purchase_dates[[i]]
    maturity_date <- maturity_dates[[i]]
    population_start <- if (is.null(population_simulation)) {
      NA_integer_
    } else {
      population_simulation$population_start[[i]]
    }
    deaths <- if (is.null(population_simulation)) {
      NA_integer_
    } else {
      population_simulation$deaths[[i]]
    }
    population_end <- if (is.null(population_simulation)) {
      NA_integer_
    } else {
      population_simulation$population_end[[i]]
    }
    realized_p_k <- if (is.null(population_simulation)) {
      NA_real_
    } else {
      population_simulation$realized_p_k[[i]]
    }

    total_contribution <- population_start * annual_contribution_per_member
    country_contribution <- total_contribution / length(curves)
    zero_rates <- vapply(curves, function(country_name) {
      tontine_extract_curve_rate(
        curve_panel = curve_panels[[country_name]],
        target_date = purchase_date,
        tenor_col = tenor_col
      ) / 100
    }, numeric(1))
    discount_factors <- exp(-zero_rates * bond_maturity_years)
    country_par <- country_contribution / discount_factors
    total_purchase_cost <- sum(rep(country_contribution, length(curves)))
    total_maturity_proceeds <- sum(country_par)
    gross_maturity_value_per_start_member <- total_maturity_proceeds / population_start
    survivor_payout_per_member <- if (!is.na(population_end) && population_end > 0) {
      total_maturity_proceeds / population_end
    } else {
      NA_real_
    }
    investment_gain_per_start_member <-
      gross_maturity_value_per_start_member - annual_contribution_per_member
    mortality_credit_per_survivor <-
      survivor_payout_per_member - gross_maturity_value_per_start_member
    total_mortality_credit <- deaths * gross_maturity_value_per_start_member

    country_holdings[[i]] <- data.frame(
      model = model_name,
      scenario = scenario_name,
      scenario_type = scenario_type,
      shock_label = shock_label,
      rung_id = i,
      purchase_year = purchase_years[[i]],
      purchase_date = purchase_date,
      maturity_date = maturity_date,
      bond_maturity_years = bond_maturity_years,
      tenor_col = tenor_col,
      country = curves,
      q_k = q[[i]],
      p_k = p_k[[i]],
      annual_contribution_per_member = annual_contribution_per_member,
      population_start = population_start,
      deaths = deaths,
      population_end = population_end,
      realized_p_k = realized_p_k,
      total_contribution = total_contribution,
      country_contribution = country_contribution,
      country_par = as.numeric(country_par[curves]),
      survival_adjusted_par = as.numeric(country_par[curves]),
      purchase_par = as.numeric(country_par[curves]),
      priced_par = as.numeric(country_par[curves]),
      zero_rate = as.numeric(zero_rates[curves]),
      zero_rate_1y = as.numeric(zero_rates[curves]),
      discount_factor = as.numeric(discount_factors[curves]),
      purchase_cost = rep(country_contribution, length(curves)),
      maturity_proceeds = as.numeric(country_par[curves]),
      country_gross_maturity_value_per_start_member =
        as.numeric(country_par[curves]) / population_start,
      country_survivor_payout_per_member =
        as.numeric(country_par[curves]) / population_end,
      stringsAsFactors = FALSE
    )

    payout_accounts[[i]] <- data.frame(
      model = model_name,
      scenario = scenario_name,
      scenario_type = scenario_type,
      shock_label = shock_label,
      rung_id = i,
      purchase_year = purchase_years[[i]],
      purchase_date = purchase_date,
      maturity_date = maturity_date,
      bond_maturity_years = bond_maturity_years,
      q_k = q[[i]],
      p_k = p_k[[i]],
      annual_contribution_per_member = annual_contribution_per_member,
      population_start = population_start,
      deaths = deaths,
      population_end = population_end,
      realized_p_k = realized_p_k,
      total_contribution = total_contribution,
      purchase_cost = total_purchase_cost,
      maturity_proceeds = total_maturity_proceeds,
      bond_profit_loss = total_maturity_proceeds - total_purchase_cost,
      gross_maturity_value_per_start_member = gross_maturity_value_per_start_member,
      investment_gain_per_start_member = investment_gain_per_start_member,
      survivor_payout_per_member = survivor_payout_per_member,
      mortality_credit_per_survivor = mortality_credit_per_survivor,
      total_mortality_credit = total_mortality_credit,
      stringsAsFactors = FALSE
    )
  }

  country_holdings <- do.call(rbind, country_holdings)
  payout_accounts <- do.call(rbind, payout_accounts)
  rownames(country_holdings) <- NULL
  rownames(payout_accounts) <- NULL

  list(
    country_holdings = country_holdings,
    payout_accounts = payout_accounts
  )
}

tontine_price_annual_contribution_models <- function(baseline_curve_models, yields,
                                                     reference_curve,
                                                     purchase_years, curves,
                                                     annual_contribution_per_member,
                                                     q,
                                                     population_simulation = NULL,
                                                     bond_maturity_years = 1,
                                                     mat_names,
                                                     scenario_type = "baseline",
                                                     scenario_name_prefix = "Baseline",
                                                     shock_label = NA_character_) {
  model_valuations <- lapply(names(baseline_curve_models), function(model_name) {
    model_obj <- baseline_curve_models[[model_name]]
    tontine_build_annual_contribution_ladder(
      curve_panels = model_obj$curve_panels,
      yields = yields,
      reference_curve = reference_curve,
      purchase_years = purchase_years,
      curves = curves,
      annual_contribution_per_member = annual_contribution_per_member,
      q = q,
      population_simulation = population_simulation,
      bond_maturity_years = bond_maturity_years,
      mat_names = mat_names,
      model_name = model_obj$model,
      scenario_name = paste(scenario_name_prefix, model_obj$model),
      scenario_type = scenario_type,
      shock_label = shock_label
    )
  })

  country_holdings <- do.call(rbind, lapply(model_valuations, `[[`, "country_holdings"))
  payout_accounts <- do.call(rbind, lapply(model_valuations, `[[`, "payout_accounts"))
  rownames(country_holdings) <- NULL
  rownames(payout_accounts) <- NULL

  list(
    country_holdings = country_holdings,
    payout_accounts = payout_accounts
  )
}

#############################################
#        Specification Placeholders         #
#############################################

tontine_build_valuation_contexts <- function(...) {
  stop("TODO: define valuation context specification.")
}

tontine_build_schedule <- function(...) {
  stop("TODO: define tontine payout schedule specification.")
}

tontine_fit_baselines <- function(...) {
  stop("TODO: define baseline model fitting specification.")
}

tontine_fit_shocks <- function(...) {
  stop("TODO: define shocked model fitting specification.")
}

tontine_price_ladder <- function(...) {
  stop("TODO: define ladder pricing specification.")
}

tontine_summarize_results <- function(...) {
  stop("TODO: define result summary specification.")
}
