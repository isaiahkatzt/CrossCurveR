source_project_file(file.path(repo_root, "applications/tontine_helpers.R"))

make_ladder_yields <- function() {
  dates <- as.Date(c(
    "2020-01-02", "2020-12-31",
    "2021-01-04", "2021-12-31",
    "2022-01-03", "2022-12-30"
  ))
  curve_stub <- data.frame(
    time = dates,
    X01Y = c(2, 2.1, 3, 3.1, 4, 4.1),
    X02Y = c(2.2, 2.3, 3.2, 3.3, 4.2, 4.3)
  )

  list(
    usa = curve_stub,
    gbr = curve_stub
  )
}

test_that("tontine_build_bond_ladder_schedule uses observed year-start dates and survival probabilities", {
  schedule <- tontine_build_bond_ladder_schedule(
    yields = make_ladder_yields(),
    reference_curve = "usa",
    purchase_years = 2020:2021,
    curves = c("usa", "gbr"),
    total_target_par = 1000,
    q = c(0.01, 0.02)
  )

  expect_equal(nrow(schedule), 4)
  expect_equal(unique(schedule$rung_id), 1:2)
  expect_equal(unique(schedule$purchase_date), as.Date(c("2020-01-02", "2021-01-04")))
  expect_equal(unique(schedule$maturity_date), as.Date(c("2021-01-04", "2022-01-03")))
  expect_equal(unique(schedule$q_k), c(0.01, 0.02))
  expect_equal(unique(schedule$p_k), cumprod(1 - c(0.01, 0.02)), tolerance = tol)

  rung_par <- stats::aggregate(country_par ~ rung_id, data = schedule, sum)
  expect_equal(rung_par$country_par, c(1000, 1000))
  expect_equal(schedule$survival_adjusted_par, schedule$country_par * schedule$p_k, tolerance = tol)
})

test_that("tontine_build_bond_ladder_schedule validates mortality vector length and maturity dates", {
  expect_error(
    tontine_build_bond_ladder_schedule(
      yields = make_ladder_yields(),
      reference_curve = "usa",
      purchase_years = 2020:2021,
      curves = c("usa", "gbr"),
      total_target_par = 1000,
      q = 0.01
    ),
    "one value per purchase year"
  )

  expect_error(
    tontine_build_bond_ladder_schedule(
      yields = make_ladder_yields(),
      reference_curve = "usa",
      purchase_years = 2022,
      curves = c("usa", "gbr"),
      total_target_par = 1000,
      q = 0.01
    ),
    "No reference-curve observations found for year 2023"
  )
})

test_that("tontine_build_valuation_windows creates fixed lookback windows", {
  windows <- tontine_build_valuation_windows(
    yields = make_ladder_yields(),
    reference_curve = "usa",
    purchase_years = 2021,
    lookback_years = 1
  )

  expect_equal(windows$purchase_date, as.Date("2021-01-04"))
  expect_equal(windows$maturity_date, as.Date("2022-01-03"))
  expect_equal(windows$estimation_start, as.Date("2020-01-04"))
  expect_equal(windows$estimation_end, as.Date("2021-01-04"))
  expect_equal(windows$n_reference_observations, 2)

  sliced <- tontine_slice_yields_for_window(make_ladder_yields(), windows)
  expect_equal(sliced$usa$time, as.Date(c("2020-12-31", "2021-01-04")))

  expect_error(
    tontine_build_valuation_windows(
      yields = make_ladder_yields(),
      reference_curve = "usa",
      purchase_years = 2021,
      lookback_years = 1,
      min_observations = 3
    ),
    "fewer than 3 reference observations"
  )
})

test_that("calendar decay shocks resolve purchase-date anchors and decay to zero", {
  yields <- make_ladder_yields()
  windows <- tontine_build_valuation_windows(
    yields = yields,
    reference_curve = "usa",
    purchase_years = 2020:2021,
    lookback_years = 1
  )
  shock_spec <- list(
    shock_profile = "linear_decay",
    shock_start_year = 2020,
    shock_start_anchor = "purchase_date",
    shock_end_year = 2021,
    shock_end_anchor = "purchase_date"
  )

  decay_dates <- tontine_resolve_decay_shock_dates(windows, shock_spec)
  expect_equal(decay_dates$start_date, as.Date("2020-01-02"))
  expect_equal(decay_dates$end_date, as.Date("2021-01-04"))

  profile <- tontine_build_calendar_decay_profile(
    dates = as.Date(c("2019-12-31", "2020-01-02", "2021-01-04", "2021-01-05")),
    shock_start_date = decay_dates$start_date,
    shock_end_date = decay_dates$end_date,
    shock_start_magnitude = 1,
    shock_end_magnitude = 0
  )

  expect_equal(profile$shock_magnitude[[1]], 0)
  expect_equal(profile$shock_magnitude[[2]], 1)
  expect_equal(profile$shock_magnitude[[3]], 0)
  expect_equal(profile$shock_magnitude[[4]], 0)
})

test_that("tontine_simulate_population is reproducible and tracks realized survivors", {
  sim_1 <- tontine_simulate_population(N_members = 10, q = c(0, 1), seed = 42)
  sim_2 <- tontine_simulate_population(N_members = 10, q = c(0, 1), seed = 42)

  expect_equal(sim_1, sim_2)
  expect_equal(sim_1$population_start, c(10, 10))
  expect_equal(sim_1$deaths, c(0, 10))
  expect_equal(sim_1$population_end, c(10, 0))
  expect_equal(sim_1$realized_p_k, c(1, 0))
  expect_equal(sim_1$expected_p_k, c(1, 0))
})

test_that("tontine_build_bond_ladder_schedule can use one realized population simulation", {
  population_simulation <- tontine_simulate_population(N_members = 10, q = c(0, 1), seed = 42)
  schedule <- tontine_build_bond_ladder_schedule(
    yields = make_ladder_yields(),
    reference_curve = "usa",
    purchase_years = 2020:2021,
    curves = c("usa", "gbr"),
    total_target_par = 1000,
    q = c(0, 1),
    population_simulation = population_simulation,
    target_rung_par_per_member = 100
  )

  expect_equal(schedule$population_start, c(10, 10, 10, 10))
  expect_equal(schedule$population_end, c(10, 10, 0, 0))
  expect_equal(schedule$realized_p_k, c(1, 1, 0, 0))
  expect_equal(schedule$purchase_par, c(500, 500, 500, 500))
  expect_equal(schedule$mortality_credit_par, c(0, 0, 500, 500))
  expect_equal(schedule$next_rung_required_par, c(1000, 1000, 0, 0))
  expect_equal(schedule$per_survivor_maturity_value, c(100, 100, NA, NA))
  expect_equal(schedule$mortality_credit_per_survivor, c(0, 0, NA, NA))
  expect_equal(unique(schedule$target_rung_par_per_member), 100)
})

test_that("tontine_price_bond_ladder prices survival-adjusted 1Y ZCB holdings", {
  yields <- make_ladder_yields()
  yields$gbr$X01Y <- yields$gbr$X01Y + 1

  schedule <- tontine_build_bond_ladder_schedule(
    yields = yields,
    reference_curve = "usa",
    purchase_years = 2020,
    curves = c("usa", "gbr"),
    total_target_par = 1000,
    q = 0.01
  )

  holdings <- tontine_price_bond_ladder(
    curve_panels = yields,
    ladder_schedule = schedule,
    mat_names = c("X01Y", "X02Y"),
    model_name = "DNS",
    scenario_name = "Baseline DNS",
    scenario_type = "baseline"
  )

  expect_equal(holdings$zero_rate_1y, c(0.02, 0.03), tolerance = tol)
  expect_equal(holdings$discount_factor, exp(-c(0.02, 0.03)), tolerance = tol)
  expect_equal(holdings$survival_adjusted_par, c(495, 495), tolerance = tol)
  expect_equal(holdings$purchase_cost, c(495 * exp(-0.02), 495 * exp(-0.03)), tolerance = tol)
  expect_equal(holdings$maturity_proceeds, c(495, 495), tolerance = tol)
})

test_that("tontine_price_bond_ladder supports multi-year ZCB maturities", {
  yields <- make_ladder_yields()
  schedule <- tontine_build_bond_ladder_schedule(
    yields = yields,
    reference_curve = "usa",
    purchase_years = 2020,
    curves = c("usa", "gbr"),
    total_target_par = 1000,
    q = 0,
    bond_maturity_years = 2
  )

  holdings <- tontine_price_bond_ladder(
    curve_panels = yields,
    ladder_schedule = schedule,
    mat_names = c("X01Y", "X02Y"),
    model_name = "DNS",
    scenario_name = "Baseline DNS",
    scenario_type = "baseline"
  )

  expect_equal(unique(schedule$maturity_date), as.Date("2022-01-03"))
  expect_equal(unique(holdings$bond_maturity_years), 2)
  expect_equal(unique(holdings$tenor_col), "X02Y")
  expect_equal(holdings$zero_rate, c(0.022, 0.022), tolerance = tol)
  expect_equal(holdings$discount_factor, rep(exp(-0.022 * 2), 2), tolerance = tol)
  expect_equal(holdings$purchase_cost, rep(500 * exp(-0.022 * 2), 2), tolerance = tol)
})

test_that("tontine_compare_mtm_shock values outstanding multi-year holdings", {
  yields <- make_ladder_yields()
  schedule <- tontine_build_bond_ladder_schedule(
    yields = yields,
    reference_curve = "usa",
    purchase_years = 2020,
    curves = c("usa", "gbr"),
    total_target_par = 1000,
    q = 0,
    bond_maturity_years = 2
  )
  holdings <- tontine_price_bond_ladder(
    curve_panels = yields,
    ladder_schedule = schedule,
    mat_names = c("X01Y", "X02Y"),
    model_name = "DNS",
    scenario_name = "Baseline DNS",
    scenario_type = "baseline"
  )

  shocked_yields <- yields
  shocked_yields$usa$X01Y <- shocked_yields$usa$X01Y + 1
  shocked_yields$usa$X02Y <- shocked_yields$usa$X02Y + 1

  mtm_delta <- tontine_compare_mtm_shock(
    baseline_holdings = holdings,
    baseline_curve_models = list(DNS = list(model = "DNS", curve_panels = yields)),
    shocked_curve_models = list(DNS = list(model = "DNS", curve_panels = shocked_yields)),
    mats = c(12, 24),
    mat_names = c("X01Y", "X02Y"),
    valuation_date = as.Date("2021-01-04"),
    shock_label = "USA +100 bps"
  )

  expect_equal(nrow(mtm_delta), 2)
  expect_equal(unique(mtm_delta$purchase_year), 2020)
  expect_true(mtm_delta$mtm_value_delta[mtm_delta$country == "usa"] < 0)
  expect_equal(mtm_delta$mtm_value_delta[mtm_delta$country == "gbr"], 0, tolerance = tol)
})

test_that("tontine_resolve_initial_cash and cash account roll proceeds forward while retaining deficits", {
  yields <- make_ladder_yields()
  schedule <- tontine_build_bond_ladder_schedule(
    yields = yields,
    reference_curve = "usa",
    purchase_years = 2020:2021,
    curves = c("usa", "gbr"),
    total_target_par = 1000,
    q = c(0, 0)
  )
  holdings <- tontine_price_bond_ladder(
    curve_panels = yields,
    ladder_schedule = schedule,
    mat_names = c("X01Y", "X02Y"),
    model_name = "DNS",
    scenario_name = "Baseline DNS",
    scenario_type = "baseline"
  )

  initial_cash <- tontine_resolve_initial_cash(holdings)
  expect_equal(initial_cash, sum(holdings$purchase_cost[holdings$rung_id == 1]), tolerance = tol)

  cash_account <- tontine_apply_ladder_cash_account(
    country_holdings = holdings,
    initial_cash = initial_cash
  )

  expect_equal(cash_account$cash_start[[1]], initial_cash, tolerance = tol)
  expect_equal(cash_account$cash_end[[1]], 0, tolerance = tol)
  expect_equal(cash_account$maturity_proceeds_from_prior[[2]], 1000, tolerance = tol)
  expect_equal(
    cash_account$cash_end[[2]],
    1000 - sum(holdings$purchase_cost[holdings$rung_id == 2]),
    tolerance = tol
  )
  expect_equal(
    cash_account$rung_profit_loss,
    cash_account$maturity_proceeds - cash_account$purchase_cost,
    tolerance = tol
  )
  expect_equal(
    cash_account$cumulative_profit_loss,
    cumsum(cash_account$rung_profit_loss),
    tolerance = tol
  )

  deficit_account <- tontine_apply_ladder_cash_account(
    country_holdings = holdings,
    initial_cash = 0
  )
  expect_lt(deficit_account$cash_end[[1]], 0)
  expect_equal(deficit_account$reserve_draw, pmax(-deficit_account$cash_end, 0), tolerance = tol)
})

test_that("same ladder and survival schedule can be reused across shocked curve panels", {
  yields <- make_ladder_yields()
  schedule <- tontine_build_bond_ladder_schedule(
    yields = yields,
    reference_curve = "usa",
    purchase_years = 2020:2021,
    curves = c("usa", "gbr"),
    total_target_par = 1000,
    q = c(0.01, 0.02)
  )

  shocked_yields <- yields
  shocked_yields$usa$X01Y <- shocked_yields$usa$X01Y + 2

  baseline_holdings <- tontine_price_bond_ladder(
    curve_panels = yields,
    ladder_schedule = schedule,
    mat_names = c("X01Y", "X02Y"),
    model_name = "DNS",
    scenario_name = "Baseline DNS",
    scenario_type = "baseline"
  )
  shocked_holdings <- tontine_price_bond_ladder(
    curve_panels = shocked_yields,
    ladder_schedule = schedule,
    mat_names = c("X01Y", "X02Y"),
    model_name = "DNS",
    scenario_name = "DNS USA L +1",
    scenario_type = "shock",
    shock_label = "USA L +1"
  )

  expect_equal(shocked_holdings[, c("rung_id", "purchase_date", "maturity_date", "country", "country_par", "q_k", "p_k")],
               baseline_holdings[, c("rung_id", "purchase_date", "maturity_date", "country", "country_par", "q_k", "p_k")])
  expect_true(any(shocked_holdings$purchase_cost != baseline_holdings$purchase_cost))
  expect_equal(
    tontine_resolve_initial_cash(shocked_holdings),
    sum(shocked_holdings$purchase_cost[shocked_holdings$rung_id == 1]),
    tolerance = tol
  )
  expect_equal(
    tontine_resolve_initial_cash(baseline_holdings),
    sum(baseline_holdings$purchase_cost[baseline_holdings$rung_id == 1]),
    tolerance = tol
  )
})

test_that("baseline curve model fitter assembles observed and delegated baseline panels", {
  yields <- make_ladder_yields()
  mat_names <- c("X01Y", "X02Y")
  fake_panel <- yields

  result <- with_temp_bindings(
    list(
      tontine_fit_dns_baseline_curves = function(...) list(
        model = "DNS",
        fit = list(id = "dns"),
        curve_panels = fake_panel
      ),
      tontine_fit_dly_baseline_curves = function(...) list(
        model = "DLY",
        fit = list(id = "dly"),
        component = "dislocation_propagated",
        curve_panels = fake_panel
      ),
      tontine_fit_mce_baseline_curves = function(...) list(
        model = "MCE",
        fit = list(id = "mce"),
        curve_panels = fake_panel
      )
    ),
    tontine_fit_baseline_curve_models(
      yields = yields,
      curves = c("usa", "gbr"),
      mats = c(12, 24),
      mat_names = mat_names,
      lambdas = 0.1,
      cutoffs = c(1, 1),
      reference_curve = "usa",
      include_models = c("Observed", "DNS", "DLY", "MCE")
    )
  )

  expect_equal(names(result), c("Observed", "DNS", "DLY", "MCE"))
  expect_equal(result$Observed$model, "Observed")
  expect_equal(colnames(result$Observed$curve_panels$usa), c("time", mat_names))
  expect_equal(result$DNS$fit$id, "dns")
  expect_equal(result$DLY$fit$id, "dly")
  expect_equal(result$MCE$fit$id, "mce")
})

test_that("rolling baseline fitter refits models inside each valuation window", {
  yields <- make_ladder_yields()
  mat_names <- c("X01Y", "X02Y")
  windows <- tontine_build_valuation_windows(
    yields = yields,
    reference_curve = "usa",
    purchase_years = 2020:2021,
    lookback_years = 1
  )

  result <- with_temp_bindings(
    list(
      tontine_fit_dns_baseline_curves = function(yields, curves, mats, mat_names,
                                                 lambdas, cutoffs) {
        shifted <- yields
        for (curve_name in curves) {
          shifted[[curve_name]][, mat_names] <- shifted[[curve_name]][, mat_names] + 10
        }

        list(
          model = "DNS",
          fit = list(
            start = min(yields$usa$time),
            end = max(yields$usa$time)
          ),
          curve_panels = shifted
        )
      }
    ),
    tontine_fit_rolling_baseline_curve_models(
      yields = yields,
      valuation_windows = windows,
      curves = c("usa", "gbr"),
      mats = c(12, 24),
      mat_names = mat_names,
      lambdas = 0.1,
      cutoffs = c(1, 1),
      reference_curve = "usa",
      include_models = c("Observed", "DNS")
    )
  )

  expect_equal(names(result), c("Observed", "DNS"))
  expect_equal(result$Observed$curve_panels$usa$time, windows$purchase_date)
  expect_equal(result$DNS$curve_panels$usa$time, windows$purchase_date)
  expect_equal(result$DNS$curve_panels$usa$X01Y, c(12, 13), tolerance = tol)
  expect_equal(length(result$DNS$fit), 2)
  expect_equal(result$DNS$fit[[1]]$fit$start, as.Date("2020-01-02"))
  expect_equal(result$DNS$fit[[2]]$fit$start, as.Date("2020-12-31"))
})

test_that("baseline valuation prices each baseline model and builds cash-account summaries", {
  yields <- make_ladder_yields()
  shocked_yields <- yields
  shocked_yields$usa$X01Y <- shocked_yields$usa$X01Y + 1

  schedule <- tontine_build_bond_ladder_schedule(
    yields = yields,
    reference_curve = "usa",
    purchase_years = 2020:2021,
    curves = c("usa", "gbr"),
    total_target_par = 1000,
    q = c(0.01, 0.02)
  )

  baseline_curve_models <- list(
    Observed = list(model = "Observed", curve_panels = yields),
    DNS = list(model = "DNS", curve_panels = shocked_yields)
  )

  valuation <- tontine_price_baseline_models(
    baseline_curve_models = baseline_curve_models,
    ladder_schedule = schedule,
    mat_names = c("X01Y", "X02Y")
  )

  expect_equal(sort(unique(valuation$country_holdings$model)), c("DNS", "Observed"))
  expect_equal(sort(unique(valuation$cash_accounts$model)), c("DNS", "Observed"))
  expect_equal(
    valuation$initial_cash[["Observed"]],
    sum(valuation$country_holdings$purchase_cost[
      valuation$country_holdings$model == "Observed" &
        valuation$country_holdings$rung_id == 1
    ]),
    tolerance = tol
  )

  summary <- tontine_summarize_baseline_ladder(valuation)
  expect_true(all(c(
    "model", "initial_cash", "total_purchase_cost", "ending_cash_balance",
    "minimum_cash_balance", "maximum_reserve_draw", "total_profit_loss",
    "total_mortality_credit", "distributable_tontine_gain", "average_1y_zero_rate"
  ) %in% colnames(summary)))
  expect_equal(sort(summary$model), c("DNS", "Observed"))
  expect_true(all(is.finite(summary$total_purchase_cost)))
})

test_that("valuation accepts paired baseline initial cash with extra unpriced models", {
  yields <- make_ladder_yields()
  schedule <- tontine_build_bond_ladder_schedule(
    yields = yields,
    reference_curve = "usa",
    purchase_years = 2020:2021,
    curves = c("usa", "gbr"),
    total_target_par = 1000,
    q = c(0.01, 0.02)
  )

  valuation <- tontine_price_baseline_models(
    baseline_curve_models = list(
      DNS = list(model = "DNS", curve_panels = yields)
    ),
    ladder_schedule = schedule,
    mat_names = c("X01Y", "X02Y"),
    initial_cash = c(Observed = 999, DNS = 1000)
  )

  expect_equal(unique(valuation$country_holdings$model), "DNS")
  expect_equal(unique(valuation$cash_accounts$model), "DNS")
  expect_equal(names(valuation$initial_cash), "DNS")
  expect_equal(valuation$initial_cash[["DNS"]], 1000)
})
