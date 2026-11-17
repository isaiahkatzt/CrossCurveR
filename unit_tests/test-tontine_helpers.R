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

make_ladder_schedule <- function(yields, purchase_years, q, ...) {
  tontine_build_bond_ladder_schedule(yields, "usa", purchase_years,
    curves = names(yields), total_target_par = 1000, q = q, ...)
}

price_ladder <- function(yields, schedule, scenario_type = "baseline") {
  tontine_price_bond_ladder(yields, schedule, c("X01Y", "X02Y"),
    model_name = "DNS", scenario_name = scenario_type, scenario_type = scenario_type)
}

test_that("CSV reads preserve valuation panel names and rate lookups", {
  folder <- tempfile("curve-data-")
  dir.create(folder)
  withr::defer(unlink(folder, recursive = TRUE))
  expected <- make_ladder_yields()$usa
  utils::write.csv(expected, file.path(folder, "usa.csv"), row.names = FALSE)

  yields <- tontine_read_yields(list(usa = "usa.csv"), data_dir = folder)
  expect_equal(yields$usa, expected)
  expect_equal(tontine_extract_curve_rate(yields$usa, expected$time[1], "X01Y"), 2)
})

test_that("tontine_build_bond_ladder_schedule uses observed year-start dates and survival probabilities", {
  schedule <- tontine_build_bond_ladder_schedule(
    yields = make_ladder_yields(),
    reference_curve = "usa",
    purchase_years = 2020:2021,
    curves = c("usa", "gbr"),
    total_target_par = 1000,
    q = c(0.01, 0.02)
  )

  expect_equal(schedule$country, rep(c("usa", "gbr"), 2))
  expect_equal(schedule$purchase_date, rep(as.Date(c("2020-01-02", "2021-01-04")), each = 2))
  expect_equal(schedule$maturity_date, rep(as.Date(c("2021-01-04", "2022-01-03")), each = 2))
  expect_equal(schedule$p_k, rep(c(.99, .9702), each = 2), tolerance = tol)
  expect_equal(schedule$country_par, rep(500, 4))
  expect_equal(schedule$survival_adjusted_par, rep(c(495, 485.1), each = 2), tolerance = tol)
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

  expect_equal(profile$shock_magnitude, c(0, 1, 0, 0))
})

test_that("seeded population simulation tracks survivor counts", {
  withr::local_preserve_seed()
  sim <- tontine_simulate_population(N_members = 10, q = c(.2, .6), seed = 42)

  expect_equal(sim$deaths, c(1, 5))
  expect_equal(sim$population_end, c(9, 4))
  expect_equal(sim$realized_p_k, c(.9, .4))
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

  expect_equal(schedule$purchase_par, c(500, 500, 500, 500))
  expect_equal(schedule$mortality_credit_par, c(0, 0, 500, 500))
  expect_equal(schedule$next_rung_required_par, c(1000, 1000, 0, 0))
  expect_equal(schedule$per_survivor_maturity_value, c(100, 100, NA, NA))
  expect_equal(schedule$mortality_credit_per_survivor, c(0, 0, NA, NA))
})

test_that("tontine_price_bond_ladder prices survival-adjusted 1Y ZCB holdings", {
  yields <- make_ladder_yields()
  yields$gbr$X01Y <- yields$gbr$X01Y + 1

  schedule <- make_ladder_schedule(yields, 2020, q = 0.01)

  holdings <- price_ladder(yields, schedule)

  expect_equal(holdings$purchase_cost, 495 * exp(-c(.02, .03)), tolerance = tol)
  expect_equal(holdings$maturity_proceeds, c(495, 495), tolerance = tol)
})

test_that("tontine_price_bond_ladder supports multi-year ZCB maturities", {
  yields <- make_ladder_yields()
  schedule <- make_ladder_schedule(yields, 2020, q = 0, bond_maturity_years = 2)

  holdings <- price_ladder(yields, schedule)

  expect_equal(unique(schedule$maturity_date), as.Date("2022-01-03"))
  expect_equal(holdings$tenor_col, rep("X02Y", 2))
  expect_equal(holdings$purchase_cost, rep(500 * exp(-0.022 * 2), 2), tolerance = tol)
})

test_that("tontine_compare_mtm_shock values outstanding multi-year holdings", {
  yields <- make_ladder_yields()
  schedule <- make_ladder_schedule(yields, 2020, q = 0, bond_maturity_years = 2)
  holdings <- price_ladder(yields, schedule)

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
  schedule <- make_ladder_schedule(yields, 2020:2021, q = c(0, 0))
  holdings <- price_ladder(yields, schedule)

  initial_cash <- tontine_resolve_initial_cash(holdings)
  expect_equal(initial_cash, sum(holdings$purchase_cost[holdings$rung_id == 1]), tolerance = tol)

  cash_account <- tontine_apply_ladder_cash_account(
    country_holdings = holdings,
    initial_cash = initial_cash
  )

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

test_that("shock pricing uses the same ladder cash flows and new discount rates", {
  yields <- make_ladder_yields()
  schedule <- make_ladder_schedule(yields, 2020:2021, q = c(.01, .02))
  baseline <- price_ladder(yields, schedule)

  yields$usa$X01Y <- yields$usa$X01Y + 2
  shocked <- price_ladder(yields, schedule, scenario_type = "shock")

  expect_equal(shocked$maturity_proceeds, baseline$maturity_proceeds)
  expect_equal(shocked$purchase_cost / baseline$purchase_cost,
               rep(c(exp(-.02), 1), 2), tolerance = tol)
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
  starts <- lapply(result$DNS$fit, function(window) window$fit$start)
  expect_equal(unname(do.call(c, starts)), as.Date(c("2020-01-02", "2020-12-31")))
})

test_that("rolling shocks reuse baselines without changing refit draws or the saved fits", {
  yields <- make_ladder_yields()
  windows <- tontine_build_valuation_windows(yields, "usa", 2020:2021, lookback_years = 1)
  args <- list(yields = yields, valuation_windows = windows, curves = c("usa", "gbr"),
    mats = c(12, 24), mat_names = c("X01Y", "X02Y"), lambdas = .1,
    cutoffs = c(1, 1), reference_curve = "usa", include_models = "MCE", fit_seed = 999)
  calls <- 0L

  with_temp_bindings(list(
    tontine_fit_mce_baseline_curves = function(yields, ...) {
      calls <<- calls + 1L
      list(model = "MCE", curve_panels = yields,
        fit = list(single_curve_state = list(time = yields$usa$time), draw = stats::rnorm(1)))
    },
    tontine_mce_shock_from_baseline_profile = function(baseline_state, ...) {
      n <- length(baseline_state$single_curve_state$time)
      list(shocked_curve = data.frame(time = baseline_state$single_curve_state$time,
        usa.01Y = stats::rnorm(n), usa.02Y = 2, gbr.01Y = 3, gbr.02Y = 4))
    }
  ), {
    baseline <- do.call(tontine_fit_rolling_baseline_curve_models, args)
    saved <- unserialize(serialize(baseline, NULL))
    for (year in c(2020, 2019)) {
      for (magnitude in c(.1, 0)) {
        shock_args <- c(args, list(shock_spec = list(shock_year = year,
          shock_curves = "usa", shock_factors = "L", shock_magnitude = magnitude)))
        expected <- do.call(tontine_fit_rolling_shocked_curve_models, shock_args)
        expected_rng <- .Random.seed
        calls_before <- calls
        set.seed(123)
        actual <- do.call(tontine_fit_rolling_shocked_curve_models,
          c(shock_args, list(baseline_curve_models = baseline)))
        expect_identical(actual, expected)
        expect_identical(.Random.seed, expected_rng)
        expect_identical(calls, calls_before)
      }
    }
    expect_identical(baseline, saved)
    cached_args <- c(shock_args, list(baseline_curve_models = baseline))
    for (field in c("yields", "valuation_windows", "mce_controls", "fit_seed")) {
      incompatible <- cached_args
      incompatible[[field]] <- switch(field,
        yields = list(usa = transform(yields$usa, X01Y = X01Y + 1), gbr = yields$gbr),
        valuation_windows = windows[2:1, ], mce_controls = list(CR_maxiter = 3), fit_seed = 1000)
      expect_error(do.call(tontine_fit_rolling_shocked_curve_models, incompatible),
                   "same data, windows, settings, and seed")
    }
  })
})

test_that("baseline valuation prices each baseline model and builds cash-account summaries", {
  yields <- make_ladder_yields()
  shocked_yields <- yields
  shocked_yields$usa$X01Y <- shocked_yields$usa$X01Y + 1

  schedule <- make_ladder_schedule(yields, 2020:2021, q = c(0.01, 0.02))

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
  expect_equal(
    summary$total_purchase_cost,
    vapply(summary$model, function(model) {
      sum(valuation$country_holdings$purchase_cost[valuation$country_holdings$model == model])
    }, numeric(1), USE.NAMES = FALSE),
    tolerance = tol
  )
})

test_that("valuation accepts paired baseline initial cash with extra unpriced models", {
  yields <- make_ladder_yields()
  schedule <- make_ladder_schedule(yields, 2020:2021, q = c(0.01, 0.02))

  valuation <- tontine_price_baseline_models(
    baseline_curve_models = list(
      DNS = list(model = "DNS", curve_panels = yields)
    ),
    ladder_schedule = schedule,
    mat_names = c("X01Y", "X02Y"),
    initial_cash = c(Observed = 999, DNS = 1000)
  )

  expect_equal(valuation$initial_cash, c(DNS = 1000))
  expect_equal(unique(valuation$cash_accounts$model), "DNS")
})
