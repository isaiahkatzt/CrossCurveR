#############################################
#           Ogden Application Helpers       #
#############################################

ogden_build_liability_cashflows <- function(M_0, q, b0, g_b, payment_years) {
  payment_index <- seq_along(q)
  survival_to_start <- c(1, cumprod(1 - q[-length(q)]))
  expected_claimants <- round(M_0 * survival_to_start)
  expected_deaths <- expected_claimants * q
  per_claimant_benefit <- b0 * (1 + g_b) ^ payment_index

  tibble::tibble(
    payment_index = payment_index,
    payment_year = payment_years,
    q = q,
    survival_to_start = survival_to_start,
    expected_claimants = expected_claimants,
    expected_deaths = expected_deaths,
    per_claimant_benefit = per_claimant_benefit,
    liability_cashflow = expected_claimants * per_claimant_benefit
  )
}

ogden_curve_rate <- function(curve_panel, target_date, horizon_years,
                             mats, mat_names) {
  interpolate_curve_rate(
    curve_panel = curve_panel,
    target_date = target_date,
    target_months = 12 * horizon_years,
    mats = mats,
    mat_names = mat_names,
    rate_scale = "decimal"
  )
}

ogden_common_curve_dates <- function(curve_panels, start_date, end_date) {
  curve_dates <- lapply(curve_panels, function(curve_panel) {
    curve_time <- as.Date(curve_panel$time)
    as.character(curve_time[curve_time >= start_date & curve_time <= end_date])
  })

  as.Date(sort(Reduce(intersect, curve_dates)))
}

ogden_box_constrained_gmv_weights <- function(cov_mat, weight_floor = 0,
                                              weight_ceiling = 1,
                                              tol = 1e-10) {
  cov_mat <- as.matrix(cov_mat)
  assets <- colnames(cov_mat)
  n_assets <- ncol(cov_mat)

  if (n_assets * weight_floor > 1 + tol || n_assets * weight_ceiling < 1 - tol) {
    stop("Infeasible GMV box constraints.")
  }

  bound_states <- c("free", "lower", "upper")
  active_sets <- expand.grid(
    replicate(n_assets, bound_states, simplify = FALSE),
    stringsAsFactors = FALSE
  )

  best_var <- Inf
  best_weights <- rep(NA_real_, n_assets)

  for (row_idx in seq_len(nrow(active_sets))) {
    state <- as.character(active_sets[row_idx, ])
    fixed_idx <- which(state != "free")
    free_idx <- which(state == "free")
    weights <- rep(NA_real_, n_assets)

    if (length(fixed_idx) > 0) {
      weights[fixed_idx] <- ifelse(
        state[fixed_idx] == "lower",
        weight_floor,
        weight_ceiling
      )
    }

    remaining_weight <- 1 - sum(weights[fixed_idx], na.rm = TRUE)
    if (remaining_weight < -tol) {
      next
    }

    if (length(free_idx) == 0) {
      if (abs(remaining_weight) > sqrt(tol)) {
        next
      }
    } else {
      if (remaining_weight < length(free_idx) * weight_floor - tol ||
          remaining_weight > length(free_idx) * weight_ceiling + tol) {
        next
      }

      sigma_ff <- cov_mat[free_idx, free_idx, drop = FALSE]
      sigma_fb <- if (length(fixed_idx) > 0) {
        cov_mat[free_idx, fixed_idx, drop = FALSE]
      } else {
        matrix(0, nrow = length(free_idx), ncol = 1)
      }
      fixed_weights <- if (length(fixed_idx) > 0) weights[fixed_idx] else 0
      rhs_shift <- if (length(fixed_idx) > 0) {
        as.numeric(sigma_fb %*% fixed_weights)
      } else {
        rep(0, length(free_idx))
      }

      unit <- rep(1, length(free_idx))
      kkt <- rbind(
        cbind(2 * sigma_ff + diag(tol, length(free_idx)), unit),
        c(unit, 0)
      )
      rhs <- c(-2 * rhs_shift, remaining_weight)
      solution <- tryCatch(solve(kkt, rhs), error = function(e) NULL)
      if (is.null(solution)) {
        next
      }

      weights[free_idx] <- solution[seq_along(free_idx)]
    }

    if (any(weights < weight_floor - sqrt(tol)) ||
        any(weights > weight_ceiling + sqrt(tol)) ||
        abs(sum(weights) - 1) > sqrt(tol)) {
      next
    }

    weights <- pmin(pmax(weights, weight_floor), weight_ceiling)
    weights <- weights / sum(weights)
    if (any(weights < weight_floor - sqrt(tol)) ||
        any(weights > weight_ceiling + sqrt(tol))) {
      next
    }

    portfolio_var <- as.numeric(t(weights) %*% cov_mat %*% weights)
    if (portfolio_var < best_var) {
      best_var <- portfolio_var
      best_weights <- weights
    }
  }

  if (!all(is.finite(best_weights))) {
    stop("No feasible GMV box-constrained weights found.")
  }

  names(best_weights) <- assets
  as.numeric(best_weights)
}

ogden_estimate_gmv_weights <- function(yields, curves, horizons,
                                       start_date, end_date,
                                       mats, mat_names,
                                       weight_floor,
                                       weight_ceiling) {
  calibration_dates <- ogden_common_curve_dates(
    curve_panels = yields[curves],
    start_date = start_date,
    end_date = end_date
  )

  purrr::map_dfr(horizons, function(horizon_years) {
    log_price_panel <- do.call(cbind, lapply(curves, function(curve_name) {
      zero_rates <- vapply(calibration_dates, function(date_value) {
        ogden_curve_rate(
          curve_panel = yields[[curve_name]],
          target_date = date_value,
          horizon_years = horizon_years,
          mats = mats,
          mat_names = mat_names
        )
      }, numeric(1))
      -zero_rates * horizon_years
    }))
    colnames(log_price_panel) <- curves

    weights <- ogden_box_constrained_gmv_weights(
      cov_mat = cov(diff(log_price_panel)),
      weight_floor = weight_floor,
      weight_ceiling = weight_ceiling
    )
    names(weights) <- curves

    tibble::tibble(
      horizon_years = horizon_years,
      country = names(weights),
      portfolio_weight = as.numeric(weights)
    )
  })
}

ogden_build_cashflow_matched_assets <- function(liability_cashflows,
                                                portfolio_weights,
                                                curve_panels,
                                                purchase_date,
                                                mats,
                                                mat_names) {
  portfolio_purchase_date <- as.Date(purchase_date)

  liability_cashflows %>%
    dplyr::select(payment_index, payment_year, liability_cashflow) %>%
    dplyr::inner_join(
      portfolio_weights,
      by = c("payment_index" = "horizon_years")
    ) %>%
    dplyr::arrange(payment_index, country) %>%
    dplyr::mutate(
      purchase_date = portfolio_purchase_date,
      horizon_years = payment_index,
      asset_id = paste0(toupper(country), "_", payment_year),
      country_par = liability_cashflow * portfolio_weight,
      purchase_zero_rate = purrr::map2_dbl(
        country,
        horizon_years,
        function(country_name, maturity_years) {
          ogden_curve_rate(
            curve_panel = curve_panels[[country_name]],
            target_date = portfolio_purchase_date,
            horizon_years = maturity_years,
            mats = mats,
            mat_names = mat_names
          )
        }
      ),
      purchase_discount_factor = exp(-purchase_zero_rate * horizon_years),
      purchase_cost = country_par * purchase_discount_factor
    ) %>%
    dplyr::select(
      asset_id,
      country,
      payment_index,
      payment_year,
      purchase_date,
      horizon_years,
      liability_cashflow,
      portfolio_weight,
      country_par,
      purchase_zero_rate,
      purchase_discount_factor,
      purchase_cost
    )
}

ogden_summarise_asset_cashflow_match <- function(asset_portfolio) {
  asset_portfolio %>%
    dplyr::group_by(payment_index, payment_year, horizon_years) %>%
    dplyr::summarise(
      liability_cashflow = dplyr::first(liability_cashflow),
      asset_cashflow = sum(country_par),
      cashflow_difference = asset_cashflow - liability_cashflow,
      purchase_cost = sum(purchase_cost),
      .groups = "drop"
    )
}

ogden_summarise_asset_purchase <- function(asset_portfolio) {
  asset_portfolio %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(
      par_value = sum(country_par),
      purchase_cost = sum(purchase_cost),
      .groups = "drop"
    )
}

ogden_build_curve_panel_sets <- function(observed_yields, baseline_curve_models,
                                         model_order = c("Observed", "DNS", "DLY")) {
  model_sets <- list(Observed = observed_yields)
  model_sets <- c(
    model_sets,
    lapply(baseline_curve_models, function(model_obj) model_obj$curve_panels)
  )

  model_sets[intersect(model_order, names(model_sets))]
}

ogden_value_asset_portfolio <- function(asset_portfolio,
                                        valuation_windows,
                                        curve_panels,
                                        curves,
                                        mats,
                                        mat_names) {
  purrr::map_dfr(seq_len(nrow(valuation_windows)), function(i) {
    valuation_row <- valuation_windows[i, , drop = FALSE]
    valuation_index <- valuation_row$rung_id[[1]]
    valuation_date_value <- as.Date(valuation_row$purchase_date[[1]])

    outstanding <- asset_portfolio[
      asset_portfolio$payment_index >= valuation_index,
      ,
      drop = FALSE
    ]

    outstanding %>%
      dplyr::mutate(
        valuation_index = valuation_index,
        valuation_year = valuation_row$purchase_year[[1]],
        valuation_date = valuation_date_value,
        remaining_horizon_years = payment_index - valuation_index,
        valuation_zero_rate = purrr::map2_dbl(
          country,
          remaining_horizon_years,
          function(country_name, horizon_years) {
            if (horizon_years == 0) {
              return(NA_real_)
            }
            if (!country_name %in% curves) {
              stop(sprintf("Asset country `%s` is not included in `curves`.", country_name))
            }
            ogden_curve_rate(
              curve_panel = curve_panels[[country_name]],
              target_date = valuation_date_value,
              horizon_years = horizon_years,
              mats = mats,
              mat_names = mat_names
            )
          }
        ),
        valuation_discount_factor = dplyr::if_else(
          remaining_horizon_years == 0,
          1,
          exp(-valuation_zero_rate * remaining_horizon_years)
        ),
        asset_present_value = country_par * valuation_discount_factor
      ) %>%
      dplyr::select(
        valuation_index,
        valuation_year,
        valuation_date,
        asset_id,
        country,
        payment_index,
        payment_year,
        remaining_horizon_years,
        country_par,
        valuation_zero_rate,
        valuation_discount_factor,
        asset_present_value
      )
  })
}

ogden_value_liability_cashflows <- function(liability_cashflows,
                                            valuation_windows,
                                            portfolio_weights,
                                            curve_panels,
                                            curves,
                                            mats,
                                            mat_names,
                                            psi_1 = 0,
                                            psi_2 = 0,
                                            psi_max_horizon = max(liability_cashflows$payment_index),
                                            psi_panel = NULL) {
  purrr::map_dfr(seq_len(nrow(valuation_windows)), function(i) {
    valuation_row <- valuation_windows[i, , drop = FALSE]
    valuation_index <- valuation_row$rung_id[[1]]
    valuation_date_value <- as.Date(valuation_row$purchase_date[[1]])
    valuation_year_value <- valuation_row$purchase_year[[1]]

    outstanding <- liability_cashflows[
      liability_cashflows$payment_index >= valuation_index,
      ,
      drop = FALSE
    ]

    outstanding %>%
      dplyr::mutate(
        valuation_index = valuation_index,
        valuation_year = valuation_year_value,
        valuation_date = valuation_date_value,
        remaining_horizon_years = payment_index - valuation_index,
        liability_market_discount_rate = purrr::map_dbl(
          remaining_horizon_years,
          function(horizon_years) {
            if (horizon_years == 0) {
              return(NA_real_)
            }
            weight_rows <- portfolio_weights[
              portfolio_weights$horizon_years == horizon_years,
              ,
              drop = FALSE
            ]
            weight_vector <- setNames(weight_rows$portfolio_weight, weight_rows$country)
            country_rates <- vapply(curves, function(curve_name) {
              ogden_curve_rate(
                curve_panel = curve_panels[[curve_name]],
                target_date = valuation_date_value,
                horizon_years = horizon_years,
                mats = mats,
                mat_names = mat_names
              )
            }, numeric(1))
            sum(weight_vector[curves] * country_rates[curves])
          }
        )
      ) %>%
      {
        if (is.null(psi_panel)) {
          dplyr::mutate(
            .,
            psi_horizon_loading = remaining_horizon_years / psi_max_horizon,
            psi_1_jt = psi_1 * psi_horizon_loading,
            psi_2_jt = psi_2 * psi_horizon_loading
          )
        } else {
          psi_rows <- psi_panel %>%
            dplyr::filter(valuation_year == valuation_year_value) %>%
            dplyr::select(
              remaining_horizon_years,
              psi_horizon_loading,
              psi_1_jt,
              psi_2_jt
            )

          dplyr::left_join(., psi_rows, by = "remaining_horizon_years") %>%
            dplyr::mutate(
              psi_horizon_loading = dplyr::coalesce(
                psi_horizon_loading,
                remaining_horizon_years / psi_max_horizon
              ),
              psi_1_jt = dplyr::coalesce(psi_1_jt, 0),
              psi_2_jt = dplyr::coalesce(psi_2_jt, 0)
            )
        }
      } %>%
      dplyr::mutate(
        psi_jt = psi_1_jt + psi_2_jt,
        liability_discount_rate = liability_market_discount_rate - psi_jt,
        liability_discount_factor = dplyr::if_else(
          remaining_horizon_years == 0,
          1,
          exp(-liability_discount_rate * remaining_horizon_years)
        ),
        liability_present_value = liability_cashflow * liability_discount_factor
      ) %>%
      dplyr::select(
        valuation_index,
        valuation_year,
        valuation_date,
        payment_index,
        payment_year,
        remaining_horizon_years,
        liability_cashflow,
        liability_market_discount_rate,
        psi_horizon_loading,
        psi_1_jt,
        psi_2_jt,
        psi_jt,
        liability_discount_rate,
        liability_discount_factor,
        liability_present_value
      )
  })
}

ogden_summarise_asset_value <- function(asset_value_panel) {
  asset_value_panel %>%
    dplyr::group_by(valuation_index, valuation_year, valuation_date) %>%
    dplyr::summarise(
      asset_cashflow_due = sum(country_par[remaining_horizon_years == 0]),
      remaining_asset_cashflow = sum(country_par),
      asset_present_value = sum(asset_present_value),
      .groups = "drop"
    ) %>%
    dplyr::arrange(valuation_index)
}

ogden_summarise_liability_value <- function(liability_value_panel) {
  liability_value_panel %>%
    dplyr::group_by(valuation_index, valuation_year, valuation_date) %>%
    dplyr::summarise(
      liability_cashflow_due = sum(liability_cashflow[remaining_horizon_years == 0]),
      remaining_liability_cashflow = sum(liability_cashflow),
      liability_present_value = sum(liability_present_value),
      .groups = "drop"
    ) %>%
    dplyr::arrange(valuation_index)
}

ogden_value_assets_by_model <- function(asset_portfolio,
                                        valuation_windows,
                                        curve_panel_sets,
                                        curves,
                                        mats,
                                        mat_names) {
  purrr::imap_dfr(curve_panel_sets, function(curve_panels, model_name) {
    ogden_summarise_asset_value(
      ogden_value_asset_portfolio(
        asset_portfolio = asset_portfolio,
        valuation_windows = valuation_windows,
        curve_panels = curve_panels,
        curves = curves,
        mats = mats,
        mat_names = mat_names
      )
    ) %>%
      dplyr::mutate(model = model_name, .before = 1)
  })
}

ogden_value_liabilities_by_model <- function(liability_cashflows,
                                             valuation_windows,
                                             portfolio_weights,
                                             curve_panel_sets,
                                             curves,
                                             mats,
                                             mat_names,
                                             psi_1 = 0,
                                             psi_2 = 0,
                                             psi_max_horizon = max(liability_cashflows$payment_index),
                                             psi_panel = NULL) {
  purrr::imap_dfr(curve_panel_sets, function(curve_panels, model_name) {
    model_psi_panel <- if (!is.null(psi_panel) && "model" %in% colnames(psi_panel)) {
      psi_panel %>% dplyr::filter(model == model_name)
    } else {
      psi_panel
    }

    ogden_summarise_liability_value(
      ogden_value_liability_cashflows(
        liability_cashflows = liability_cashflows,
        valuation_windows = valuation_windows,
        portfolio_weights = portfolio_weights,
        curve_panels = curve_panels,
        curves = curves,
        mats = mats,
        mat_names = mat_names,
        psi_1 = psi_1,
        psi_2 = psi_2,
        psi_max_horizon = psi_max_horizon,
        psi_panel = model_psi_panel
      )
    ) %>%
      dplyr::mutate(model = model_name, .before = 1)
  })
}

ogden_summarise_model_valuations <- function(asset_valuation_summary,
                                             liability_valuation_summary) {
  asset_valuation_summary %>%
    dplyr::inner_join(
      liability_valuation_summary,
      by = c("model", "valuation_index", "valuation_year", "valuation_date")
    ) %>%
    dplyr::mutate(
      surplus_loss = asset_present_value - liability_present_value,
      funding_ratio = dplyr::if_else(
        liability_present_value == 0,
        NA_real_,
        asset_present_value / liability_present_value
      ),
      current_cashflow_difference = asset_cashflow_due - liability_cashflow_due
    ) %>%
    dplyr::arrange(model, valuation_index)
}

ogden_interpolate_dislocation <- function(dislocation_panel, target_horizon,
                                          country_name) {
  country_rows <- dislocation_panel[
    dislocation_panel$country == country_name,
    ,
    drop = FALSE
  ]

  if (target_horizon == 0) {
    return(0)
  }

  stats::approx(
    x = country_rows$tenor_years,
    y = country_rows$residual_decimal,
    xout = target_horizon,
    rule = 2,
    ties = "ordered"
  )$y
}

ogden_extract_mce_residual_transmission <- function(mce_model,
                                                    valuation_windows,
                                                    curves,
                                                    mats,
                                                    mat_names) {
  target_years <- valuation_windows$purchase_year
  model_windows <- mce_model$valuation_windows
  model_indices <- which(model_windows$purchase_year %in% target_years)

  purrr::map_dfr(model_indices, function(window_index) {
    window_fit <- mce_model$fit[[window_index]]$fit
    window_row <- model_windows[window_index, , drop = FALSE]
    final_x <- as.numeric(window_fit$xt_state$Xt[nrow(window_fit$xt_state$Xt), ])

    tenor_dislocations <- purrr::map2_dfr(
      window_fit$BS0,
      seq_along(window_fit$BS0),
      function(BS0_j, tenor_index) {
        B_j <- as.matrix(BS0_j$B)
        row_names <- rownames(B_j)
        if (is.null(row_names)) {
          row_names <- curves
        }
        residuals <- as.numeric(B_j %*% final_x)

        tibble::tibble(
          tenor_name = mat_names[[tenor_index]],
          tenor_years = mats[[tenor_index]] / 12,
          country = row_names,
          residual_percent = residuals,
          residual_decimal = residuals / 100
        )
      }
    )

    tenor_dislocations %>%
      dplyr::mutate(
        valuation_index = window_row$rung_id[[1]],
        valuation_year = window_row$purchase_year[[1]],
        valuation_date = as.Date(window_row$purchase_date[[1]]),
        .before = 1
      )
  }) %>%
    dplyr::arrange(valuation_year, tenor_years, country)
}

ogden_shock_psi <- function(baseline_mce_model,
                            shocked_mce_model,
                            valuation_windows,
                            curves,
                            mats,
                            mat_names,
                            reference_curve,
                            liability_horizons,
                            dislocation_scale = 1,
                            cross_country_scale = 1) {
  shock_transmission <- ogden_extract_mce_residual_transmission(
    mce_model = shocked_mce_model,
    valuation_windows = valuation_windows,
    curves = curves,
    mats = mats,
    mat_names = mat_names
  ) %>%
    dplyr::rename(
      shocked_residual_percent = residual_percent,
      shocked_residual_decimal = residual_decimal
    )

  purrr::map_dfr(unique(shock_transmission$valuation_year), function(year_value) {
    year_row <- shock_transmission %>%
      dplyr::filter(valuation_year == year_value)
    valuation_row <- year_row[1, , drop = FALSE]

    purrr::map_dfr(liability_horizons, function(horizon_years) {
      psi_horizon_loading <- horizon_years / max(liability_horizons)
      reference_residual <- ogden_interpolate_dislocation(
        dislocation_panel = year_row %>%
          dplyr::select(country, tenor_years, residual_decimal = shocked_residual_decimal),
        target_horizon = horizon_years,
        country_name = reference_curve
      )
      cross_residual <- sum(vapply(
        setdiff(curves, reference_curve),
        function(curve_name) {
          ogden_interpolate_dislocation(
            dislocation_panel = year_row %>%
              dplyr::select(country, tenor_years, residual_decimal = shocked_residual_decimal),
            target_horizon = horizon_years,
            country_name = curve_name
          )
        },
        numeric(1)
      ))

      tibble::tibble(
        model = "MCE",
        valuation_index = valuation_row$valuation_index[[1]],
        valuation_year = valuation_row$valuation_year[[1]],
        valuation_date = valuation_row$valuation_date[[1]],
        remaining_horizon_years = horizon_years,
        psi_horizon_loading = psi_horizon_loading,
        reference_residual_decimal = reference_residual,
        cross_residual_decimal = cross_residual,
        psi_1_jt = dislocation_scale * psi_horizon_loading * reference_residual,
        psi_2_jt = dislocation_scale * cross_country_scale *
          psi_horizon_loading * cross_residual,
        psi_jt = psi_1_jt + psi_2_jt
      )
    })
  }) %>%
    dplyr::arrange(valuation_year, remaining_horizon_years)
}

ogden_hybrid_shock_psi <- function(mce_psi_panel,
                                   exogenous_psi_1,
                                   exogenous_psi_2 = 0,
                                   psi_max_horizon = max(mce_psi_panel$remaining_horizon_years)) {
  mce_psi_panel %>%
    dplyr::mutate(
      exogenous_psi_1_jt = exogenous_psi_1 *
        remaining_horizon_years / psi_max_horizon,
      exogenous_psi_2_jt = exogenous_psi_2 *
        remaining_horizon_years / psi_max_horizon,
      mce_psi_1_jt = psi_1_jt,
      mce_psi_2_jt = psi_2_jt,
      psi_1_jt = exogenous_psi_1_jt + mce_psi_1_jt,
      psi_2_jt = exogenous_psi_2_jt + mce_psi_2_jt,
      psi_jt = psi_1_jt + psi_2_jt
    )
}
