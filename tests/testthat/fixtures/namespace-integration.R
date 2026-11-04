local({
  setwd(tempdir())
  initial_search <- search()
  initial_globals <- ls(globalenv(), all.names = TRUE)
  initial_directory <- getwd()
  initial_environment <- Sys.getenv(c("LC_ALL", "LANGUAGE", "RSTUDIO_PANDOC"))

  library(CrossCurveR)

  stopifnot(
    identical(setdiff(search(), initial_search), "package:CrossCurveR"),
    identical(ls(globalenv(), all.names = TRUE), initial_globals),
    identical(getwd(), initial_directory),
    identical(Sys.getenv(c("LC_ALL", "LANGUAGE", "RSTUDIO_PANDOC")),
              initial_environment)
  )

  # Functions in a user's workspace must not replace package dependencies.
  for (name in c("bdiag", "transpose", "NS_loadings")) {
    assign(name, function(...) stop("Unexpected workspace dependency"), globalenv())
  }

  set.seed(321)
  mats <- c(1, 3, 6, 12, 24, 60)
  n <- 100L
  time <- as.Date("2020-01-01") + seq_len(n) - 1L
  betas <- cbind(2 + cumsum(stats::rnorm(n, sd = .03)),
                 -.8 + cumsum(stats::rnorm(n, sd = .01)),
                 .4 + stats::rnorm(n, sd = .02))
  yields <- stats::setNames(lapply(c(0, .2, .4), function(shift) {
    panel <- data.frame(time, betas %*% CrossCurveR::NS_loadings(.35, mats) +
                          shift + matrix(stats::rnorm(n * length(mats), sd = .01), n))
    colnames(panel) <- c("time", vapply(mats, CrossCurveR::numeric_to_matname,
                                       character(1)))
    panel
  }), c("usa", "gbr", "jpn"))

  dly <- CrossCurveR::dly_fit(yields, mats, lambdas = .35, reference = "usa")
  stopifnot(all(is.finite(as.matrix(dly$yields$dislocation_propagated$gbr[, -1]))))

  mce <- CrossCurveR::mc_fit_exo(
    yields, .35, c(2, 4), "usa", Xt = cbind(1, seq(-1, 1, length.out = n)),
    curves = names(yields), mats = mats, CR_algo = CrossCurveR::zero_mean_covreg_em,
    CR_maxiter = 3, CR_S0_shrink_diag = .4, CR_store_path = FALSE
  )
  stopifnot(all(is.finite(as.matrix(mce$curve[, -1]))), nrow(mce$curve) == n)

  shock <- CrossCurveR::st_apply_factor_shock_profile(
    curves = names(yields), curve_dns_factor = dly$country_betas$raw,
    shock_curves = "gbr", shock_factors = "L", shock_magnitudes = rep(0, n)
  )
  stopifnot(isTRUE(all.equal(shock$shocked_factors$gbr,
                            dly$country_betas$raw$gbr)))

  formatted <- CrossCurveR::dr_yield_format(
    data.frame(date = c("2020-01-01", "2020-01-02", "2020-01-03"),
               short = c(1, NaN, 3), long = c(2, 4, 6)),
    mats = c(1, 12), start_date = "2020-01-02", end_date = "2020-01-03"
  )
  stopifnot(identical(formatted$time, c("2020-01-02", "2020-01-03")),
            is.na(formatted[[2]][1]))

  cat("Namespace integration check passed\n")
})
