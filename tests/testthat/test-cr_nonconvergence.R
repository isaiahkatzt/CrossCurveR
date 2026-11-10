test_that("regularized covariance limits emit detailed nonconvergence warnings", {
  W <- matrix(c(-1, .2, .5, -.4, 1.2, .8, -.3, 1.1, .7, -.6), ncol = 2, byrow = TRUE)
  X <- cbind(1, seq(-1, 1, length.out = nrow(W)))
  estimators <- list(zero_mean_covreg_em = zero_mean_covreg_em,
                     covreg_comp_estimator = covreg_comp_estimator)
  for (name in names(estimators)) for (logging in c(FALSE, TRUE)) for (keep in c(FALSE, TRUE)) {
    result <- capture_covreg_warnings(estimators[[name]](W, X,
      B = matrix(c(.3, .2, -.1, .25), 2), S0 = diag(c(.4, .6)),
      max_iter = 1, tol = 1e-12, use_loglik = logging, store_path = keep))
    expect_length(result$warnings, 1)
    w <- result$warnings[[1]]
    expect_s3_class(w, "crosscurver_nonconvergence")
    expect_identical(w$estimator, name)
    expect_equal(w$iter, 1)
    expect_equal(w$max_iter, 1)
    expect_equal(w$tol, 1e-12)
    expect_equal(w$dimensions, c(N = 5, P = 2, Q = 2))
    expect_equal(w$last_check, 1)
    expect_false(result$value$converged)
    expect_match(conditionMessage(w), "max_iter=1 without confirmed convergence")
    expect_match(conditionMessage(w), "dB=")
    if (logging && keep) expect_true(is.finite(w$loglik_change)) else expect_null(w$loglik_change)
    if (logging) expect_true(is.finite(w$loglik)) else expect_null(w$loglik)
  }
})

test_that("convergence on the final allowed iteration does not warn", {
  W <- matrix(c(-1, .2, .5, -.4, 1.2, .8, -.3, 1.1, .7, -.6), ncol = 2, byrow = TRUE)
  X <- cbind(1, seq(-1, 1, length.out = nrow(W)))
  for (estimator in list(zero_mean_covreg_em, covreg_comp_estimator)) {
    expect_silent(fit <- estimator(W, X, B = matrix(0, 2, 2),
      S0 = Matrix::crossprod(W) / nrow(W), ridge_S0 = 0, max_iter = 1))
    expect_true(fit$converged)
    expect_equal(fit$iter, 1)
  }
  expect_silent(fit <- zero_mean_covreg(W, X, B = matrix(.01, 2, 2),
                                      max_iter = 10, tol = 1e8, term = TRUE))
  expect_equal(fit$iter, 10)
  expect_true(fit$converged)
})

test_that("unchecked iterations and nonfinite scores are identified", {
  W <- matrix(c(-1, .2, .5, -.4, 1.2, .8, -.3, 1.1, .7, -.6), ncol = 2, byrow = TRUE)
  X <- cbind(1, seq(-1, 1, length.out = nrow(W)))
  for (estimator in list(zero_mean_covreg_em, covreg_comp_estimator)) {
    result <- capture_covreg_warnings(estimator(W, X, B = matrix(.01, 2, 2),
      max_iter = 1, check_every = 2))
    expect_length(result$warnings, 1)
    expect_equal(result$warnings[[1]]$last_check, 0)
    expect_match(conditionMessage(result$warnings[[1]]), "no convergence check was performed")
    result <- capture_covreg_warnings(estimator(W, X, B = matrix(.01, 2, 2),
      max_iter = 3, check_every = 2, tol = 1e-100))
    expect_equal(result$warnings[[1]]$last_check, 2)
    expect_match(conditionMessage(result$warnings[[1]]), "between convergence checks")
    result <- with_temp_bindings(list(.cr_loglik = function(...) NA_real_),
      capture_covreg_warnings(estimator(W, X, B = matrix(.01, 2, 2), max_iter = 2)))
    expect_match(conditionMessage(result$warnings[[1]]), "likelihood score is nonfinite")
  }
})

test_that("legacy warnings distinguish history checks from fixed-iteration runs", {
  W <- matrix(c(-1, .2, .5, -.4, 1.2, .8, -.3, 1.1, .7, -.6), ncol = 2, byrow = TRUE)
  X <- cbind(1, seq(-1, 1, length.out = nrow(W)))
  result <- capture_covreg_warnings(zero_mean_covreg(W, X, B = matrix(.01, 2, 2),
    max_iter = 20, tol = 1e-100, term = TRUE))
  expect_length(result$warnings, 1)
  w <- result$warnings[[1]]
  expect_identical(w$estimator, "zero_mean_covreg")
  expect_equal(result$value$iter, w$iter)
  expect_false(result$value$converged)
  expect_equal(w$last_check, 20)
  expect_equal(w$parameter_changes, c(
    dB = Matrix::norm(result$value$b_est[[20]] - result$value$b_est[[11]], "F"),
    dS0 = Matrix::norm(result$value$s_est[[20]] - result$value$s_est[[11]], "F")))
  expect_match(conditionMessage(w), "legacy history comparison")
  result <- capture_covreg_warnings(zero_mean_covreg(W, X, B = matrix(.01, 2, 2),
    max_iter = 1, term = TRUE))
  expect_equal(result$warnings[[1]]$last_check, 0)
  expect_identical(result$value$termination_reason, "max_iter")
  expect_length(result$warnings[[1]]$parameter_changes, 0)
  result <- capture_covreg_warnings(zero_mean_covreg(W, X, B = matrix(.01, 2, 2),
    max_iter = 1, term = FALSE))
  expect_match(conditionMessage(result$warnings[[1]]), "not assessed \\(term=FALSE\\)")
  expect_identical(result$value$termination_reason, "not_assessed")
  expect_equal(result$warnings[[1]]$last_check, 0)
})

test_that("tenor wrappers add context and emit one warning per exhausted fit", {
  W <- cbind(usa = c(-1, .5, 1.2, -.3, .7), gbr = c(.2, -.4, .8, 1.1, -.6))
  X <- cbind(1, seq(-1, 1, length.out = nrow(W)))
  Wt <- list(X01M = W, X03M = W)
  for (named in c(TRUE, FALSE)) {
    result <- capture_covreg_warnings(blcc_fe(X, if (named) Wt else unname(Wt),
      mats = c(1, 3), B = matrix(.01, 2, 2), max_iter = 1))
    expect_length(result$warnings, 2)
    expect_identical(names(result$value), if (named) names(Wt) else NULL)
    expect_identical(vapply(result$warnings, `[[`, character(1), "tenor"),
                     if (named) c("X01M", "X03M") else c("01M", "03M"))
    expect_equal(vapply(result$warnings, `[[`, numeric(1), "maturity_months"), c(1, 3))
    expect_identical(result$warnings[[1]]$response_names, c("usa", "gbr"))
    expect_match(conditionMessage(result$warnings[[1]]), "Tenor .*01M.*responses \\[usa, gbr\\]")
  }
  backend <- function(W, X, ...) {
    warning("unrelated backend warning", call. = FALSE)
    list(S0 = diag(ncol(W)), B = matrix(0, ncol(W), ncol(X)))
  }
  expect_warning(blcc_fe(X, Wt[1], mats = 1, covreg = backend),
                 "^unrelated backend warning$")
})
