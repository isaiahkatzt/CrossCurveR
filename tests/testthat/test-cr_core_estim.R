test_that("zero_mean_covreg returns correctly shaped covariance-regression outputs", {
  set.seed(123)
  W <- cbind(seq(1, 10), seq(1, 10)^2)
  X <- cbind(1, seq(-1, 1, length.out = 10))

  fit <- quiet_covreg(zero_mean_covreg(W, X, init = "static", max_iter = 5))

  expect_true(is.matrix(fit$S0))
  expect_true(is.matrix(fit$B))
  expect_equal(dim(fit$S0), c(ncol(W), ncol(W)))
  expect_equal(dim(fit$B), c(ncol(W), ncol(X)))
})

test_that("zero_mean_covreg defaults to adaptive initialization", {
  W <- cbind(seq(1, 10), seq(1, 10)^2)
  X <- cbind(1, seq(-1, 1, length.out = 10))

  set.seed(123)
  default_fit <- quiet_covreg(zero_mean_covreg(W, X, max_iter = 5))
  set.seed(123)
  adaptive_fit <- quiet_covreg(zero_mean_covreg(W, X, init = "adaptive", max_iter = 5))

  expect_identical(default_fit, adaptive_fit)
})

test_that("zero_mean_covreg supports adaptive initialization and early termination mode", {
  set.seed(123)
  W <- cbind(seq(1, 12), seq(1, 12)^2)
  X <- cbind(1, seq(-1, 1, length.out = 12))

  fit <- zero_mean_covreg(W, X, init = "adaptive", max_iter = 20, tol = 1e8, term = TRUE)

  expect_true(is.list(fit$b_est))
  expect_true(is.list(fit$s_est))
  expect_lte(fit$iter, 20)
})

test_that("legacy covariance outputs share core fields across termination paths", {
  W <- matrix(c(-1, .2, .5, -.4, 1.2, .8, -.3, 1.1, .7, -.6), ncol = 2, byrow = TRUE)
  X <- cbind(1, seq(-1, 1, length.out = nrow(W)))
  B <- matrix(.01, 2, 2)
  core <- quiet_covreg(zero_mean_covreg_em(W, X, B = B, max_iter = 1, use_loglik = FALSE))
  for (term in c(FALSE, TRUE)) for (limit in c(1, 20)) for (tolerance in c(1e-100, 1e8)) {
    fit <- quiet_covreg(zero_mean_covreg(W, X, B = B, max_iter = limit,
                                        tol = tolerance, term = term))
    converged <- term && limit >= 10 && tolerance == 1e8
    expect_true(all(names(core) %in% names(fit)))
    expect_identical(fit$converged, converged)
    expect_equal(fit$iter, if (converged) 10 else limit)
    expect_identical(fit$termination_reason,
                     if (converged) "converged" else if (term) "max_iter" else "not_assessed")
    expect_null(fit$loglik)
    expect_identical(fit$B_path, fit$b_est[seq_len(fit$iter + 1)])
    expect_identical(fit$S0_path, fit$s_est[seq_len(fit$iter + 1)])
    expect_identical(fit$B_path[[fit$iter + 1]], fit$B)
    expect_identical(fit$S0_path[[fit$iter + 1]], fit$S0)
    expect_true(all(is.finite(c(fit$mz, fit$vz, fit$sz))))
    expect_length(fit$mz, nrow(W))
    expect_equal(fit$sz^2, fit$vz)
  }
})

test_that("covariance regression shrinkage can diagonalize S0", {
  set.seed(123)
  W <- cbind(seq(1, 10), seq(1, 10) + c(rep(0, 5), rep(3, 5)))
  X <- cbind(1, seq(-1, 1, length.out = 10))

  fit_dense <- quiet_covreg(zero_mean_covreg_em(W, X, init = "static", max_iter = 3, S0_shrink_diag = 0))
  fit_diag <- quiet_covreg(zero_mean_covreg_em(W, X, init = "static", max_iter = 3, S0_shrink_diag = 1))

  expect_gt(abs(fit_dense$S0[1, 2]), 0)
  expect_equal(fit_diag$S0[1, 2], 0, tolerance = 1e-10)
  expect_equal(fit_diag$S0_shrink_diag, 1)
})

test_that("residual sign reversal preserves covariances with later convergence checks", {
  W <- matrix(c(-1, .2, .5, -.4, 1.2, .8, -.3, 1.1, .7, -.6),
              ncol = 2, byrow = TRUE)
  X <- cbind(1, seq(-1, 1, length.out = nrow(W)))

  for (estimator in list(zero_mean_covreg_em, zero_mean_covreg, covreg_comp_estimator)) {
    for (init in c("adaptive", "static")) {
      args <- list(X = X, init = init, max_iter = 20, tol = 1e-8,
                   S0_shrink_diag = 0.4)
      legacy <- identical(estimator, zero_mean_covreg)
      # The sign-normalized EM initialization can affect the first stopping check.
      # Later checks isolate the covariance and coefficient-orientation contracts.
      if (legacy) args$term <- TRUE else args$check_every <- 2

      set.seed(451)
      original <- quiet_covreg(do.call(estimator, c(list(W = W), args)))
      original_rng <- .Random.seed
      set.seed(451)
      reversed <- quiet_covreg(do.call(estimator, c(list(W = -W), args)))

      expect_equal(reversed$S0, original$S0, tolerance = 1e-12)
      expect_equal(cr_build(reversed, X), cr_build(original, X), tolerance = 1e-12)
      expect_identical(.Random.seed, original_rng)
      expect_identical(reversed$iter, original$iter)
      expected_B <- if (legacy && init == "adaptive") -original$B else original$B
      expect_equal(reversed$B, expected_B, tolerance = 1e-12)
    }
  }
})

test_that("regularized estimators preserve the historical inversion ridge and score", {
  W <- matrix(c(-1, .2, .5, -.4, 1.2, .8, -.3, 1.1, .7, -.6), ncol = 2, byrow = TRUE)
  X <- cbind(1, seq(-1, 1, length.out = nrow(W)))
  B0 <- matrix(c(.3, .2, -.1, .25), nrow = 2)
  S0 <- diag(c(.4, .6))
  ridge <- .2
  initial_covariance <- S0 + diag(2 * ridge, 2)
  u <- X %*% t(B0)
  precision <- solve(initial_covariance)
  expected_vz <- 1 / (1 + rowSums((u %*% precision) * u))
  expected_mz <- expected_vz * rowSums((u %*% precision) * W)

  for (estimator in list(zero_mean_covreg_em, covreg_comp_estimator)) {
    fit <- quiet_covreg(estimator(W, X, B = B0, S0 = S0, ridge_S0 = ridge, max_iter = 1))
    expect_equal(fit$vz, expected_vz, tolerance = 1e-12)
    expect_equal(fit$mz, expected_mz, tolerance = 1e-12)

    # Evaluate the historical score using full matrix solves and determinants.
    # It combines the stored-covariance determinant with a ridged precision.
    inverse_covariance <- fit$S0 + diag(ridge, 2)
    logdet_offset <- as.numeric(determinant(fit$S0, logarithm = TRUE)$modulus) -
      as.numeric(determinant(inverse_covariance, logarithm = TRUE)$modulus)
    direct_loglik <- sum(vapply(seq_len(nrow(W)), function(i) {
      covariance <- inverse_covariance + tcrossprod(fit$B %*% X[i, ])
      -.5 * (ncol(W) * log(2 * pi) +
        as.numeric(determinant(covariance, logarithm = TRUE)$modulus) + logdet_offset +
        as.numeric(crossprod(W[i, ], solve(covariance, W[i, ]))))
    }, numeric(1)))
    expect_equal(fit$loglik, direct_loglik, tolerance = 1e-12)
  }
})

test_that("fixed-iteration EM updates match the pre-correction reference", {
  # Reference values were generated from aee15d9^ before the ridge correction.
  # Compare the same number of updates; convergence thresholds are tested separately.
  W <- matrix(c(-1, .2, .5, -.4, 1.2, .8, -.3, 1.1, .7, -.6), ncol = 2, byrow = TRUE)
  X <- cbind(1, seq(-1, 1, length.out = nrow(W)))
  expected <- list(zero_mean_covreg_em = list(adaptive = list(B = structure(c(-0.0854461763877651,
-0.018112690300922, 0.109287550150037, -0.0372771409652622), dim = c(2L,
2L)), S0 = structure(c(0.84083942244831, -0.0375166751094538,
-0.0375166751094538, 0.680978650390849), dim = c(2L, 2L)), loglik_path = c(-10.6829353454911,
-10.6879031362601, -10.7113087192477, -10.7237466859472, -10.730516790528
), iter = 4, converged = FALSE), static = list(B = structure(c(-7.38210084477577e-05,
-3.91668656321759e-06, 0.000109876306859027, -3.32797632593783e-05
), dim = c(2L, 2L)), S0 = structure(c(0.853999988514057, -0.037999998460805,
-0.037999998460805, 0.681999999430889), dim = c(2L, 2L)), loglik_path = c(-10.6032124398426,
-10.7460118166549, -10.7460118270088, -10.7460118324504), iter = 3,
    converged = TRUE)), covreg_comp_estimator = list(adaptive = list(
    B = structure(c(-0.149977929981612, -0.0153849796912038,
    0.170467461016641, -0.0455513552116462), dim = c(2L, 2L)),
    S0 = structure(c(0.6, 0, 0, 0.8), dim = c(2L, 2L)), loglik_path = c(-10.6829353454911,
    -10.6266475540998, -10.6136648757209, -10.6075900726767,
    -10.6041675193585), iter = 4, converged = FALSE), static = list(
    B = structure(c(-0.000128572037750332, 9.26622999894188e-06,
    0.000182120856736871, -4.34163988228258e-05), dim = c(2L,
    2L)), S0 = structure(c(0.6, 0, 0, 0.8), dim = c(2L, 2L)),
    loglik_path = c(-10.6032124398426, -10.6032123954475, -10.6032123901774
    ), iter = 2, converged = TRUE)))
  estimators <- list(zero_mean_covreg_em = zero_mean_covreg_em,
                     covreg_comp_estimator = covreg_comp_estimator)
  for (name in names(estimators)) {
    for (init in c("adaptive", "static")) {
      reference <- expected[[name]][[init]]
      set.seed(451)
      fit <- quiet_covreg(estimators[[name]](W, X, init = init, S0 = diag(c(.4, .6)),
                               ridge_S0 = .2, max_iter = reference$iter,
                               check_every = reference$iter + 1))
      reference$converged <- NULL
      expect_equal(fit[names(reference)], reference, tolerance = 1e-12)
      expect_false(fit$converged)
    }
  }
})

test_that("stored EM paths require tol-sized parameter changes", {
  W <- matrix(c(-1, .2, .5, -.4, 1.2, .8, -.3, 1.1, .7, -.6),
              ncol = 2, byrow = TRUE)
  X <- cbind(1, seq(-1, 1, length.out = nrow(W)))
  tolerance <- 1e-8

  for (estimator in list(zero_mean_covreg_em, covreg_comp_estimator)) {
    args <- list(W = W, X = X, init = "static", S0 = diag(c(.4, .6)),
                 ridge_S0 = .2, max_iter = 1000, tol = tolerance, use_loglik = TRUE)
    set.seed(451)
    stored <- do.call(estimator, c(args, list(store_path = TRUE)))
    stored_rng <- .Random.seed
    set.seed(451)
    unstored <- do.call(estimator, c(args, list(store_path = FALSE)))

    expect_true(stored$converged)
    expect_true(unstored$converged)
    expect_identical(stored$iter, unstored$iter)
    expect_equal(stored$B, unstored$B, tolerance = 1e-12)
    expect_equal(stored$S0, unstored$S0, tolerance = 1e-12)
    expect_identical(.Random.seed, stored_rng)
    last <- length(stored$B_path)
    expect_lt(Matrix::norm(stored$B_path[[last]] - stored$B_path[[last - 1]], "F"), tolerance)
    expect_lt(Matrix::norm(stored$S0_path[[last]] - stored$S0_path[[last - 1]], "F"), tolerance)
    expect_lt(abs(stored$loglik_path[last] - stored$loglik_path[last - 1]), tolerance)
  }
})

test_that("covariance sample errors identify the requirement and received count", {
  for (estimator in list(zero_mean_covreg_em, covreg_comp_estimator)) {
    expect_error(estimator(matrix(c(1, 2), nrow = 1), matrix(1, 1, 1)),
                 "Covariance regression: at least two observations are required; received 1")
  }
})

test_that("regularized estimators reject invalid iteration controls before fitting", {
  invalid <- list(
    max_iter = list(0, -1, 1.5, Inf, NA_real_, numeric(0), c(1, 2), "1"),
    check_every = list(0, -1, 1.5, Inf, NA_real_),
    tol = list(0, -1, Inf, NA_real_)
  )
  for (estimator in list(zero_mean_covreg_em, covreg_comp_estimator)) {
    for (control in names(invalid)) {
      for (value in invalid[[control]]) {
        args <- c(list(W = NULL, X = NULL), stats::setNames(list(value), control))
        expect_error(do.call(estimator, args), paste0("`", control, "` must be"))
      }
    }
  }
})

test_that("one iteration remains valid without an intervening convergence check", {
  W <- matrix(c(-1, .2, .5, -.4, 1.2, .8, -.3, 1.1, .7, -.6), ncol = 2, byrow = TRUE)
  X <- cbind(1, seq(-1, 1, length.out = nrow(W)))
  for (estimator in list(zero_mean_covreg_em, covreg_comp_estimator)) {
    fit <- quiet_covreg(estimator(W, X, B = matrix(0, 2, 2), S0 = diag(2), max_iter = 1L,
                     check_every = 2L, tol = 1e-8, store_path = FALSE))
    expect_equal(fit$iter, 1)
    expect_false(fit$converged)
    expect_true(all(is.finite(c(fit$mz, fit$vz, fit$loglik))))
  }
})
