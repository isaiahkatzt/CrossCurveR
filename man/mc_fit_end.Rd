\name{mc_fit_end}
\alias{mc_fit_end}
\alias{mc_fit_exo}
\title{Fit multi-curve endogenous or exogenous yield models}
\description{
Calibrate independent DNS baselines, estimate maturity-specific
conditional covariance, and reconstruct covariance-transformed curves.
}
\usage{
mc_fit_end(yields, lambdas, cutoffs, reference, ECM_estim = "ML",
    ECM_type = "eigen", ECM_alpha = 0.1, X_normalize = TRUE,
    X_trunc = FALSE, CR_algo = zero_mean_covreg_em, CR_init = "adaptive",
    Xt_smooth = FALSE, smoother = c("ns", "bs", "rm", "hp", "henderson"),
    knot_count = 5, k_count = 9, k_pass = 3, hp_lambda = 1600,
    henderson_k = 13, CR_maxiter = 1000, CR_tol = 1e-08, CR_Binit = NULL,
    CR_S0init = NULL, CR_S0_shrink_diag = 0, CR_check_every = 1,
    CR_use_loglik = TRUE, CR_store_path = TRUE, CR_verb = FALSE,
    CR_term = TRUE, curves, mats)

mc_fit_exo(yields, lambdas, cutoffs, reference, Xt, X_normalize = TRUE,
    CR_algo = zero_mean_covreg, CR_init = "adaptive", Xt_smooth = FALSE,
    smoother = c("ns", "bs", "rm", "hp", "henderson"), knot_count = 5,
    k_count = 9, k_pass = 3, hp_lambda = 1600, henderson_k = 13,
    CR_maxiter = 1000, CR_tol = 1e-08, CR_Binit = NULL, CR_S0init = NULL,
    CR_S0_shrink_diag = 0, CR_check_every = 1, CR_use_loglik = TRUE,
    CR_store_path = TRUE, CR_verb = FALSE, CR_term = TRUE, curves,
    mats)
}
\arguments{
\item{yields}{Named list of finite yield panels, each with a common unique,
increasing observation index first, followed by maturity columns. Labels are
matched to \code{mats} and curve identities to \code{curves}.}
\item{lambdas}{Positive candidate NS decays in inverse months. A scalar fixes
the decay; otherwise each curve selects its own grid minimum of yield SSE.}
\item{cutoffs}{Two integer maturity-column indices defining short, medium
and long residual covariance blocks; see \code{\link{bl_single_curve}}.}
\item{reference}{Reference curve for endogenous country-factor relationships.
Accepted but unused by \code{mc_fit_exo()}.}
\item{ECM_estim}{Estimation method passed to \code{tsDyn::VECM()}.}
\item{ECM_type}{Rank-test type passed to \code{tsDyn::rank.test()}.}
\item{ECM_alpha}{Critical-value/significance setting passed as \code{cval}
to the rank test. It must be supported by that dependency.}
\item{X_normalize}{Standardize endogenous projected spreads using sample
means and standard deviations. Accepted but inactive in \code{mc_fit_exo()}.}
\item{X_trunc}{Keep one spread per usable country ECM, choosing the column
with the largest adjustment-coefficient column norm. Otherwise keep all.}
\item{CR_algo}{Covariance estimator function, run separately by maturity.
It must accept \code{W}, \code{X}, \code{init}, \code{max_iter}, \code{tol},
\code{S0}, \code{B}, \code{S0_shrink_diag}, and return at least \code{B} and
\code{S0}. Additional named controls are forwarded only when declared in its
formals. The MCE default is \code{zero_mean_covreg_em}; the MCX default is
\code{zero_mean_covreg}. These algorithms are not interchangeable numerically.}
\item{CR_init}{Covariance loading initialization rule.}
\item{Xt_smooth}{Smooth feature columns before covariance estimation.}
\item{smoother}{Natural spline (\code{"ns"}), B-spline (\code{"bs"}),
repeated running mean (\code{"rm"}), Hodrick--Prescott (\code{"hp"}), or
Henderson smoothing (\code{"henderson"}).}
\item{knot_count}{Number of interior quantile knots for spline smoothing.}
\item{k_count, k_pass}{Running-mean window length and number of passes.}
\item{hp_lambda}{Nonnegative HP smoothing penalty.}
\item{henderson_k}{Positive odd Henderson filter length.}
\item{CR_maxiter, CR_tol}{Covariance iteration budget and convergence tolerance.}
\item{CR_Binit, CR_S0init}{Optional covariance loading and intercept covariance
initializations, passed to each maturity's estimator. Shapes are
country-by-feature and country-by-country, respectively.}
\item{CR_S0_shrink_diag}{Covariance diagonal shrinkage fraction in \eqn{[0,1]}.}
\item{CR_check_every}{Regularized covariance convergence-check interval.}
\item{CR_use_loglik}{Compute the regularized covariance likelihood score.}
\item{CR_store_path}{Store regularized covariance parameter paths. Together
with \code{CR_use_loglik}, enables the additional score-change stopping check.}
\item{CR_verb}{Request regularized covariance diagnostics through its
\code{verbose} argument. The legacy estimator does not print diagnostics.}
\item{CR_term}{Enable the legacy covariance estimator's termination checks.
Ignored by the regularized estimator, which always evaluates its own criterion.}
\item{curves}{At least two unique curve names, specifying output order.}
\item{mats}{Positive, unique, increasing maturities in months.}

\item{Xt}{For MCX, a finite numeric observation-by-feature matrix. Include
an intercept explicitly if wanted. Rows must align with the yield observation
index; a time column is not part of this numeric design.}
}
\details{
Both workflows select curve-specific NS decays, fit level/slope/curvature,
and collect observed-minus-fitted residuals across countries at each maturity.
For maturity \eqn{j} and feature vector \eqn{x_t}, covariance regression gives
\deqn{\Sigma_{j,t}=S_{0,j}+(B_j x_t)(B_j x_t)' .}
The tenor-specific covariances are assembled in tenor order for whitening,
transformed NS projection and rescaling back to yield units. Whitening requires
numerically invertible conditional covariances; loading regressions require
invertible cross-products. No pseudoinverse is substituted.

MCE constructs features from reference/country NS-factor pairs. Its lag order
is selected using the Schwarz criterion from \code{vars::VARselect()}, with
a minimum VAR order of two, before estimating and rank-testing the VECM.
Only ranks strictly between zero and the number of factor variables are used.
No usable relationships, zero-variance normalization, or failed ECM estimation
cause errors. Optional truncation selects by the adjustment-column rule above.

MCX uses the supplied \code{Xt} directly, apart from optional smoothing.
It does not estimate ECM relationships or normalize \code{Xt} automatically;
\code{X_normalize} and \code{reference} are currently inactive in this path.
Choosing \code{CR_algo} explicitly is necessary when comparing the same
covariance algorithm across MCE and MCX, because their defaults differ.

Spline smoothing is fitted over the supplied observation window. The running
mean, HP and Henderson paths also use the supplied window, including future
rows relative to an interior observation. These outputs are in-sample
reconstructions, not independently evaluated out-of-sample forecasts.

Finite, aligned inputs and percentage-yield/month units are described in
\code{\link{estimation_inputs}}. Nonconvergence warnings are annotated by
tenor; inspect the returned per-tenor status. Controls not supported by a
chosen covariance estimator have no effect on that estimator. See
\code{\link{zero_mean_covreg}} for stopping rules, ridge behavior and latent
moment timing. No fitting seeds are set by these wrappers.
}
\value{
A list with \code{tenor} and \code{curve} fitted yield panels (index first,
then maturity-major or curve-major series), \code{sigma_JT} containing
\code{tSigma} and maturity-specific \code{jSigma} covariance paths,
\code{BS0} containing the per-tenor estimator results,
the feature matrix actually used in \code{Xt}, per-tenor residual samples
\code{W}, and \code{input_spec} recording units, alignment and transformations.
MCE additionally returns \code{curveECM} and baseline factor estimates in
\code{ns_factor}; MCX does not return these ECM/factor fields.
}
\examples{
set.seed(21)
mats <- c(1, 3, 6, 12, 24, 60)
n <- 120
time <- as.Date("2020-01-01") + seq_len(n) - 1
beta <- cbind(2 + cumsum(rnorm(n, sd = .03)),
              -.5 + cumsum(rnorm(n, sd = .02)),
              .2 + cumsum(rnorm(n, sd = .01)))
yields <- setNames(lapply(c(0, .2), function(shift) {
  panel <- data.frame(time, (beta + shift) \%*\% NS_loadings(.35, mats) +
                      matrix(rnorm(n * length(mats), sd = .02), n))
  names(panel)[-1] <- c("01M", "03M", "06M", "01Y", "02Y", "05Y")
  panel
}), c("usa", "gbr"))
Xt <- cbind(intercept = 1, driver = seq(-1, 1, length.out = n))
fit <- mc_fit_exo(yields, .35, c(2, 4), "usa", Xt,
                  curves = names(yields), mats = mats,
                  CR_algo = zero_mean_covreg_em, CR_maxiter = 500,
                  CR_tol = 1e-5, CR_use_loglik = FALSE, CR_store_path = FALSE)
head(fit$curve)
vapply(fit$BS0, function(x) x$converged, logical(1))
\donttest{
# Endogenous fitting also requires usable estimated cointegration features.
endogenous <- mc_fit_end(yields, .35, c(2, 4), "usa",
                         curves = names(yields), mats = mats,
                         CR_maxiter = 500, CR_tol = 1e-5)
}
}
\seealso{
\code{\link{sc_fit}}, \code{\link{zero_mean_covreg}}, \code{\link{mc_fit_shock_base}}, \code{\link{estimation_inputs}}
}
