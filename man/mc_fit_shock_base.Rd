\name{mc_fit_shock_base}
\alias{mc_fit_shock_base}
\alias{mc_fit_shock_baseline}
\alias{mc_fit_shock_cached}
\alias{mc_fit_shock_from_baseline}
\alias{mc_fit_shock_refit}
\alias{mc_fit_shock_refit_from_baseline}
\alias{mc_fit_shock}
\alias{mc_fit_shocked}
\title{Prepare and reuse MCE stress-test baselines}
\description{
Fit a reusable endogenous multi-curve baseline, reconstruct factor shocks,
and independently reuse or refit ECM and covariance coefficients.
}
\usage{
mc_fit_shock_base(yields, lambdas, cutoffs, reference, ECM_estim = "ML",
    ECM_type = "eigen", ECM_alpha = 0.1, X_normalize = TRUE,
    X_trunc = FALSE, CR_algo = zero_mean_covreg_em, CR_init = "adaptive",
    Xt_smooth = FALSE, smoother = c("ns", "bs", "rm", "hp", "henderson"),
    knot_count = 5, k_count = 9, k_pass = 3, hp_lambda = 1600,
    henderson_k = 13, CR_maxiter = 1000, CR_tol = 1e-08, CR_Binit = NULL,
    CR_S0init = NULL, CR_S0_shrink_diag = 0, CR_check_every = 1,
    CR_use_loglik = TRUE, CR_store_path = TRUE, CR_verb = FALSE,
    CR_term = TRUE, curves, mats)

mc_fit_shock_cached(baseline_state, shock_curves, shock_factors, shock_magnitude,
    shock_type = c("additive", "multiplicative"), shock_window = NULL,
    shock_window_position = c("first", "last"), reuse_BS0 = TRUE,
    reuse_ECM = TRUE, ECM_estim = NULL, ECM_type = NULL, ECM_alpha = NULL,
    ECM_fallback = c("error", "baseline"), X_normalize = NULL,
    X_trunc = NULL, CR_algo = NULL, CR_init = NULL, Xt_smooth = NULL,
    smoother = NULL, knot_count = NULL, k_count = NULL, k_pass = NULL,
    hp_lambda = NULL, henderson_k = NULL, CR_maxiter = NULL,
    CR_tol = NULL, CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = NULL,
    CR_check_every = NULL, CR_use_loglik = NULL, CR_store_path = NULL,
    CR_verb = NULL, CR_term = NULL)

mc_fit_shock_refit(baseline_state, shock_curves, shock_factors, shock_magnitude,
    shock_type = c("additive", "multiplicative"), shock_window = NULL,
    shock_window_position = c("first", "last"), ECM_estim = NULL,
    ECM_type = NULL, ECM_alpha = NULL, ECM_fallback = c("error",
        "baseline"), X_normalize = NULL, X_trunc = NULL, CR_algo = NULL,
    CR_init = NULL, Xt_smooth = NULL, smoother = NULL, knot_count = NULL,
    k_count = NULL, k_pass = NULL, hp_lambda = NULL, henderson_k = NULL,
    CR_maxiter = NULL, CR_tol = NULL, CR_Binit = NULL, CR_S0init = NULL,
    CR_S0_shrink_diag = NULL, CR_check_every = NULL, CR_use_loglik = NULL,
    CR_store_path = NULL, CR_verb = NULL, CR_term = NULL)

mc_fit_shock(yields = NULL, lambdas = NULL, cutoffs = NULL, reference = NULL,
    shock_curves, shock_factors, shock_magnitude, shock_type = c("additive",
        "multiplicative"), shock_window = NULL, shock_window_position = c("first",
        "last"), reuse_BS0 = TRUE, reuse_ECM = TRUE, ECM_estim = "ML",
    ECM_type = "eigen", ECM_alpha = 0.1, ECM_fallback = c("error",
        "baseline"), X_normalize = TRUE, X_trunc = FALSE, CR_algo = zero_mean_covreg_em,
    CR_init = "adaptive", Xt_smooth = FALSE, smoother = c("ns",
        "bs", "rm", "hp", "henderson"), knot_count = 5, k_count = 9,
    k_pass = 3, hp_lambda = 1600, henderson_k = 13, CR_maxiter = 1000,
    CR_tol = 1e-08, CR_Binit = NULL, CR_S0init = NULL, CR_S0_shrink_diag = 0,
    CR_check_every = 1, CR_use_loglik = TRUE, CR_store_path = TRUE,
    CR_verb = FALSE, CR_term = TRUE, curves = NULL, mats = NULL,
    baseline_state = NULL)

mc_fit_shocked(baseline_state, shocked_factors, reuse_BS0 = TRUE,
    reuse_ECM = TRUE, ECM_fallback = c("error", "baseline"),
    controls = list())
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

\item{baseline_state}{Reusable baseline from \code{mc_fit_shock_base()}.
Treat its calibration contents as a cache and pass it back unchanged.
\code{mc_fit_shock()} may instead construct one from raw inputs.}
\item{shock_curves}{Names of curves to shock.}
\item{shock_factors}{Identities \code{"L"}, \code{"S"}, and/or \code{"C"}.}
\item{shock_magnitude}{Finite scalar for one factor, or a finite vector
matching \code{shock_factors}, optionally named by factor. Additive magnitudes
use the factor/yield units; multiplicative magnitudes are proportional changes.}
\item{shock_type}{Add factor magnitudes or multiply affected values by
\code{1 + shock_magnitude}.}
\item{shock_window}{Positive number of successive observations to shock.
\code{NULL} shocks the whole estimation window. This counts rows, not days.}
\item{shock_window_position}{Apply a partial shock window to the first or
last rows of the saved sample.}
\item{reuse_BS0, reuse_ECM}{Independently retain covariance coefficients
and baseline ECM relationships. \code{mc_fit_shock_refit()} sets both to
\code{FALSE} while retaining the baseline single-curve calibration.}
\item{ECM_fallback}{When ECM refitting or feature construction fails,
\code{"error"} propagates the failure; \code{"baseline"} warns and rebuilds
features using the saved ECM relationships and the scenario's factors.}
\item{shocked_factors}{For \code{mc_fit_shocked()}, a named list of
observation-by-factor country matrices matching baseline curves and rows.
These are already-shocked factor values, not magnitudes.}
\item{controls}{For \code{mc_fit_shocked()}, a named list of ECM, X, CR and
smoothing controls documented here. Omitted or \code{NULL} controls inherit
the cache. Raw calibration inputs are not control overrides.}
}
\details{
\code{mc_fit_shock_base()} prepares the independent DNS, ECM and covariance
calibration and saves the fitted baseline. This stress API is endogenous MCE;
it does not accept an exogenous \code{Xt} design.

\code{mc_fit_shock_cached()} applies a scalar/vector factor shock over the
selected rows and reconstructs through the saved NS loadings. Baseline
observed-minus-fitted residual samples are retained. Covariance reuse and ECM
reuse are separate choices; neither implies the other. When covariance
coefficients are retained, feature count and identities must remain compatible.

\code{mc_fit_shock_refit()} requests ECM and covariance refitting from that
same cached single-curve calibration. It does not re-select NS decays or refit
the original independent curve parameters.

\code{mc_fit_shock()} prepares a baseline if absent, then dispatches to the
cached workflow. With a supplied cache, any supplied raw calibration inputs
must match it. Omitted fitting controls inherit the cache even where the
wrapper's displayed defaults differ; explicit \code{NULL} also means inherit.
Cached/refit functions use \code{NULL} defaults for this purpose. See
\code{\link{feature_reuse}} for feature compatibility and conditioning.

\code{mc_fit_shocked()} is the shared workflow for already-shocked factor
matrices, including observation-specific profile adapters. Factor matrices
follow baseline observation order; no dates are matched or observations
removed. Labelled factors are aligned by identity to L/S/C, recognizing
internal \code{l1/l2/l3} labels. Unnamed matrices use positional columns.

Endogenous normalization uses the scenario sample's mean and standard
deviation; it does not freeze baseline centering/scaling. Conditioning on
saved ECM and covariance coefficients therefore preserves their identities
without requiring statistically identical feature transformations.
Zero shocks do not in general guarantee equality of baseline and shocked
yield panels: the baseline fits observed yields, whereas shock reconstruction
starts from the independent DNS fitted path and retains baseline residuals.

The original \code{mc_fit_shock_baseline()},
\code{mc_fit_shock_from_baseline()} and
\code{mc_fit_shock_refit_from_baseline()} names are compatibility aliases for
the base, cached and refit functions, with identical arguments and results.
Intermediate state-building helpers are internal.
}
\value{
Baseline construction returns a list of class \code{mc_shock_baseline},
containing curve/maturity/reference metadata, independent curve calibration,
feature state, per-tenor \code{BS0}, \code{baseline_fit}, \code{config},
and \code{input_spec}. This is a calibration cache, not a forecast object.

Shock fitting returns baseline \code{tenor}/\code{curve} panels;
\code{shocked_tenor}/\code{shocked_curve} panels; shocked independent
\code{shocked_yields} and \code{ns_factor}; \code{shocked_factors};
\code{sigma_JT}, \code{BS0}, \code{Xt}, \code{W}; and the saved baseline,
scenario \code{xt_state}, \code{shocked_single_curve_state}, and
\code{input_spec}. Scalar adapters additionally attach \code{shock_matrix}.
Inspect per-tenor \code{BS0} for convergence status and \code{xt_state} for
the transformations actually used. Nested state fields support the current
workflow; applications should prefer these public fitting functions to
direct state construction.
}
\examples{
\dontrun{
# With aligned observed yield panels and usable country-factor relationships:
baseline <- mc_fit_shock_base(yields, lambdas = c(.03, .05, .07),
                              cutoffs = c(2, 4), reference = "usa",
                              curves = names(yields), mats = mats)
scenario <- mc_fit_shock_cached(baseline, shock_curves = "gbr",
                                shock_factors = "L", shock_magnitude = .1)
refitted <- mc_fit_shock_refit(baseline, "gbr", "L", .1)

# An observation-specific profile uses the same MCE fitting workflow:
profile <- seq(0, .1, length.out = nrow(baseline$baseline_fit$curve))
factors <- st_factor_profile(baseline$curves,
                             baseline$single_curve_state$ns_factor,
                             "gbr", "L", profile)
profile_fit <- mc_fit_shocked(baseline, factors$shocked_factors)
}
}
\seealso{
\code{\link{mc_fit_end}}, \code{\link{st_factor_profile}}, \code{\link{feature_reuse}}, \code{\link{shock_fits}}
}
