\name{estimation_inputs}
\alias{estimation_inputs}
\title{Alignment and units for estimation inputs}
\description{
Shared input conventions for single-curve, DLY, MCE, MCX and reusable
stress-baseline estimation. These checks do not impute observations, change
the estimation sample, normalize covariates, or convert rate conventions.
}
\details{
Multi-curve inputs are lists of data frames. Each panel has an observation index
first and finite numeric yields thereafter. Names identify curves; \code{curves}
specifies their order. Lists with no names can be labelled positionally by an
explicit \code{curves} argument. Partial, duplicate or conflicting names fail.
Named panels are ordered by identity, never relabelled as another country.

Observation indices must be nonmissing, unique and strictly increasing, and
must match across curves. Supported indices are \code{Date}, \code{POSIXct},
ISO date strings and numeric observation indices. Numeric indices are not
implicitly interpreted as Unix or spreadsheet dates. Date representations are
compared without replacing the supplied index in outputs. Market-calendar gaps
are allowed. Dynamic lags mean successive observations, not fixed calendar days.
Missing or mismatched dates must be resolved explicitly before fitting.
Natural- and B-spline feature smoothing use numeric observation indices
directly, retaining their supplied spacing without requiring a date origin.
Date-based spline smoothing retains its existing calendar-day coordinates.

Numeric \code{mats} are positive, unique, increasing maturities in months.
Each labelled panel must contain exactly those maturities. Equivalent labels,
including \code{12M}, \code{01Y} and \code{X01Y}, identify the same maturity;
duplicate economic tenors fail. Labelled columns may be reordered to match
\code{mats}, retaining their labels.
Unlabelled matrices accepted by lower-level helpers use positional columns.

The package workflow convention is percentage yields: \code{5} means five
percent. Zero and negative yields are allowed. Country factors and yield
residuals retain that scale; yield covariances use squared percentage points.
There is no magnitude-based unit inference or automatic rescaling. Standalone
covariance regression uses response and feature units as supplied.
Covariance regression requires at least two response columns in \code{W};
univariate responses are not supported. A single covariate column in \code{X}
remains valid.

Independent-curve yield residuals returned as \code{W} use observed minus
fitted yields, for both the three-factor NS and two-factor DLY country fits.
Shocked states retain the baseline residual sample; they do not recompute
residuals against the shocked yield panel.
A simultaneous sign reversal of all responses leaves the zero-mean covariance
model unchanged. With adaptive initialization, the legacy
\code{zero_mean_covreg()} backend can reverse \code{B} and signed \code{B x}
while retaining the same covariance. The regularized EM estimators normalize
the sign of \code{B}; adaptive initialization can nevertheless affect their
first convergence check, potentially changing the stopping iteration and final
covariance estimates when termination occurs at that check.
Signed latent moments and initialization histories
also depend on the residual convention. This convention change does not alter
the covariance algorithms, stopping rules or random-number seeds.
Previously saved fits retain their stored residual convention. Refit to obtain
the new convention, and do not combine residual samples using opposite signs.

The regularized EM estimators use \code{tol} for Frobenius-norm changes in
\code{B} and, when estimated, \code{S0}, with or without path storage. When
\code{use_loglik = TRUE} and \code{store_path = TRUE}, convergence additionally
requires finite successive likelihood scores whose absolute difference is below
\code{tol}. This branch previously used \code{sqrt(tol)} for parameter changes;
the common parameter threshold can change its stopping iteration and final
estimates. The additional likelihood check is retained, so matching parameter
thresholds alone does not guarantee identical stopping for both storage choices.

Covariance estimators warn with class \code{crosscurver_nonconvergence} when
the iteration limit is reached without confirmed convergence. The final fit is
still returned, without changing estimation or adding convergence checks.
Warnings identify the estimator, iteration count and limit, tolerance, sample
dimensions (\code{N} observations, \code{P} responses and \code{Q} features), active stopping criterion, checking
interval, last checked iteration, and available parameter-change and likelihood
diagnostics. Nonfinite scores and limits reached between checks are identified.
The legacy estimator reports changes from its existing history comparison when
\code{term = TRUE}; with \code{term = FALSE}, it explicitly warns that convergence
was not assessed. A limit before the first scheduled check is also distinguished
from a failed check. Fits confirmed converged on their last allowed iteration
do not warn.

In \code{zero_mean_covreg_em()}, \code{covreg_comp_estimator()} and the
legacy \code{zero_mean_covreg()}, \code{mz}, \code{vz} and \code{sz} are
the conditional latent mean, variance and standard deviation from the last
executed E step. The conditional second moment is \code{vz + mz^2}.
For \code{iter = k}, with \code{k >= 1}, these moments use the parameters
preceding the final M step; returned \code{B} and \code{S0} are the parameters
after that update. No additional E step refreshes the moments before return.
When parameter paths are returned, the moments correspond to
\code{B_path[[k]]} and \code{S0_path[[k]]}, whereas the returned parameters
occupy \code{[[k + 1]]}. In \code{covreg_comp_estimator()}, \code{Psi}
(also returned as \code{S0}) is fixed, so only \code{B} changes between these
points. Regularized E steps use the inverse of stored \code{S0} plus
\code{ridge_S0} times the identity; the legacy E step inverts stored \code{S0}
without that additional ridge.

This timing applies to convergence, iteration-limit termination and legacy
fixed-iteration runs. After one iteration, the moments use the initialized
parameters. \code{converged = TRUE} reports the implemented stopping criterion;
it does not mean the moments have been recomputed using returned parameters.
Differences can be material in short runs or with loose tolerances.
When calculated, \code{loglik} is evaluated using returned parameters and the
existing scoring convention.

The regularized estimators can normalize the sign of \code{B} after the
E step without changing the already calculated \code{mz}. Signed latent means
therefore retain their preceding E-step orientation and need not match the
orientation of returned \code{B}. The sign change itself does not affect
\code{vz} or \code{sz}. These are iteration diagnostics, rather than refreshed
conditional moments for the returned fit.

The legacy \code{zero_mean_covreg()} estimator returns the common fields
\code{B}, \code{S0}, \code{mz}, \code{vz}, \code{sz}, \code{iter},
\code{converged}, \code{loglik}, \code{S0_shrink_diag} and \code{input_spec}.
Its latent moments follow the shared timing convention above.
\code{loglik} is \code{NULL}: no likelihood score is computed.
\code{B_path} and \code{S0_path} contain initialization followed by completed
updates, with length \code{iter + 1}. The original \code{b_est} and \code{s_est}
histories remain available unchanged, including unused preallocated entries.
\code{termination_reason} is \code{"converged"}, \code{"max_iter"} or
\code{"not_assessed"}; \code{converged = FALSE} means convergence was not
confirmed, including runs with \code{term = FALSE}.
These output conventions preserve the legacy zero-mean model, initialization,
updates, ten-iteration history-check schedule and stopping criterion. Legacy
fitting still stores its histories and does not implement the regularized
estimator's \code{check_every}, \code{use_loglik} or \code{store_path} controls.

Tenor-specific fitting through \code{blcc_fe()} adds the tenor, maturity in months
and named responses to the same warning. Conditions retain these as
\code{tenor}, \code{maturity_months} and \code{response_names}; the underlying
diagnostics are available as \code{estimator}, \code{iter}, \code{max_iter},
\code{tol}, \code{dimensions}, \code{parameter_changes}, \code{check_every},
\code{last_check}, \code{criterion}, \code{loglik} and \code{loglik_change}.
The class can be handled with \code{withCallingHandlers()} without suppressing
unrelated warnings. Emitting a warning does not alter fitted estimates.

Within-curve \code{H} construction warns with class \code{crosscurver_h_not_psd}
when the assembled matrix has a minimum eigenvalue below minus
\code{1e-8} times its largest absolute eigenvalue. This relative tolerance avoids
warnings for negligible negative eigenvalues without imposing an absolute
covariance scale. PSD and singular PSD matrices do not warn. A nonfinite matrix
or unavailable eigenvalues instead produces a warning that PSD could not be
assessed. Warnings include the eigenvalue, tolerance, cutoffs, maturity grid and
observation count; multi-curve estimation also identifies the curve.
Conditions retain these as \code{minimum_eigenvalue}, \code{psd_tolerance},
\code{psd_assessed}, \code{cutoffs}, \code{mats}, \code{n_observations} and,
where available, \code{curve}. The original matrix is returned unchanged: this
is a diagnostic, not a PSD projection or a change to the estimation procedure.

Shock helpers align named factor columns to \code{L, S, C}, or \code{L, S}
for two-factor models, moving values with their identities. Internal basis names
\code{l1, l2, l3} are recognized aliases. Duplicate, blank, \code{NA}, unknown
or incompatible named factor layouts fail. Unnamed matrices retain their historical positional order;
the average-level profile helper continues to require named factors.
Supplied shock magnitudes, profile magnitudes and profile fractions must be
finite: \code{NA}, \code{NaN}, \code{Inf} and \code{-Inf} raise errors.
Finite zero and negative shocks remain supported.

DLY global and country residual dynamic regressions stop if estimated
coefficients are nonfinite, including unidentified coefficients from
rank-deficient designs. Such coefficients are not replaced with zero.
Finite-coefficient regressions retain their existing specification and estimates.

NS, Svensson, DLY country and transformed NS loading regressions require an
invertible normal-equation matrix. Nonfinite matrix entries, failed inversions
and nonfinite inverse entries raise errors identifying the regression and maturity grid, together
with decay parameters or the observation number. Failed inversions also report
the reciprocal condition estimate when available and the original solver error.
The existing \code{Matrix::solve()} call and its singularity criteria are
retained; no pseudoinverse, additional regularization or new condition-number
cutoff is introduced. A successful finite inversion does not guarantee that
the regression is well conditioned.

NS decay parameters multiply maturity and are in inverse months. Svensson
parameters named \code{lambda1}, \code{lambda2} or supplied in \code{slambdas}
are time constants in months. Both must be finite and positive. Changing
maturities from months to years requires multiplying NS decay by twelve and
dividing Svensson time constants by twelve. Fixed decay inputs are validated
just like grids. \code{cutoffs} are maturity-column indices, not durations.

No estimation function converts par yields, yields to maturity or compounding
conventions. The application pricing examples interpret rates as continuously
compounded zero yields, convert percentage rates to decimals and months to
years, then discount. Users must establish that their source data support
that interpretation.

\code{Xt} contains one finite numeric row per yield observation. Nondefault
date row names allow the row index to be checked; numeric data-frame row labels
are positional even after subsetting. Unindexed matrices are
explicitly positional. Feature names, when supplied, must be unique. Covariance
prediction orders named features against the fitted design and rejects
conflicting names. Initial \code{B} and \code{S0} matrices are dimension-checked
and matched to named response and feature identities where available.
Intercepts and normalization are never inserted automatically. Regularized
covariance estimators require full-column-rank \code{X} before
ridge stabilization; increasing \code{ridge_X} does not bypass that check.
The existing
\code{mc_fit_exo()} arguments \code{reference} and \code{X_normalize} remain
inactive compatibility arguments.

Complete fits record an \code{input_spec} containing the observation index,
curve and maturity order, declared units, decay choices and feature layout. Endogenous fitting
also records the feature bases, centering, scaling and smoothing specification.
These describe transformations actually used. Scenario standardization uses
scenario-sample means and standard deviations. ECM and covariance-regression
conditioning remain independently selectable; feature-column alignment checks
do not establish statistical equivalence between different feature bases.
Spread normalization errors identify columns with zero or nonfinite standard
deviations; unnormalized constant spreads remain permitted at that stage.
ECM errors retain the original dependency message and add observation and
variable counts, plus the selected VECM lag when available. These diagnostics
do not impose a new sample-size threshold or change lag selection or the
configured ECM fallback policy.
See \code{\link{feature_reuse}}.
Precomputed DLY inputs must satisfy the same alignment and matrix-layout checks.
}
\seealso{\code{\link{dly_fit}}, \code{\link{dr_read_yields}}}
