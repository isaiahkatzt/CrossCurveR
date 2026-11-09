\name{feature_reuse}
\alias{feature_reuse}
\title{Cached MCE controls and stress conditioning}
\description{
Control inheritance, structural coefficient checks, and independent conditioning
choices for cached MCE stress scenarios.
}
\details{
\code{mc_fit_shock(baseline_state = ...)} inherits omitted fitting controls from
the saved baseline, even when they differ from the wrapper's defaults. Explicit
controls are forwarded to \code{mc_fit_shock_from_baseline()}; its existing
\code{NULL} convention means inheritance. Calls without a cache retain the
ordinary baseline-fitting defaults. Supplied baseline-only inputs (yields,
decay grid, cutoffs, curve order, maturity grid and reference) must agree with
the cache. A new baseline is needed to change these calibration inputs.

\code{reuse_ECM} and \code{reuse_BS0} independently choose whether to condition
on baseline cointegration/reversion estimates and covariance-regression
coefficients. In particular, \code{reuse_ECM = FALSE, reuse_BS0 = TRUE}
re-estimates the ECM while retaining covariance coefficients. Reused
coefficients require the same feature count and alignable feature names;
a changed rank or selected-column layout can make this combination undefined.
These are structural checks, not a claim that different feature bases have
the same statistical interpretation.

Normalization, selection and smoothing remain explicit scenario controls.
With \code{X_normalize = TRUE}, projected cointegration spreads are
standardized using scenario-sample statistics:
\deqn{x_t^{scenario} = (u_t^{scenario} - \bar{u}^{scenario}) / s^{scenario}.}
This preserves the application convention preceding the feature-reuse pass.
The recent baseline-relative standardization and restrictions on feature
controls have been reverted. Sample standardization is separate from the
identification normalization of the cointegration matrix, and the manuscript's
projected-spread equations do not specify this sample-standardization policy.

Scenario results expose the covariance parameters actually used in \code{BS0}
and the corresponding \code{xt_state}, including its \code{feature_spec}.
Caches do not need transformation metadata merely to reuse coefficients.
Available baseline metadata is checked for internal consistency. An explicitly
requested baseline ECM fallback also covers errors during ECM estimation.

This reversion preserves the pre-existing observed-versus-DNS-fitted yield
convention, baseline residual reuse, covariance-intercept shrinkage, and
spread-selection rule. It does not resolve those separate paper/code issues.
Zero shocks preserve features and conditional covariances when coefficients
are reused; equality of baseline and shocked fitted yields is a separate issue.
}
\seealso{\code{\link{estimation_inputs}}}
