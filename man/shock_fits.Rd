\name{shock_fits}
\alias{shock_fits}
\title{Yield-model stress-test workflows}
\description{
Overview of the current DNS, PCA-based DLY and MCE shock interfaces.
}
\details{
\code{\link{sc_fit_shock}} prepares independent DNS curves and reconstructs
factor profiles with optional reuse of its saved baseline.

\code{\link{mc_fit_shock_base}} documents MCE baseline construction,
cached shocks, refitting, the convenience \code{mc_fit_shock()} wrapper,
and the shared \code{mc_fit_shocked()} workflow for already-shocked factors.
ECM and covariance reuse are independent controls.

\code{\link{st_factor_profile}} applies observation-specific factor shocks
and converts fractions of cross-country mean factors to additive magnitudes.

\code{\link{st_dly_shock}} documents the distinct country and global DLY
reconstruction rules, including which baseline quantities remain fixed.

All profiles follow the supplied observation order. Additive shocks have the
same units as their factors; multiplicative shocks are proportional changes.
These are in-sample reconstruction and stress scenarios, not independently
evaluated forecasts. Public compatibility aliases are documented with their
corresponding functions. State construction and matrix-layout helpers are
internal, and the current export set is not a finalized publication API.
}
\seealso{
\code{\link{feature_reuse}}, \code{\link{estimation_inputs}}
}
