\name{sc_fit_shock}
\alias{sc_fit_shock}
\title{Fit and shock independent DNS curves}
\description{
Prepare an independent DNS baseline and reconstruct curves under a factor
shock profile, or reuse a baseline returned by an earlier call.
}
\usage{
sc_fit_shock(yields = NULL, lambdas = NULL, cutoffs = NULL,
  curves = NULL, mats = NULL, baseline_state = NULL,
  shock_curves = NULL, shock_factors = "L", shock_profile = NULL,
  shock_type = c("additive", "multiplicative"))
}
\arguments{
\item{yields}{Named list of yield panels, each containing an observation index
and maturity columns. Required when \code{baseline_state} is absent.}
\item{lambdas}{Positive candidate NS decay parameters in inverse months.}
\item{cutoffs}{Two maturity-column indices defining the baseline residual
covariance blocks.}
\item{curves}{Curve names in the required output order.}
\item{mats}{Maturities in months.}
\item{baseline_state}{The \code{baseline_state} returned by this function.
When supplied, the saved calibration is reused and raw fitting inputs are ignored.}
\item{shock_curves}{Names of curves to shock. With \code{NULL}, no curves
are targeted even if a profile is supplied.}
\item{shock_factors}{One factor identity: \code{"L"}, \code{"S"}, or \code{"C"}.}
\item{shock_profile}{A finite numeric vector with one magnitude per observation,
or a data frame containing that vector in \code{shock_magnitude}. Rows correspond
positionally to the baseline observations; dates are not matched or reordered.
\code{NULL} returns unshocked DNS curves and a reusable baseline.}
\item{shock_type}{Additive magnitudes or proportional changes to factor values.}
}
\details{
This workflow uses the existing independent NS calibration, factor-shock and
curve-reconstruction calculations. It does not estimate cross-curve covariance
regression or DLY dynamics. Baseline residuals and residual covariance estimates
are retained when applying shocks; decay parameters are not re-estimated.
See \code{\link{bl_single_curve}} for the cutoff and residual covariance
construction and \code{\link{st_factor_profile}} for profile conventions.
}
\value{
A list containing reconstructed \code{yields}, \code{ns_factor}, the supplied
\code{shock_profile}, \code{shocked_factors}, and a reusable \code{baseline_state}.
The baseline is an opaque calibration cache to pass back to this function;
its internal components are not a separately supported interface.
}
\examples{
mats <- c(1, 3, 6, 12, 24, 60)
t <- seq_len(20)
beta <- cbind(2 + t / 100, -.5 + sin(t) / 100, .2 + cos(t) / 100)
panel <- data.frame(time = as.Date("2020-01-01") + t - 1,
                    beta \%*\% NS_loadings(.35, mats))
names(panel)[-1] <- c("01M", "03M", "06M", "01Y", "02Y", "05Y")
fit <- sc_fit_shock(list(usa = panel), .35, c(2, 4), "usa", mats)
shock <- sc_fit_shock(baseline_state = fit$baseline_state,
                     shock_curves = "usa", shock_profile = rep(.1, length(t)))
}
\seealso{\code{\link{shock_fits}}}
