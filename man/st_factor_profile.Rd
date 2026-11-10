\name{st_factor_profile}
\alias{st_factor_profile}
\alias{st_apply_factor_shock_profile}
\alias{st_mean_profile}
\alias{st_average_level_shock_profile}
\title{Construct observation-specific factor shock profiles}
\description{
Apply an aligned shock profile to country factor matrices, or convert
fractions of cross-country mean factor values to additive magnitudes.
}
\usage{
st_factor_profile(curves, curve_dns_factor, shock_curves, shock_factors,
    shock_magnitudes, shock_type = c("additive", "multiplicative"))

st_mean_profile(curve_factor, shock_fraction_profile, shock_factor = "L",
    curves = names(curve_factor))
}
\arguments{
\item{curves}{Model curve names. For \code{st_mean_profile()}, the country
subset used to compute the cross-sectional mean.}
\item{curve_dns_factor, curve_factor}{Named lists of observation-by-factor
matrices in common row order. Profile application supports L/S or L/S/C;
mean conversion also accepts a single selected factor.}
\item{shock_curves}{Subset of model curves whose factors are changed.}
\item{shock_factors}{Exactly one factor identity compatible with the matrices.}
\item{shock_magnitudes}{Finite numeric vector, one value per observation.
Magnitudes are additive factor units or proportional changes, depending on
\code{shock_type}.}
\item{shock_type}{Add magnitudes, or multiply affected factors by
\code{1 + shock_magnitudes}. For example, \code{.1} means a ten-percent
multiplicative change, or 0.1 yield-percentage-point additive change under
the package's percentage-yield convention.}
\item{shock_fraction_profile}{Finite numeric vector of proportional fractions
to multiply by the selected countries' mean factor value on each observation.}
\item{shock_factor}{The factor used for mean conversion: L, S or C.}
}
\details{
Rows are aligned positionally; these helpers do not match dates, interpolate
profiles, or estimate model parameters. Named factor columns are aligned to
L/S/C or L/S, including recognition of internal \code{l1/l2/l3} names.
Ambiguous labels fail. Unnamed matrices use the canonical positional order.

\code{st_factor_profile()} changes one factor on selected curves and retains
the remaining factor values. It does not reconstruct yields or fit covariance
regression. Pass its \code{shocked_factors} to \code{mc_fit_shocked()} for
the shared MCE scenario workflow.

\code{st_mean_profile()} computes
\deqn{\delta_t=q_t\frac{1}{|C|}\sum_{i\in C} f_{i,t},}
where \eqn{q_t} is the supplied fraction and \eqn{f_{i,t}} the selected
factor. This is a cross-country mean at each observation, not a time-window
mean. Negative mean factor values retain their sign.

\code{st_apply_factor_shock_profile()} and
\code{st_average_level_shock_profile()} are compatibility aliases with the
same arguments and results as the profile and mean functions, respectively.
}
\value{
Profile application returns a list with numeric \code{shock_profile} and
a named list \code{shocked_factors}. Mean conversion returns a numeric
vector of additive magnitudes, one per observation.
}
\examples{
factors <- list(usa = cbind(L = c(2, 2.1, 2.2), S = c(-.5, -.4, -.3)),
                gbr = cbind(L = c(3, 3.1, 3.2), S = c(-.6, -.5, -.4)))
magnitude <- st_mean_profile(factors, c(0, .05, .1), shock_factor = "L")
scenario <- st_factor_profile(names(factors), factors, "gbr", "L", magnitude)
scenario$shocked_factors$gbr
}
\seealso{
\code{\link{mc_fit_shocked}}, \code{\link{sc_fit_shock}}, \code{\link{st_dly_shock}}
}
