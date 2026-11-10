\name{st_dly_shock}
\alias{st_dly_shock}
\alias{st_dly_shock_from_fit_profile}
\alias{st_dly_global_shock}
\alias{st_dly_global_shock_from_fit_profile}
\title{Reconstruct DLY country or global factor shocks}
\description{
Apply factor profiles with the fitted PCA basis and country loadings
held fixed. Country and global shocks use distinct reconstruction rules.
}
\usage{
st_dly_shock(dly_fit_obj, shock_curves, shock_factors, shock_profile,
    shock_type = c("additive", "multiplicative"))

st_dly_global_shock(dly_fit_obj, shock_factors, shock_profile, shock_type = c("additive",
    "multiplicative"))
}
\arguments{
\item{dly_fit_obj}{A complete baseline fit returned by \code{dly_fit()} or
\code{dly_fit_from_sc()}, including the saved PCA objects and country models.}
\item{shock_curves}{Subset of baseline countries directly shocked by
\code{st_dly_shock()}.}
\item{shock_factors}{One of L or S for country shocks. Global shocks accept
L and/or S; the same magnitude profile is applied to each requested factor.}
\item{shock_profile}{Finite numeric magnitudes with one value per baseline
observation, or a data frame containing \code{shock_magnitude}. Data-frame
dates and metadata are returned but are not used to match or reorder rows.}
\item{shock_type}{Country shocks support additive or multiplicative factor
changes. Global shocks support additive changes only; requesting multiplicative
global shocks is an error.}
}
\details{
Write baseline country factors as
\eqn{f_{i,t}=a_i+b_i g_t+u_{i,t}}, separately for level and slope.

Country shocks first modify the requested raw country factors and project
the modified country panel through the saved PCA centering, scaling and
rotation to obtain \eqn{g'_t}. They keep the baseline country intercepts and
loadings. For directly shocked countries the idiosyncratic component becomes
\eqn{u'_{i,t}=u_{i,t}+(f'_{i,t}-f_{i,t})}; it remains unchanged elsewhere.
The reconstructed combined factor is \eqn{a_i+b_i g'_t+u'_{i,t}}.
Thus the combined reconstruction is a specified stress rule, rather than a
new residual decomposition that necessarily reproduces every shocked raw factor.

Global shocks convert the additive country-average magnitude \eqn{\delta_t}
to a global-factor shift \eqn{\delta_t/\bar b}, with \eqn{\bar b} the mean
saved country loading for that factor. A nonfinite or near-zero mean loading
is rejected. Baseline raw country curves and idiosyncratic components are
retained, while global and combined yield components are reconstructed.

Neither workflow re-estimates PCA, decay, country loadings, or dynamic
parameters. The country-shock fit is a static reconstruction without dynamic
fields. The global-shock helper copies the baseline fit and updates its static
components; retained dynamics and propagated-path fields remain baseline
quantities. These helpers do not regenerate a dynamic stress path or an
out-of-sample forecast.

\code{st_dly_shock_from_fit_profile()} and
\code{st_dly_global_shock_from_fit_profile()} are compatibility aliases for
the country and global functions, preserving arguments and return layout.
}
\value{
A list with reconstructed \code{fit}, the original \code{shock_profile},
and \code{baseline_fit}. Country shocks also return \code{shocked_factors}.
Within \code{fit$yields}, use the requested \code{raw}, \code{global},
\code{idiosyncratic}, or \code{global_plus_idiosyncratic} static component.
Country and global shock results do not have identical nested fit structures;
in particular, only the global path retains baseline dynamic fields.
}
\examples{
set.seed(31)
mats <- c(1, 3, 6, 12, 24, 60)
n <- 80
time <- as.Date("2020-01-01") + seq_len(n) - 1
global <- cbind(2 + cumsum(rnorm(n, sd = .03)),
                -.5 + cumsum(rnorm(n, sd = .02)))
yields <- setNames(lapply(c(0, .2, .4), function(shift) {
  beta <- global + shift + matrix(rnorm(2 * n, sd = .015), n)
  panel <- data.frame(time, beta \%*\% NS_loadings_ls(.35, mats))
  names(panel)[-1] <- c("01M", "03M", "06M", "01Y", "02Y", "05Y")
  panel
}), c("usa", "gbr", "jpn"))
fit <- dly_fit(yields, mats, lambdas = .35, reference = "usa")

profile <- data.frame(time, shock_magnitude = seq(0, .1, length.out = n))
country <- st_dly_shock(fit, "gbr", "L", profile)
global_shock <- st_dly_global_shock(fit, "L", profile)
head(country$fit$yields$global_plus_idiosyncratic$gbr)
}
\seealso{
\code{\link{dly_fit}}, \code{\link{st_factor_profile}}, \code{\link{shock_fits}}
}
