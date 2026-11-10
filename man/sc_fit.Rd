\name{sc_fit}
\alias{sc_fit}
\title{Fit an independent NS or Svensson yield curve}
\description{
Estimate factors for each observation using a window-selected NS decay,
daily nonlinear NS estimation, or supplied Svensson time constants.
}
\usage{
sc_fit(yields, mats, type = c("nelson", "svensson", "daily"),
    ts = TRUE, lambdas = seq(from = 0.001, to = 1, length.out = 100),
    slambdas = NULL)
}
\arguments{
\item{yields}{A data frame with an observation index first and finite numeric
yield columns. Labelled maturities are matched to \code{mats}.}
\item{mats}{Positive, unique, increasing maturities in months.}
\item{type}{\code{"nelson"} selects an NS decay over the supplied window;
\code{"daily"} estimates NS parameters separately each day;
\code{"svensson"} uses the supplied time constants.}
\item{ts}{Whether to fit a separate \code{stats::ar()} model to each estimated
factor series. These are in-sample fits with the dependency's default controls.}
\item{lambdas}{Positive candidate NS decays in inverse months. A scalar fixes
the decay. Used only for \code{type = "nelson"}.}
\item{slambdas}{Two positive Svensson time constants in months, required for
\code{type = "svensson"}.}
}
\details{
The \code{"nelson"} path minimizes the sum of squared yield residuals over
all supplied observations for each candidate decay, then fits three factors
by ordinary least squares for that decay. The \code{"svensson"} path fits
four factors by ordinary least squares with fixed time constants. Both require
invertible loading cross-products; no pseudoinverse is substituted.

The \code{"daily"} path uses \code{YieldCurve::Nelson.Siegel()} on the
observation-by-maturity yield matrix. Its nonlinear calibration differs from
the fixed-grid NS path. Yields and factors retain their input scale; the
package workflow convention is percentage yields.
See \code{\link{estimation_inputs}} for alignment and units.
}
\value{
A list with fitted \code{yields} and \code{betas} matrices,
\code{dynamic} containing the per-factor autoregressive fits (or \code{NA}
when \code{ts = FALSE}), and \code{input_spec} recording input and selected
decay conventions. Observation indices are stored in \code{input_spec}; the
two fitted matrices have one row per input observation.
}
\examples{
mats <- c(1, 3, 6, 12, 24, 60)
t <- seq_len(30)
beta <- cbind(2 + t / 100, -.5 + sin(t) / 100, .2 + cos(t) / 100)
panel <- data.frame(time = as.Date("2020-01-01") + t - 1,
                    beta \%*\% NS_loadings(.35, mats))
names(panel)[-1] <- c("01M", "03M", "06M", "01Y", "02Y", "05Y")
fit <- sc_fit(panel, mats, lambdas = .35, ts = FALSE)
head(fit$betas)
}
\seealso{
\code{\link{NS_loadings}}, \code{\link{bl_single_curve}}, \code{\link{sc_fit_shock}}
}
