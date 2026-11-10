\name{bl_single_curve}
\alias{bl_single_curve}
\title{Fit an NS baseline and its structured residual covariance}
\description{
Fit three NS factors at a fixed decay and construct the within-curve
residual covariance used in multi-curve baseline preparation.
}
\usage{
bl_single_curve(lambda_select, yields, cutoffs, mats)
}
\arguments{
\item{lambda_select}{A positive finite NS decay in inverse months.}
\item{yields}{A data frame with an observation index first and finite numeric
yield columns matched to \code{mats}.}
\item{cutoffs}{Two integer maturity-column indices \eqn{c_1<c_2}, defining
short, medium and long blocks. Require \eqn{1\leq c_1<c_2<p}, where
\eqn{p} is the number of maturities. These are indices, not maturity values.}
\item{mats}{Positive, unique, increasing maturities in months.}
}
\details{
Coefficients are ordinary least squares estimates using
\code{NS_loadings(lambda_select, mats)}. Residuals are observed minus fitted.
Loading cross-products must be invertible.

\code{H} contains sample covariance blocks for the three maturity groups and
the adjacent cross-block covariances at the two cutoffs. Other cross-block
entries are zero. Although each sample covariance block is positive
semidefinite, the resulting structured \code{H} need not be. A materially
indefinite or unassessable \code{H} produces a
\code{crosscurver_h_not_psd} warning and is returned unchanged.
At least two observations are needed for finite sample covariances.
}
\value{
A list with \code{phi} (loadings \code{phi} and the OLS map \code{cross_phi}),
\code{betas} (\code{NSbetas} and observed \code{m_yield}),
\code{nsfit} (fitted \code{NSfit} and residuals \code{W}), and
the maturity-by-maturity covariance matrix \code{H}.
}
\examples{
mats <- c(1, 3, 6, 12, 24, 60)
t <- seq_len(20)
beta <- cbind(2 + t / 100, -.5 + sin(t) / 100, .2 + cos(t) / 100)
panel <- data.frame(time = t, beta \%*\% NS_loadings(.35, mats))
names(panel)[-1] <- c("01M", "03M", "06M", "01Y", "02Y", "05Y")
baseline <- bl_single_curve(.35, panel, c(2, 4), mats)
max(abs(baseline$nsfit$W))
}
\seealso{
\code{\link{sc_fit}}, \code{\link{mc_fit_end}}, \code{\link{estimation_inputs}}
}
