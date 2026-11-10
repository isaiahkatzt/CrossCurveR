\name{NS_loadings}
\alias{NS_loadings}
\alias{NS_loadings_ls}
\alias{Svensson_loadings}
\title{Nelson--Siegel and Svensson factor loadings}
\description{
Construct loading matrices for three-factor NS, two-factor level/slope NS,
and four-factor Svensson yield curves.
}
\usage{
NS_loadings(lambda, mats)

NS_loadings_ls(lambda, mats)

Svensson_loadings(lambda1, lambda2, mats)
}
\arguments{
\item{lambda}{A positive finite NS decay parameter in inverse months.}
\item{lambda1, lambda2}{Positive finite Svensson time constants in months.}
\item{mats}{Positive, unique, increasing maturities in months.}
}
\details{
For maturity \eqn{m}, define \eqn{a(z)=(1-\exp(-z))/z}.
The three NS loadings are \eqn{1}, \eqn{a(\lambda m)}, and
\eqn{a(\lambda m)-\exp(-\lambda m)}.
\code{NS_loadings_ls()} returns only the first two.

Svensson uses \eqn{1}, \eqn{a(m/\lambda_1)},
\eqn{a(m/\lambda_1)-\exp(-m/\lambda_1)}, and
\eqn{a(m/\lambda_2)-\exp(-m/\lambda_2)}. Its decay arguments are time
constants, whereas NS takes an inverse time constant. They are not
interchangeable without converting this parameterization.

For an observation-by-factor coefficient matrix \code{beta}, fitted yields
are \code{beta \%*\% loadings}, with factors in the stated row order.
These functions construct loadings only; they do not fit coefficients.
}
\value{
A numeric matrix with one column per maturity and two, three or four rows.
Rows represent level, slope, then the applicable curvature factors, and are
named \code{l1}, \code{l2}, \code{l3}, and \code{l4} as applicable.
}
\examples{
mats <- c(1, 3, 6, 12, 24, 60, 120)
NS_loadings(.0609, mats)
NS_loadings_ls(.0609, mats)
Svensson_loadings(12, 48, mats)
}
\seealso{
\code{\link{sc_fit}}, \code{\link{dly_fit}}
}
