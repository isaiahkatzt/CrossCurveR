\name{cr_build}
\alias{cr_build}
\title{Construct conditional covariance matrices}
\description{
Evaluate fitted covariance regression loadings and intercept covariance
at supplied feature rows.
}
\usage{
cr_build(object, X_new = NULL, return_array = TRUE, drop_single = TRUE)
}
\arguments{
\item{object}{A fitted covariance regression object or a list containing
\code{B} (response-by-feature) and \code{S0} (response-by-response).}
\item{X_new}{Finite numeric matrix with one feature row per requested
covariance and one column per column of \code{B}. If omitted,
\code{object$X} is used when present. The package estimators do not store
\code{X}, so normally supply it explicitly.}
\item{return_array}{Return an array when \code{TRUE}, or a list when \code{FALSE}.}
\item{drop_single}{For a one-row request in list mode, return the sole matrix
instead of a one-element list. Array mode retains all three dimensions,
including when this argument is \code{TRUE}.}
}
\details{
Each feature row gives \eqn{\Sigma_t=S_0+(Bx_t)(Bx_t)'}. This function does
not estimate parameters, center features, add an intercept, or regularize
the supplied covariance. Named feature columns are matched to saved feature
identities, falling back to \code{colnames(B)} when no saved specification
is available. Unnamed designs use positional columns.
}
\value{
By default, a numeric \eqn{P\times P\times N} array in feature-row order.
For \code{return_array = FALSE}, a length-\eqn{N} list of covariance matrices,
or a single matrix if \eqn{N=1} and \code{drop_single = TRUE}.
}
\examples{
model <- list(B = matrix(c(.2, -.1), 2, 1), S0 = diag(c(.4, .6)))
cr_build(model, X_new = matrix(c(-1, 0, 1), ncol = 1))
}
\seealso{
\code{\link{zero_mean_covreg}}, \code{\link{cr_synthetic}}
}
