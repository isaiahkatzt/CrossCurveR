\name{cr_synthetic}
\alias{cr_synthetic}
\title{Simulate zero-mean covariance regression}
\description{
Generate latent factors, responses and true conditional covariances
for covariance regression examples and numerical checks.
}
\usage{
cr_synthetic(N = 500, B, S0, X = NULL, intercept = TRUE, seed = NULL)
}
\arguments{
\item{N}{Number of observations.}
\item{B}{Numeric \eqn{P\times Q} covariance loading matrix.}
\item{S0}{Positive definite \eqn{P\times P} covariance matrix; simulation
uses its Cholesky factor.}
\item{X}{Optional \eqn{N\times Q} design matrix. If absent, generate standard
normal features, with an optional leading column of ones.}
\item{intercept}{Include a leading column of ones in a generated design.
Ignored when \code{X} is supplied.}
\item{seed}{Optional seed passed to \code{set.seed()}. Simulation advances
the session RNG state; it does not restore the preceding state.}
}
\details{
Responses follow \eqn{w_t=\gamma_t Bx_t+\epsilon_t}, with independent standard
normal latent \eqn{\gamma_t} and Gaussian errors of covariance \eqn{S_0}.
True conditional covariances are \eqn{S_0+(Bx_t)(Bx_t)'}. No parameters are
estimated. Responses and features use the supplied units.
}
\value{
A list with observation-by-response \code{W}, design \code{X}, latent
\code{gamma}, Gaussian errors \code{eps}, supplied \code{B_true} and
\code{S0_true}, and a \eqn{P\times P\times N} array \code{Sigma_true}.
}
\examples{
sim <- cr_synthetic(N = 40, B = matrix(c(.2, .1, -.1, .3), 2),
                    S0 = diag(c(.4, .6)), seed = 12)
dim(sim$W)
dim(sim$Sigma_true)
}
\seealso{
\code{\link{zero_mean_covreg}}, \code{\link{cr_build}}
}
