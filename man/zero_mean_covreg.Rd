\name{zero_mean_covreg}
\alias{zero_mean_covreg}
\alias{zero_mean_covreg_em}
\title{Estimate zero-mean covariance regression by EM}
\description{
Estimate feature-dependent response covariance using either the legacy
zero-mean iteration or the regularized EM estimator.
}
\usage{
zero_mean_covreg(W, X, init = c("adaptive", "static"), max_iter = 1000,
    tol = 1e-10, S0 = NULL, B = NULL, S0_shrink_diag = 0, verb = FALSE,
    term = FALSE)

zero_mean_covreg_em(W, X, init = c("adaptive", "static"), max_iter = 500,
    tol = 1e-08, S0 = NULL, B = NULL, ridge_S0 = 1e-08, ridge_X = 1e-10,
    S0_shrink_diag = 0, check_every = 1, use_loglik = TRUE, verbose = FALSE,
    store_path = TRUE)
}
\arguments{
\item{W}{Finite numeric response matrix with \eqn{N} observations and at
least two response columns. The conditional mean is assumed zero; responses
are not centered automatically.}
\item{X}{Finite numeric \eqn{N\times Q} feature matrix in the same observation
order. Include an intercept column explicitly if desired. Meaningful row
indices must match \code{W}; regularized EM requires full column rank.}
\item{init}{Adaptive randomized latent initialization or small random
\code{"static"} loadings. Used when \code{B} is absent. Set a seed externally
for reproducibility; neither estimator sets one.}
\item{max_iter}{Maximum number of EM updates.}
\item{tol}{Positive convergence tolerance. The two estimators use different
stopping rules as described below.}
\item{S0}{Optional symmetric positive semidefinite \eqn{P\times P} covariance
initialization, where \eqn{P=\operatorname{ncol}(W)}. Named dimensions are aligned
to response identities. The legacy estimator requires an invertible covariance.}
\item{B}{Optional \eqn{P\times Q} loading initialization. Named rows and columns
are aligned to the response and feature identities.}
\item{ridge_S0}{Nonnegative covariance ridge for regularized EM. The historical
additional inversion ridge and scoring convention are retained.}
\item{ridge_X}{Nonnegative regression ridge for regularized EM. It does not
bypass the initial rejection of rank-deficient \code{X}.}
\item{S0_shrink_diag}{Diagonal shrinkage fraction in \eqn{[0,1]}. Off-diagonal
covariance entries are multiplied by one minus this fraction at each update.}
\item{verb}{Legacy interface argument; currently has no effect on output.}
\item{term}{Whether the legacy estimator evaluates its historical convergence
criterion. With \code{FALSE}, it performs the full iteration budget.}
\item{check_every}{Positive integer interval between regularized convergence
checks and scheduled verbose diagnostics.}
\item{use_loglik}{Whether regularized EM computes the historical likelihood
score. A likelihood-change stopping check is additionally used when
\code{store_path = TRUE}.}
\item{verbose}{Print regularized EM diagnostics.}
\item{store_path}{Store regularized parameter paths, and the score path when
\code{use_loglik = TRUE}. See the stopping-rule distinction below.}
}
\details{
For a feature vector \eqn{x_t}, the model is
\deqn{w_t=\gamma_t Bx_t+\epsilon_t,\qquad
\gamma_t\sim N(0,1),\quad\epsilon_t\sim N(0,S_0),}
with conditional covariance
\deqn{\Sigma_t=S_0+(Bx_t)(Bx_t)' .}
No mean regression is fitted. Simultaneously reversing the sign of \eqn{B} leaves covariance
unchanged. Regularized EM normalizes its largest-magnitude loading to be
nonnegative after updates; the legacy estimator retains its signed path.

Regularized EM checks Frobenius changes in \code{B} and \code{S0} against
\code{tol}. With both likelihood computation and path storage enabled, it
also requires the absolute score change to be less than \code{tol}.
Without stored paths, convergence is based on parameter changes even when
the score is computed. Storage therefore retains an additional stopping
condition; it is not purely a memory control in the current implementation.

The legacy estimator has no ridge controls or likelihood score. If
\code{term = TRUE}, every ten iterations it compares stored history entries
\code{iter} and \code{iter - 9}, rather than the newest update at
\code{iter + 1}. It retains that criterion between checks.
These estimator-specific algorithms and defaults are intentionally preserved.

Exhausting the iteration budget without confirmed convergence emits a
\code{crosscurver_nonconvergence} warning, including when legacy checks are
disabled. Inspect \code{converged} and \code{iter} rather than treating a
returned object as evidence of convergence.

Returned latent moments come from the last E-step, before the final M-step.
They are not recomputed under the returned final parameters, and regularized
sign normalization can also change their orientation relative to final
\code{B}. See \code{\link{estimation_inputs}} for precise timing and paths.
}
\value{
A list with final \code{B}, \code{S0}, last-E-step moments \code{mz},
\code{vz}, \code{sz}, \code{iter}, \code{converged}, \code{loglik},
\code{S0_shrink_diag}, and \code{input_spec}.
Neither estimator stores the original \code{X}.

Regularized EM returns \code{B_path} and \code{S0_path} when requested, plus
\code{loglik_path} when scores are enabled. Paths start at initialization and
contain \code{iter + 1} entries. The legacy estimator always returns these
parameter paths, the original \code{b_est}/\code{s_est} history lists,
\code{loglik = NULL}, and \code{termination_reason} equal to
\code{"converged"}, \code{"max_iter"}, or \code{"not_assessed"}.
}
\examples{
set.seed(12)
X <- cbind(intercept = 1, driver = seq(-1, 1, length.out = 100))
B <- matrix(c(.2, .1, .15, -.1), 2)
sim <- cr_synthetic(N = nrow(X), B = B, S0 = diag(c(.4, .6)), X = X)
fit <- zero_mean_covreg_em(sim$W, X, B = B, tol = 1e-5,
                           max_iter = 500, use_loglik = FALSE, store_path = FALSE)
fit$converged
Sigma <- cr_build(fit, X_new = X)
dim(Sigma)
}
\seealso{
\code{\link{cr_build}}, \code{\link{cr_synthetic}}, \code{\link{estimation_inputs}}
}
