\name{dly_fit}
\alias{dly_fit}
\alias{dly_sc_fit_all}
\alias{dly_all_sc}
\alias{dly_fit_from_sc}
\alias{dly_fit_from_single_curve}
\title{Fit a PCA-based Diebold--Li--Yue multi-curve model}
\description{
Estimate country level and slope factors, decompose them using global principal
components, and estimate and recursively propagate their dynamic states.
}
\usage{
dly_fit(yields, mats,
        lambdas = seq(from = 0.001, to = 1, length.out = 100),
        curves = names(yields), reference = curves[[1]],
        center = TRUE, scale. = FALSE)

dly_sc_fit_all(lambdas, yields, mats, curves = names(yields))

dly_fit_from_sc(dly_sc, actual_yields = NULL,
                curves = names(dly_sc), reference = curves[[1]],
                center = TRUE, scale. = FALSE, mats = NULL)
}
\arguments{
  \item{yields}{A named list of country yield data frames. Each has a time
  column first, followed by finite numeric yields in the order of \code{mats}.
  Country observations and maturities must already be aligned.}
  \item{mats}{Numeric maturities, conventionally in months. When omitted from
  \code{dly_fit_from_sc()}, maturities are inferred from the yield column names.}
  \item{lambdas}{Candidate positive decay parameters in inverse maturity units.
  A single value fixes the decay; otherwise each country minimizes its
  two-factor yield sum of squared errors over the supplied grid.}
  \item{curves}{Country names in the same order as the yield list.}
  \item{reference}{Country whose first principal-component loadings determine
  the positive orientation of the global factors.}
  \item{center, scale.}{Centering and scaling controls passed to
  \code{stats::prcomp()}.}
  \item{dly_sc}{Country fits returned by \code{dly_sc_fit_all()} or its
  compatibility alias \code{dly_all_sc()}.}
  \item{actual_yields}{Aligned observed yield data for residual calculation.
  If \code{NULL}, residuals are measured against the reconstructed country fits.}
}
\details{
Country curves use the Nelson--Siegel level and slope loadings. Separate first
principal components of the country level and slope panels provide the global
factors. Country factors are regressed on their matching global factor with an
intercept. Global states follow an unrestricted VAR(1) with an intercept;
country residual states additionally depend on the lagged global states.

The \code{raw} yield component contains independent country fits.
\code{global} includes country intercepts, and
\code{global_plus_idiosyncratic} reconstructs \code{raw}.
\code{dislocation_propagated} is a deterministic recursive path initialized at
the first extracted state, with subsequent innovations set to zero. Transition
parameters are estimated from the complete input sample.

\code{dly_all_sc()} is an alias for \code{dly_sc_fit_all()}, and
\code{dly_fit_from_single_curve()} is an alias for \code{dly_fit_from_sc()}.
The aliases preserve the original arguments, defaults, and results.
}
\value{
\code{dly_sc_fit_all()} returns a named list of country fits containing the
selected decay, time index, yield names, loadings, factors, and fitted yields.

\code{dly_fit()} and \code{dly_fit_from_sc()} return lists containing metadata,
\code{single_curve}, \code{country_factors}, \code{global_factors},
\code{factor_models}, \code{country_betas}, \code{yields}, \code{residuals},
and \code{dynamics}. \code{dly_fit()} also records \code{lambda_grid}.
}
