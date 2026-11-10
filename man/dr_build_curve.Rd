\name{dr_build_curve}
\alias{dr_build_curve}
\title{Reconstruct yield curves across maturities}
\description{
Fill missing maturity yields by spline, log-yield interpolation, or
parametric NS/Svensson fitting for each observation.
}
\usage{
dr_build_curve(data, nmats = NULL, tkey = "time", build = c("spline",
    "log", "parametric"), par_method = c("nelson", "svensson"),
    spline_method = c("fmm", "natural", "hermite", "periodic"),
    complete = FALSE, precision = 4)
}
\arguments{
\item{data}{Data frame with an observation-index column and maturity-labelled
yield columns. Missing yields are permitted during reconstruction.}
\item{nmats}{Additional maturities in months. The output uses the sorted union
of input and requested maturities; this is not a column-selection argument.}
\item{tkey}{Observation-index column name.}
\item{build}{Interpolation or parametric reconstruction method.}
\item{par_method}{Parametric NS or Svensson calibration through
\pkg{YieldCurve}, used only for \code{build = "parametric"}.}
\item{spline_method}{\code{stats::splinefun()} method. \code{"hermite"} maps
to \code{"monoH.FC"}. Used only for spline construction.}
\item{complete}{For parametric rows containing missing values, replace all
maturity yields by the fitted curve when \code{TRUE}; otherwise preserve
observed values and fill gaps. Rows without missing values are returned
unchanged before rounding. This argument is inactive for spline/log methods.}
\item{precision}{Number of decimal places used to round output yields.}
}
\details{
Maturity labels accept M/Y suffixes with an optional X prefix. Unknown or
duplicate economic tenors are rejected. Reconstruction retains observation
order and does not add time points.

Spline construction interpolates observed yields using the selected spline
method, including that method's extrapolation behavior. It can overshoot
observed yields. Log construction linearly interpolates the logarithm of
observed yields, exponentiates the result, and uses constant extrapolation.
Finite log construction therefore requires positive observed yields. This
branch operates on log yields, not log discount factors.

Parametric construction estimates parameters from the available maturity
observations. Nonlinear calibration still requires sufficient informative
tenors for the selected model. Rows with fewer than two available values
produce a warning and are returned without a successful reconstruction.
No solver substitution or automatic regularization is performed here.
}
\value{
A data frame with the supplied index first and reconstructed yield columns
in increasing maturity order, rounded to \code{precision} digits. Maturity
names use unprefixed canonical M/Y labels.
}
\examples{
panel <- data.frame(time = as.Date(c("2020-01-01", "2020-01-02")),
                    X06M = c(1, 1.1), X01Y = c(2, 2.1))
dr_build_curve(panel, nmats = c(6, 9, 12), spline_method = "natural")
}
\seealso{
\code{\link{dr_read_yields}}, \code{\link{dr_yield_format}}, \code{\link{sc_fit}}
}
