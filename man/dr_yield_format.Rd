\name{dr_yield_format}
\alias{dr_yield_format}
\alias{dr_daily}
\title{Format raw yield panels and interpolate isolated missing observations}
\description{
Assign maturity labels and a date window to a raw panel, or fill eligible
isolated missing values using neighboring observations.
}
\usage{
dr_yield_format(data = "", mats = c(), start_date = "", end_date = "")

dr_daily(data, reference = NULL, tkey = "time", threshold = 2)
}
\arguments{
\item{data}{A data frame. Formatting assumes the timestamp is first and
assigns maturity labels positionally. Daily interpolation expects a column
named by \code{tkey}, with the remaining columns coercible to numeric yields.}
\item{mats}{Maturities in months, one per raw yield column, in supplied order.}
\item{start_date, end_date}{Date-compatible window boundaries. Supply both:
the historical empty-string defaults do not specify a valid date window.}
\item{reference}{Optional ordered vector of dates defining the interpolation
rows. Existing observations are matched by date; unmatched rows start missing.
When absent, retain the existing rows.}
\item{tkey}{Observation-index column name for \code{dr_daily()}.}
\item{threshold}{Maximum count of observed yields in a row for that row's
missing entries to be eligible for interpolation. This counts observed entries,
not missing entries.}
}
\details{
\code{dr_yield_format()} parses supported timestamps in UTC, assigns
unprefixed maturity names, filters between the two boundaries, returns the
time column as ISO date strings, and converts NaN entries to NA. Bounds are
converted to UTC midnight, so an intraday observation after midnight on
\code{end_date} is outside the current filter. Unlike numeric observation
indices in the fitting interfaces, numeric timestamps here represent seconds
since the Unix epoch. This function does not infer the original maturity
identities or interpolate values.

\code{dr_daily()} replaces an eligible missing entry by the arithmetic mean
of the preceding and succeeding row's original yields at that maturity.
Both neighbors must be available. The operation is one pass: a run of missing
entries and endpoint gaps are not filled through recursive interpolation.
Rows correspond to the supplied reference dates or existing observation order;
no trading calendar is inferred. Existing nonmissing entries are retained.
Remaining missing data must be handled before estimation.
}
\value{
Formatting returns a data frame with a \code{time} column first and unprefixed
M/Y yield labels. Interpolation retains its input column order, yield labels
and \code{tkey}.
}
\examples{
raw <- data.frame(date = c("2020-01-01", "2020-01-02", "2020-01-03"),
                  short = c(1, NA, 3), long = c(2, NA, 4))
panel <- dr_yield_format(raw, c(6, 12), "2020-01-01", "2020-01-03")
dr_daily(panel, threshold = 0)
}
\seealso{
\code{\link{dr_read_yields}}, \code{\link{dr_build_curve}}, \code{\link{estimation_inputs}}
}
