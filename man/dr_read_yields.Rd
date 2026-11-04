\name{dr_read_yields}
\alias{dr_read_yields}
\alias{dr_normalize_maturity_names}
\title{Read Yield Panels and Normalize Maturity Names}
\description{
Read a CSV yield panel or normalize the maturity names of an existing data frame.
}
\usage{
dr_read_yields(file, tkey = "time")
dr_normalize_maturity_names(data, tkey = "time")
}
\arguments{
  \item{file}{A CSV file path or connection accepted by \code{utils::read.csv()}.}
  \item{data}{A data frame with a time column and maturity columns.}
  \item{tkey}{The input time column name.}
}
\details{
Labels accept an optional X prefix, a positive number, and an M (months) or
Y (years) suffix. Case and surrounding whitespace are normalized. Unknown
labels and duplicate tenors, including aliases such as \code{12M} and
\code{01Y}, cause an error before reconstruction.

The time column is placed first, with maturity columns in their original order.
Values are unchanged. The CSV reader names its time column \code{time}; the
normalizer retains \code{tkey}. Both return \code{X01M}/\code{X01Y} labels for
compatibility with downstream workflows. \code{dr_build_curve()} accepts either
prefixed or unprefixed input and retains its unprefixed output convention.
}
\value{A data frame with normalized maturity names.}
\examples{
panel <- data.frame(time = "2020-01-01", a = 1, b = 2)
names(panel)[-1] <- c(" x6m ", "12M")
dr_normalize_maturity_names(panel)
dr_build_curve(panel, nmats = c(6, 9, 12))
}
