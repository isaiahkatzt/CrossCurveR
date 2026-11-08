\name{estimation_inputs}
\alias{estimation_inputs}
\title{Alignment and units for estimation inputs}
\description{
Shared input conventions for single-curve, DLY, MCE, MCX and reusable
stress-baseline estimation. These checks do not impute observations, change
the estimation sample, normalize covariates, or convert rate conventions.
}
\details{
Multi-curve inputs are lists of data frames. Each panel has an observation index
first and finite numeric yields thereafter. Names identify curves; \code{curves}
specifies their order. Lists with no names can be labelled positionally by an
explicit \code{curves} argument. Partial, duplicate or conflicting names fail.
Named panels are ordered by identity, never relabelled as another country.

Observation indices must be nonmissing, unique and strictly increasing, and
must match across curves. Supported indices are \code{Date}, \code{POSIXct},
ISO date strings and numeric observation indices. Numeric indices are not
implicitly interpreted as Unix or spreadsheet dates. Date representations are
compared without replacing the supplied index in outputs. Market-calendar gaps
are allowed. Dynamic lags mean successive observations, not fixed calendar days.
Missing or mismatched dates must be resolved explicitly before fitting.

Numeric \code{mats} are positive, unique, increasing maturities in months.
Each labelled panel must contain exactly those maturities. Equivalent labels,
including \code{12M}, \code{01Y} and \code{X01Y}, identify the same maturity;
duplicate economic tenors fail. Labelled columns may be reordered to match
\code{mats}, retaining their labels. The public \code{matname_to_numeric()}
continues to return years, while \code{numeric_to_matname()} takes months.
Unlabelled matrices accepted by lower-level helpers use positional columns.

The package workflow convention is percentage yields: \code{5} means five
percent. Zero and negative yields are allowed. Country factors and yield
residuals retain that scale; yield covariances use squared percentage points.
There is no magnitude-based unit inference or automatic rescaling. Standalone
covariance regression uses response and feature units as supplied.

NS decay parameters multiply maturity and are in inverse months. Svensson
parameters named \code{lambda1}, \code{lambda2} or supplied in \code{slambdas}
are time constants in months. Both must be finite and positive. Changing
maturities from months to years requires multiplying NS decay by twelve and
dividing Svensson time constants by twelve. Fixed decay inputs are validated
just like grids. \code{cutoffs} are maturity-column indices, not durations.

No estimation function converts par yields, yields to maturity or compounding
conventions. The application pricing examples interpret rates as continuously
compounded zero yields, convert percentage rates to decimals and months to
years, then discount. Users must establish that their source data support
that interpretation.

\code{Xt} contains one finite numeric row per yield observation. Nondefault
date row names allow the row index to be checked; numeric data-frame row labels
are positional even after subsetting. Unindexed matrices are
explicitly positional. Feature names, when supplied, must be unique. Covariance
prediction orders named features against the fitted design and rejects
conflicting names. Initial \code{B} and \code{S0} matrices are dimension-checked
and matched to named response and feature identities where available.
Intercepts and normalization are never inserted automatically. The existing
\code{mc_fit_exo()} arguments \code{reference} and \code{X_normalize} remain
inactive compatibility arguments.

Complete fits record an \code{input_spec} containing the observation index,
curve and maturity order, declared units, decay choices and feature layout. Endogenous fitting
also records the feature bases, centering, scaling and smoothing specification.
These describe transformations actually used; they do not change scenario
normalization or make re-estimated bases interchangeable with frozen coefficients.
Precomputed DLY inputs must satisfy the same alignment and matrix-layout checks.
}
\seealso{\code{\link{dly_fit}}, \code{\link{dr_read_yields}}}
