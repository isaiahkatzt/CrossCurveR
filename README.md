# CrossCurveR

Tools for estimating, examining, and stress-testing multiple yield curves.
The core models are Diebold-Li-Yue (DLY) and covariance-regression-based
multi-curve estimation (MCE).

## Development status

Package development is in progress on the `crosscurveR-package` branch.
Reusable code now lives in `R/` and loads through the package namespace.
API documentation, numerical validation, and publication examples remain
work in progress.

## Installation

From R, install the development branch and its runtime dependencies:

```r
# Run once if remotes is not already installed:
install.packages("remotes")
remotes::install_github("isaiahkatzt/CrossCurveR", ref = "crosscurveR-package")
library(CrossCurveR)
```

For a local checkout with the declared dependencies already installed:

```sh
R CMD INSTALL .
```

Package loading never installs or attaches third-party packages. The contributed
runtime dependencies are `dplyr`, `Matrix`, `purrr`, `rlang`, `tsDyn`, `vars`, and
`YieldCurve`, plus R's `splines` and `stats` packages. R >= 4.1.0 is required for
native pipes. The dependency packages may impose additional R version requirements.
Plotting, notebook, application, and test tools are declared separately in
`Suggests`; install them explicitly when using those workflows.

## Core entrypoints

```r
library(CrossCurveR)

# Call with an aligned yield panel and the appropriate model controls:
# dly_fit(...)
# mc_fit_end(...)
# mc_fit_exo(...)
# mc_fit_shock_baseline(...)
# mc_fit_shock_from_baseline(...)
```

TODO: add a small reproducible example covering aligned curve inputs,
MCE and DLY fitting, diagnostics, and stress testing.

## Yield data input

Read labelled CSV panels through the validated reader:

```r
yields <- dr_read_yields("curve.csv")
rebuilt <- dr_build_curve(yields, nmats = c(1, 3, 6, 12, 24, 60))
```

Reads retain the downstream `X01M`, `X01Y` naming convention. The reader accepts
prefixed or unprefixed M/Y labels, normalizes aliases such as `12M` to `X01Y`,
and rejects unknown or duplicate tenors. `dr_build_curve()` accepts either style
and keeps its existing unprefixed output names. Date and yield values are preserved.

For an existing data frame, use `dr_normalize_maturity_names()`. For raw files
with arbitrary rate-column headers, retain the original CSV names with
`check.names = FALSE` and map columns explicitly through `dr_yield_format(mats = ...)`.

## Development and tests

From R in the checkout, explicitly install the declared development and test
dependencies, including the supported `testthat` and `waldo` versions:

```r
remotes::install_deps(dependencies = TRUE)
```

Then install the current checkout before exercising its installed namespace:

```sh
R CMD INSTALL .
Rscript --vanilla unit_tests/run_unit_tests.R
```

Core package tests live in `tests/testthat/` and also run under `R CMD check`.
The repository runner additionally checks application helpers in `unit_tests/`.
Tests include a separate R session with no default packages attached, running
outside the repository. Test-specific locale changes are scoped to the run.

Refactors also require before/after application comparisons and reviews of code
size and organization. Standalone validation scripts, saved references, and
reports remain local.

Research notebooks load `CrossCurveR` and their analysis packages explicitly.
Application helpers are sourced into the notebook's own environment.
The package does not activate `renv`, change language settings, select a local
Pandoc installation, or change the working directory at startup. The existing
`renv.lock` remains available for an explicit `renv::restore()` when reproducing
the archived research environment.

## Models and applications

- **MCE:** TODO: document estimation, controls, and outputs.
- **DLY:** TODO: document factor extraction, dynamics, and outputs.
- **Actuarial applications:** TODO: link reproducible tontine and Ogden examples.

## Citation

Citation metadata is provided in [CITATION.cff](CITATION.cff) for GitHub
and [inst/CITATION](inst/CITATION) for R:

```r
citation("CrossCurveR")
```

The maintainer identity follows the existing Git author record. The publication
author list and release DOI will be reviewed before publication.

## Licence

GNU General Public License, version 3 (GPL-3).
