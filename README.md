# CrossCurveR

Tools for estimating, examining, and stress-testing multiple yield curves.
The core models are Diebold-Li-Yue (DLY) and covariance-regression-based
multi-curve estimation (MCE).

## Development status

Package development is in progress on the `crosscurveR-package` branch.
Reusable code now lives in `R/` and loads through the package namespace.
Reference help is provisional documentation scaffolding for the current public
functions, including units, controls, return values and examples. Each component
will be written and reviewed individually before submission. The documentation
scheme, final API selection, further numerical validation, and publication
examples remain work in progress.

Source files in `R/` retain module prefixes to make their responsibilities visible:

| Prefix | Module |
|---|---|
| `bl_` | Single-curve baselines, cross-curve features, and MCE fitting |
| `cr_` | Covariance regression, covariance construction, and related evaluation |
| `dr_` | Yield-data input, formatting, and reconstruction |
| `st_` | Stress testing and shock propagation |
| `dly_` | Diebold-Li-Yue modelling |

Preserve these prefixes when adding or reorganizing scripts. Shared formatting,
curve-layout, and matrix utilities retain descriptive filenames.

DLY code is organized into `dly_sc.R` (country-curve estimation),
`dly_static.R` (global decomposition and reconstruction), `dly_dynamics.R`
(state estimation and propagation), and `dly_fit.R` (fitting orchestration).
Use `dly_fit()` for a complete fit, `dly_sc_fit_all()` for country fits, and
`dly_fit_from_sc()` to build a complete fit from those country fits.
The existing `dly_all_sc()` and `dly_fit_from_single_curve()` names remain
available as compatibility aliases with the same arguments and results.

Shock code is organized into `st_shock_state.R` (shared states, features and
projection), `st_dly_shock.R` (DLY reconstruction), and `st_fit_shock.R`
(DNS and MCE fitting orchestration). Use `mc_fit_shock_base()` to prepare a cache,
`mc_fit_shock_cached()` to apply a shock, and `mc_fit_shock_refit()` to request
ECM and covariance refits. `mc_fit_shock()` remains the convenience entrypoint.
The public fitting and scenario functions retain their longer compatibility
aliases; see `?shock_fits` for the name mapping.
Scalar shocks and application profiles share `mc_fit_shocked()` for MCE
reconstruction and fitting after factor shocks have been constructed.
Use `sc_fit_shock()` to calibrate independent DNS curves and apply factor
profiles through the same public workflow, optionally reusing its baseline.
Intermediate estimation, matrix-layout and shock-state helpers are internal.
Repository-only numerical diagnostics explicitly access private calculations;
those helpers are not supported package interfaces.

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

## Estimation input contracts

Each yield panel is a data frame with its observation index first, followed by
finite numeric yields. Multi-curve panels must have the same unique, increasing
observation index. Curve names identify countries; `curves` requests their order.
Maturity labels are matched to the positive, unique, increasing `mats` grid.
Equivalent prefixed or unprefixed M/Y labels are accepted. Ordering named curves
or labelled maturities moves their values with them. Fitting never intersects
dates, fills gaps, drops observations, or guesses units. Ordinary market-calendar
gaps are allowed; dynamic steps mean successive observations.

The existing workflow uses maturities in **months** and yields in **percent**
(`5` means 5%). NS decay is in inverse months; Svensson decay arguments are
month-based time constants. Zero and negative yields are valid. Rate definition and compounding
are supplied by the caller: estimation does not turn par yields into zero rates.
Pricing examples interpret inputs as continuously compounded zero yields.

Independent-curve yield residuals (`W`) use **observed minus fitted** for NS and
DLY country fits. Shock reconstruction retains the baseline residual sample.
Residual sign changes can affect signed coefficients, latent diagnostics and
the first EM convergence check; see `?estimation_inputs` for these distinctions.

Regularized EM fits use `CR_tol` for parameter-change thresholds with or without
path storage. The likelihood-and-storage branch also requires the likelihood
score to stabilize, so its full stopping rule can still differ from the unstored
branch. See `?estimation_inputs`.

Covariance fits warn with class `crosscurver_nonconvergence` when the iteration
limit is reached without confirmed convergence. Warnings include the estimator,
iteration budget, tolerance, stopping diagnostics and, for multi-curve fitting,
tenor and response identities. Legacy `term = FALSE` runs explicitly report that
convergence was not assessed. Final estimates and return structures are retained.

Within-curve `H` construction warns with class `crosscurver_h_not_psd` for a
materially negative eigenvalue, including cutoffs, maturities and curve context.
Negligible negative eigenvalues within a relative tolerance are accepted.
The matrix is returned unchanged; no PSD repair is applied.

`Xt` must be a finite numeric matrix with one row per yield observation. Date
row names are checked when supplied; numeric data-frame row labels and matrices
without row names use positional alignment.
Feature order is part of the fitted design. Named covariance initializers and
prediction features are matched by identity. No intercept or normalization is
added. In particular, `mc_fit_exo()` retains its existing behavior: `reference`
and `X_normalize` are inactive compatibility arguments.

Complete fits add `input_spec` metadata describing dates, curve and maturity
order, units, decay choices, and feature layout. Endogenous fits also record the feature basis,
centering, scaling, and smoothing actually used. See `?estimation_inputs` for the full contract.

## Cached MCE scenarios

When `baseline_state` is supplied to `mc_fit_shock()`, omitted fitting controls
inherit the saved configuration. Fresh calls retain the ordinary fitting defaults.
Explicit baseline-only inputs must match the cache.

`reuse_BS0` and `reuse_ECM` independently control conditioning on baseline
covariance-regression coefficients and cointegration estimates. Re-estimating
the ECM while retaining covariance coefficients is supported when the resulting
feature columns can be aligned with the fitted design. Feature normalization,
selection, and smoothing remain scenario controls. With `X_normalize = TRUE`,
spreads are centered and scaled within the scenario sample, preserving the
pre-existing application convention. This sample standardization is separate
from identification of the cointegration basis.

Returned `BS0` and `xt_state` identify the actual coefficient/feature pair.
Baseline-relative centering/scaling and the associated restrictions introduced
in the recent feature-reuse pass have been reverted. See `?feature_reuse`.

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
