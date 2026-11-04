# CrossCurveR: package and publication readiness assessment

## 1. Current position

**CrossCurveR has the computational foundation for a publishable R package. The main work is to establish a supported interface, strengthen numerical validation, and make the research reproducible outside your development environment.**

Following your clarification, the package should centre on **multi-curve modelling through both MCE and DLY**. Reconstruction and stress testing support those models; tontine and Ogden analyses serve as reproducible actuarial applications.

I reviewed the current working tree, including local changes, and ran the existing tests plus additional numerical and integration checks. No files were changed during that assessment.

| Area | Current state |
|---|---|
| Implementation | Distinct modules for reconstruction, single-curve fitting, MCE, DLY, and stress propagation |
| Applications | Substantial tontine and Ogden implementations and notebooks |
| Tests | 14 test files; 92 test blocks; **431 assertions passed, one skipped**, no test failures |
| Packaging | No `DESCRIPTION`, `NAMESPACE`, conventional `R/` directory, package help pages, or installation guide |
| Reproducibility | A lockfile exists, but execution depends on repository paths, script sourcing, and incompletely documented data |
| Verification limits | Tests ran locally on R 4.2.3; I did not measure coverage or fully render every research notebook |

AAS is a plausible venue: its Actuarial Software stream permits existing methods and new methodology, but expects a substantive software contribution with worked actuarial applications. It specifies GPL-2 or GPL-3 licensing, component tests summarised under **“Unit Testing and Validation”**, replication materials, and DOI-backed archival distribution. **CRAN submission is optional.** [AAS submission guidance](https://www.cambridge.org/core/journals/annals-of-actuarial-science/information/author-instructions/preparing-your-materials#actuarial-software)

## 2. What needs adding, removing, or cleaning up

### Establish the package structure

Convert the repository root into a conventional R package:

- Add `DESCRIPTION`: purpose, authors, maintainer, version, supported R version, licence, dependencies, and project/issue URLs.
- Move reusable implementation into `R/`; generate `NAMESPACE` and `man/` documentation using roxygen2.
- Move package tests into `tests/testthat/` and add the standard test runner.
- Add a README with installation instructions and a short, executable multi-curve example.
- Add a quick-start vignette, model-methodology documentation, `NEWS.md`, and citation metadata.
- Add `.Rbuildignore` so research outputs, development environments, and full replication materials stay outside the installation archive.

These are conventional R packaging mechanisms and apply equally to a GitHub-distributed package. [Writing R Extensions](https://stat.ethz.ch/R-manual/R-devel/doc/manual/R-exts.html#Package-structure)

### Remove dependence on a preconfigured session

The current bootstrap attaches 24 packages and installs missing ones automatically. Tests source implementation files into `.GlobalEnv`, and notebooks repeat similar loading sequences.

Replace that arrangement with:

- Explicit dependency declarations and qualified calls such as `Matrix::bdiag()`.
- A small runtime dependency set; plotting, benchmarking, and manuscript tools remain development or replication dependencies.
- Package loading through `library(CrossCurveR)`.
- Tests that exercise the package namespace.
- Local restoration of temporary environment changes, rather than persistent changes to language settings or working directories.

Keep `renv` for reproducing the paper’s environment. Package users should obtain dependencies through normal package installation.

### Fix confirmed defects before freezing a release

| Finding | Required change |
|---|---|
| Reading a supplied CSV with ordinary `read.csv()` introduces `X` prefixes in maturity names. Passing that result to `dr_build_curve()` returned **all missing yields** in my probe. | Normalise maturity names at input and reject unrecognised names before reconstruction. See [reconstruction code](data_formatting/dr_formatting.R), line 120. |
| `zero_mean_covreg()` fails with its default `init` because its two-element default reaches a scalar `if`. Existing tests supply `init` explicitly. | Apply `match.arg()` and test default calls. See [estimator code](yield_estimation/cr_core_estim.R), line 5. |
| The regularised EM likelihood uses an additionally ridged inverse but a determinant without that additional ridge. An independent calculation disagreed with the reported likelihood. | Use one consistently defined covariance in estimation, likelihood evaluation, and diagnostics. See [likelihood calculation](yield_estimation/cr_eval.R), line 55. |
| `symmPSD_sqrt(matrix(4, 1, 1))` fails; inverse square roots fail for singular matrices. All negative eigenvalues are currently clipped indiscriminately. | Handle scalar matrices; distinguish rounding error from materially invalid covariance matrices; document the singular-matrix policy. |
| Zero iterations can produce an undefined-variable error in the EM estimator. | Validate iteration controls and return informative input errors. |

### Consolidate implementations and research material

- Move the production EM estimator out of the evaluation-labelled file.
- Retain one authoritative implementation for each deliberately supported estimator variant. The fixed-covariance comparator should be explicitly identified as such.
- Remove unfinished placeholder functions from the supported release; preserve future work in issues or development branches.
- Consolidate duplicated fitting and rescaling operations used by ordinary and shocked MCE fits.
- **Audit the role of the transformed residual covariance `H` against the intended mathematics.** It is computed but does not enter the final projection in the current MCE paths. Changing maturity cutoffs produced identical fitted yields in my probe. Resolve the theoretical intent before deleting these calculations.
- Preserve exploratory notebooks in the research archive, while selecting and repairing a smaller publication replication set. Existing examples include missing raw-data dependencies, obsolete `overlap` arguments, and references to objects whose construction is commented out.

## 3. Additional testing needed

**The largest gap is independent numerical and integration testing.** The current suite includes useful mathematical checks, particularly for NS reconstruction and DLY decomposition, but several top-level MCE tests replace almost the entire computation with mocks.

| Priority | Test area | Specific additions |
|---|---|---|
| Essential | Installed-package operation | Install the built archive into a temporary library; run examples from an unrelated directory and clean R session without sourcing files or attaching undeclared dependencies. |
| Essential | Input contracts | Misaligned dates, reordered curves, duplicate/missing dates, inconsistent maturity columns, `X` prefixes, nonnumeric yields, missing values, short samples, and incompatible exogenous designs. |
| Essential | Covariance regression | Default arguments, invalid controls, singular designs, covariance symmetry/positive definiteness, independent likelihood calculations, shrinkage boundaries, convergence, and iteration exhaustion. |
| Essential | MCE integration | Genuine endogenous and exogenous fits on small fixtures; correct ordering through permutation, whitening, fitting, and rescaling; explicit no-cointegration behaviour. |
| Essential | DLY integration | Independently calculated PCA/OLS examples, factor orientation, dynamic transitions, rank-deficient inputs, and supported shock propagation. Retain the existing decomposition checks. |
| Essential | Stress behaviour | Zero-shock identity, targeted factor changes, window boundaries, reuse versus refit behaviour, baseline immutability, and incompatible cached states. Exercise both model families. |
| Supporting | Reconstruction | Observed-value preservation, interpolation boundaries, insufficient observations, zero/negative rates, and parameter/unit consistency. Log-yield interpolation needs a clear positive-yield restriction. |
| Supporting | Replication applications | Known-answer discounting and cashflow examples, tontine wealth conservation, and Ogden asset/liability identities for the applications included in the paper. |

Additional testing principles:

- Compare covariance-regression recovery using fitted covariance matrices and sign-invariant quantities, because `B` and `-B` imply the same covariance.
- Test ordinary EM likelihood ascent under suitable unregularised conditions; validate regularised fitting against its specified objective.
- Separate fast deterministic tests from longer simulation studies and benchmarks.
- Ensure disabling path storage does not inadvertently change the convergence rule.
- Use fixed seeds and numerical tolerances; avoid exact floating-point snapshots.
- Measure coverage to locate untested branches, without treating a percentage as evidence of statistical correctness.

My additional checks successfully exercised DLY, exogenous MCE, and a four-curve MCE baseline with a zero shock. A separate two-curve endogenous example found no usable cointegration relationships. That is an important supported outcome to document and test.

## 4. How the project’s formulation should change

### Present one multi-curve toolbox

A suitable organising description is:

> CrossCurveR provides tools for estimating, examining, and stress-testing multiple yield curves through DLY and covariance-regression-based MCE models, with reproducible actuarial applications.

Build the user workflow around:

**Prepare aligned curves → fit MCE or DLY → inspect fitted curves and diagnostics → apply shocks → use outputs in an actuarial analysis.**

Retain the established fitting entrypoints initially. Add documented model classes and consistent `print()`, `summary()`, and `fitted()` access, while preserving model-specific components. Expose supported DLY shock operations alongside MCE stress functionality; keep matrix manipulation and orchestration helpers internal.

Standardise and document:

- Maturities in months and decay-parameter units.
- Yield units and compounding assumptions.
- Date alignment and missing-data requirements.
- Fitted curves, reconstructed curves, and dynamically propagated paths as distinct outputs.
- Convergence status and numerical adjustments.
- Which quantities remain fixed or are re-estimated under each shock operation.

### Connect implementation to methodological claims

For DLY, document precisely the implemented PCA/OLS factor extraction and dynamic specification, including differences from the cited formulation.

For MCE, explain the covariance model, endogenous feature construction, identification, regularisation, optimisation, and final projection. Justify defaults and expose nonconvergence clearly.

The paper should distinguish:

- Numerical implementation correctness.
- Statistical estimation performance.
- Practical usefulness of joint curve modelling and stress transmission.

Add simulation recovery studies and sensitivity analysis. Evaluate reconstruction using held-out observed points where possible. If claiming predictive advantages, add time-separated evaluation with preprocessing and tuning confined to training data. Current in-sample comparisons alone cannot establish those advantages.

## 5. Data, release sequence, and acceptance criteria

Keep package source and paper replication in the same repository, but separate their execution requirements.

- Bundle a small documented synthetic example for installation, tests, and tutorials.
- Archive full empirical data separately with provenance, retrieval dates, units, reconstruction settings, and redistribution terms.
- Provide scripts linking source data to the exact analysis panels.
- Review local data variants before including them: several contain missing timestamps. Only seven reconstructed CSVs are currently tracked.
- Convert the selected applications to consume the installed package.
- Provide one replication command, recorded seeds, session information, expected outputs, and runtime guidance.

Complete the transition in this order:

1. **Numerical audit:** resolve the confirmed defects and estimator/theory inconsistencies.
2. **Package conversion:** establish dependencies, namespace, documented interfaces, and conventional tests.
3. **Validation:** add independent numerical tests, genuine integration tests, and continuous checks on Linux, macOS, and Windows.
4. **Research replication:** curate the data workflow, synthetic studies, and actuarial applications.
5. **GitHub release:** build an installable source archive, tag the version, and archive that exact release with a DOI.

The first release is ready when a new user can install it, run both model examples, inspect diagnostics, and reproduce the selected application without knowing the source tree. Require clean `R CMD check` results, no unexplained test failures or skips, and reproducible manuscript outputs.

**Chosen defaults:** GitHub distribution first; both MCE and DLY are core; actuarial applications remain replication examples; GPL-3 for the AAS software target; CRAN submission deferred.
