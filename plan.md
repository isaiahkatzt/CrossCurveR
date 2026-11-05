# CrossCurveR: package and publication readiness assessment

## Working requirements and verification gates

These requirements govern each subsequent change on `crosscurveR-package`.

- **Preserve functionality and application results.** Refactors must retain supported entrypoints, defaults, return structures, naming/order, estimator variants, and reuse/refit behavior. Numerical outputs must be identical or near-identical under the same data, seeds, settings, and software environment. Methodological changes require separate agreement.
- **Verify every implementation change.** Install the candidate package, run the existing package/namespace and application-helper tests, and compare the affected application outputs against a reference captured before the change. Documentation-only changes require document/diff review.
- **Keep the implementation lean.** Review the diff for duplicated logic, unnecessary wrappers, unused code, extra dependencies, and redundant tests. Retain focused tests of meaningful behavior and previously confirmed defects. Temporary validation artifacts stay outside package source.
- **Keep responsibilities clear.** Estimators fit parameters; covariance builders construct covariance paths; fitting entrypoints coordinate whitening/projection/rescaling; applications express actuarial scenarios and valuations. Review dependency direction and public-interface compatibility before accepting a refactor.
- **Proceed in small steps.** Obtain the user's sign-off before additional structural/organizational changes, as requested. Report verification results and remaining limits before the next step. Commit and push when requested.

The verification pipeline for an implementation change is:

1. **Capture the reference.** Install the currently accepted implementation into an isolated library and record application outputs before editing. Keep input hashes, configuration, convergence diagnostics, and session information with the reference. Preserve relevant local application settings.
2. **Make one scoped change.** Keep algorithm/default changes separate from restructuring so differences can be attributed to a specific cause.
3. **Install and test the candidate.** Run `R CMD INSTALL --library=<candidate-library> .`, then `Rscript --vanilla unit_tests/run_unit_tests.R` with that library first on `R_LIBS`. Resolve new failures or warnings.
4. **Compare application behavior.** Re-execute the same application selection using the candidate package and compare it with the captured reference. The current repeatable subset covers 2015–2016, Observed/DNS/DLY/MCE baselines, the 2015 level shock, fitted valuation curves, return factors, wealth, and country bond profit/loss tables. Supplement it with affected fitting paths and other `applications/` workflows; rerun full affected notebooks for broader changes and before release. A subset pass does not verify unexecuted applications or scenarios.
5. **Review code quality and structure.** Run `git diff --check` and inspect the diff for unnecessary additions, misplaced responsibilities, broken interfaces, and unintended edits. These are explicit review gates; numerical tests alone cannot establish them.
6. **Record the outcome.** Report the reference used, workflows exercised, numerical comparison, test status, code/organization review, and any unverified scope. Accept the change only when all applicable gates pass.

For the application comparison, data hashes, recorded configuration, and iteration/convergence diagnostics must match exactly. Numerical tables and fitted curves use `all.equal(..., tolerance = 1e-10, check.attributes = TRUE)` as a tight comparison threshold, not a blanket definition of actuarial significance. Preserve units and qualitative conclusions; investigate differences and agree any output-specific tolerance before changing the acceptance rule. Do not loosen tolerances or replace the reference merely to make a failing comparison pass.

Application validation uses the installed package and the notebook's actual computational chunks. Preserve existing references when comparing candidates. Standalone validation scripts, reference snapshots, candidate results, and reports remain local and are excluded from commits and package archives. The computational subset does not render plots; inspect representative rendered outputs when presentation is affected.

The assessment below records the initial findings; completed work is checked off as it is verified.

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

Completed on `crosscurveR-package`:

- [x] **Normalise maturity names before reconstruction.** Reject unrecognised or duplicate tenors and preserve downstream `X[TENOR]` names. ([9cf9b38](https://github.com/isaiahkatzt/CrossCurveR/commit/9cf9b38))
- [x] **Fix default covariance-regression initialisation.** Apply `match.arg(init)` in `zero_mean_covreg()` and test the default against explicit adaptive initialisation. ([f2020f1](https://github.com/isaiahkatzt/CrossCurveR/commit/f2020f1))
- [x] **Use a consistent regularised covariance in the E-step and likelihood.** Remove the extra ridge from inverses in both regularised estimators; verify independent likelihood calculations and comparable full tontine application results. ([aee15d9](https://github.com/isaiahkatzt/CrossCurveR/commit/aee15d9))
- [x] **Correct covariance square roots and define the whitening policy.** Support scalar matrices, tolerate negligible asymmetry and negative eigenvalues, and reject material violations. Singular PSD square roots are allowed; whitening requires numerical invertibility. Positive eigenvalues remain unchanged. ([babd366](https://github.com/isaiahkatzt/CrossCurveR/commit/babd366))
- [x] **Validate EM iteration controls.** Require positive whole-number iteration counts and checking intervals, and finite positive convergence tolerances, in the regularised EM estimator and comparator. Confirm informative errors and unchanged valid-input results across fitting paths and a tontine notebook subset. ([8d1289f](https://github.com/isaiahkatzt/CrossCurveR/commit/8d1289f))

### Consolidate implementations and research material

- [x] Move the production EM estimator from `R/cr_eval.R` into `R/cr_core_estim.R`, preserving its implementation and interface.
- [x] Share private numerical helpers for inversion, diagonal shrinkage, covariance regularisation, and likelihood evaluation within `R/cr_core_estim.R`; preserve estimator-specific loops, defaults, and interfaces.
- [x] Review tenor-fitting call paths and retain `blcc_fe()` in `R/bl_feature_extract.R`. Share the pointwise covariance calculation across reconstruction, tenor paths, and simulation while preserving output shapes, sparse assembly, tenor order, and random draws.
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
