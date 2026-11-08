# CrossCurveR (development version)

- Split DLY estimation, static decomposition, dynamics, and fitting orchestration into four modules; shorten function names while retaining existing names as compatibility aliases and preserving fitting behavior and notebook outputs.
- Add initial package metadata, README, citation, and build-ignore skeletons.
- Move reusable code into a package namespace with explicit runtime dependencies.
- Replace automatic installation and project sourcing with `library(CrossCurveR)`.
- Run core tests against the package namespace and add clean-session checks.
- Make timestamp filtering independent of optional date-comparison methods.
- Normalize maturity labels before reconstruction and reject invalid or duplicate tenors.
- Add `dr_read_yields()` for validated CSV input without automatic header repair.
- Resolve `zero_mean_covreg()` initialization choices so the default uses adaptive initialization.
- Apply the covariance ridge once so regularized E-steps and likelihoods use the returned covariance.
- Handle scalar covariance square roots, reject material asymmetry or negative eigenvalues, and require numerically invertible covariances for whitening.
- Reject invalid iteration counts, checking intervals, and convergence tolerances before regularized covariance fitting.
- Move the production EM estimator into the core covariance-estimation module without changing its implementation or interface.
- Establish functionality, code-quality, and organization review gates with local before/after application comparisons.
- Share inversion, diagonal shrinkage, covariance regularization, and likelihood helpers within the estimator module while preserving estimator-specific algorithms.
- Share the covariance formula across reconstruction and simulation while preserving output shapes, tenor ordering, and simulation draws.
- Share the whitening, normalized fitting, and rescaling calculation through a private baseline helper while preserving existing H transformations and public fit outputs.
- Share lambda selection, single-curve NS fitting, loading construction, and residual extraction across ordinary and shocked fits through a private baseline helper.
- Reuse saved rolling baselines in the tontine shock applications, preserving post-baseline random draws for refits and rejecting incompatible baseline inputs.

<!-- TODO: record changes here as each implementation step is completed. -->
