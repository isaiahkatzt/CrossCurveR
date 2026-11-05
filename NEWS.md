# CrossCurveR (development version)

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

<!-- TODO: record changes here as each implementation step is completed. -->
