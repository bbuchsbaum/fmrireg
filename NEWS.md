# fmrireg 0.2.0

## HRF Estimation

* `estimate_hrf()` is now a vectorized, condition-level smooth FIR estimator.
  It constructs an explicit event-aligned spline basis, removes baseline and
  fixed nuisance designs once, fits all voxels with one penalized
  multiresponse solve, and can choose a shared smoothing strength by
  scale-normalized GCV. The new `fmri_hrf_estimate` result preserves curve and
  voxel labels and provides standard errors, confidence intervals, `tidy()`,
  `predict()`, `coef()`, and `as.matrix()` methods. This replaces a
  voxel-by-voxel GAM path that treated convolved design values as if they were
  post-stimulus time and failed before prediction.

## Statistical Corrections

* **AR degrees of freedom.** `fmri_lm()` no longer deflates the residual
  degrees of freedom when an AR structure is used. The previous adjustment
  multiplied `n - p` by `1 / (1 + 2 * sum(1 - k/n))`, which is the variance
  inflation one would obtain if *every* autocorrelation equalled 1. It never
  depended on the fitted AR coefficients — at `n = 200`, `p = 12`, AR(1) it
  returned 62.9 whether the true autocorrelation was 0.05 or 0.9 — and it was
  applied on top of prewhitening, double-counting a correction already made.

  **This changes reported AR p-values.** They become less conservative;
  degrees of freedom roughly triple for AR(1). The previous behaviour cost
  statistical power rather than inflating false positives, so results that
  were significant before remain significant, but p-values and any
  power-sensitive analyses will differ from earlier versions.

  `calculate_effective_df()` gains a working `method = "satterthwaite"`
  (previously a verbatim duplicate of `method = "simple"`) that computes
  `tr(RV)^2 / tr(RVRV)` from a supplied design and post-whitening covariance,
  reducing exactly to `n - p` when the errors are uncorrelated. Residual
  correlation surviving the filter is what legitimately costs degrees of
  freedom; AR order by itself does not.

## Bug Fixes

* Oversized residual-bootstrap blocks now remain one contiguous temporal block
  instead of silently becoming interleaved odd/even samples.
* Mixed-model solver failures now warn and return `NA` coefficients rather than
  plausible zeros; malformed response dimensions fail before solver dispatch.

* `robust_psi` now has an effect. It was previously unreachable: `robust`
  defaulted to `FALSE` rather than `NULL`, so `robust_options$type` was always
  claimed before `robust_psi` was consulted. `robust`'s default is now `NULL`,
  which also makes "unspecified" distinguishable from "explicitly off".
* An explicit `robust = FALSE` is no longer silently overridden by
  `robust_options = list(type = "huber")`; the combination is now an error.
* Supplying `cfg` alongside other configuration arguments is now an error
  instead of silently discarding them. Previously `cor_struct = "ar2"` with a
  `cfg` naming `"iid"` ran IID without warning.
* AR shorthands (`cor_struct`, `cor_iter`, `cor_global`, `ar1_exact_first`,
  `ar_p`, `ar_voxelwise`) that disagree with the corresponding `ar_options`
  entry now error rather than being silently dropped. Shorthands that agree
  are still accepted.
* `engine` now warns about the execution arguments it ignores (`strategy`,
  `nchunks`, `use_fast_path`, `progress`, `parallel_voxels`,
  `parallel_chunks`), which were previously accepted and silently discarded.
* `compute_sandwich_variance()` computed its meat matrix as
  `X' diag(e^4) X` instead of `X' diag(e^2) X`, giving standard errors about
  3.5x too large. Both sandwich helpers are now checked against
  `sandwich::vcovHC()`.

## Performance

* Voxelwise AR fitting no longer recomputes the run-level design projection
  for every voxel. `.fast_preproject()` performs an `n x n` solve and was
  being called once per voxel on a design that does not vary by voxel;
  hoisting it gives roughly a 1.7x speedup on the voxelwise AR path
  (4.95s to 2.98s for 2000 voxels at 300 timepoints).

## Documentation

* `ar_voxelwise` was documented as overriding `ar_options$voxelwise`. It never
  did — the options list took precedence — and disagreement is now an error.
* `use_fast_path` was documented as defaulting to `FALSE` in the runwise and
  chunkwise fitters; it defaults to `TRUE`.
* Corrected claims that sandwich variance estimation is applied automatically.
  It is not used by any `fmri_lm()` fitting path; the robust path reports a
  model-based weighted-least-squares variance. The accompanying degrees-of-
  freedom formula was also wrong and has been rewritten.

## Changes

* `fmriAR` is now declared with a minimum version (`>= 0.3.3`) and listed in
  `Remotes:`; it previously had neither.

## New Features

* Reexported `feature()` (from fmridesign) and `feature_regressor()` (from
  fmrihrf) so mixed event formulas such as
  `onset ~ hrf(condition) + feature(rms, dt = 0.1)` work after
  `library(fmrireg)`. Feature terms are sampled series, not trials; helpers
  that assume `hrfspec` trial structure skip or unwrap them instead of
  erroring (`fitted_hrf()`, `shortnames()`, `design_plot()`, `preflight()`,
  pair-contrast resolution).

* `write_results()` now exports statistical maps as NIfTI volumes via
  `format = "nifti"` (or `format = c("h5", "nifti")`), reusing the same
  BIDS entity and filename machinery as the HDF5 backend. This removes the
  need to hand-roll loops over `coef_image()` + `neuroim2::write_vol()`.
* Added `coef_images()`, a plural companion to `coef_image()` that returns a
  named list of `NeuroVol`s for every coefficient of a given `type`
  (`"estimates"`, `"contrasts"`, or `"F"`) and `statistic`.

## Bug Fixes

* `write_results()` now labels F-contrast maps as `fstat`/`fpval` instead of
  reusing the `tstat`/`pval` labels. F statistics follow an F distribution
  (not Student's t), so mixing them under a `tstat` label produced
  mislabeled maps and, under `strategy = "by_stat"`, could collide t- and
  F-contrast outputs in a single file. t- and F-contrasts are now written to
  separate, correctly labeled files.

# fmrireg 0.1.2

## Changes

* Moved PDF report generation (`report.fmri_lm`) to the standalone
  `fmrireport` package. Use `fmrireport::report()` for report rendering.

## Bug Fixes

* Fixed low-rank/sketch residual degrees of freedom to use original timepoints (`T - p`) instead of sketch rows.
* Replaced placeholder `fit_contrasts.fmri_lm()` outputs with computed standard errors, test statistics, and p-values.
* Consolidated duplicate meta-analysis pathways to shared implementations and removed inconsistent behavior across runwise/chunkwise pooling.
* Fixed F-contrast meta pooling to combine evidence via p-values rather than averaging F statistics.
* Fixed thread-safety issue in C++ meta kernels by deferring Paule-Mandel non-convergence warnings until after OpenMP regions.
* Added bounds validation and safer inverse fallback in sketch kernels to prevent unsafe indexing and hard failures on ill-conditioned data.
* Improved AR pipeline consistency for whitened covariance/fitted/residual handling and unified AR effective-df calculations.
* Added regression coverage for `fit_contrasts.fmri_lm()` to prevent placeholder-stat regressions.

---

# fmrireg 0.1.1

## Bug Fixes

*   Fixed `glm_lss()` tests that incorrectly expected specific Cholesky decomposition errors.
*   Fixed convolution test with non-strictly-increasing onsets within blocks.
*   Fixed testthat API usage (`expect_lt`/`expect_gt` no longer use deprecated `info` argument).
*   Fixed `latent_dataset` API usage: now uses `get_latent_scores()` instead of deprecated `get_data()`.
*   Changed `glm_lss()` `use_cpp` parameter default from `TRUE` to `FALSE` (C++ implementation retired; fmrilss package now used).
*   Added `fmrireg.suppress_deprecation` option check to all deprecated functions for cleaner test output.
*   Added `tests/testthat/setup.R` to suppress expected deprecation warnings during testing.
*   Suppressed expected kmeans convergence warnings in landmark SRHT tests.

## Internal

*   Test warnings reduced from 129 to 31 (remaining warnings are from external packages).

---

# fmrireg 0.1.0

## Breaking Changes

*   **Design Matrix Column Naming:** The naming scheme for columns in design matrices generated by `event_model()` has been completely revised for consistency and clarity.
    *   All column names now strictly follow the format: `term_tag` + `_` + `condition_tag` + [`_b##` basis suffix].
    *   `term_tag`: Automatically generated from variable names (e.g., `var1_var2`) or user-provided `id=` in `hrf()`, sanitized (dots become underscores), and made unique with `#` suffix if needed (e.g., `cond`, `cond#1`).
    *   `condition_tag`: Represents factor levels (e.g., `Factor.Level`), continuous basis columns (e.g., `poly_RT_01`, `z_RT`), or interactions joined by `_` (e.g., `Factor.Level_poly_RT_01`).
    *   `_b##`: Optional suffix added *only* when the HRF has multiple basis functions (e.g., `_b01`, `_b02`).
    *   The previous `style` argument (`"compact"`, `"qualified"`, `"uid"`) in `design_matrix()` is removed. Only the single canonical format is produced.
    *   Scripts or analyses that relied on matching previous column name formats (e.g., using `Var[Level]`, `Var:Level`, `:basis[]`) **will need to be updated** to use the new `term_tag_Condition.Tag_b##` format.

## Major Changes

*   **Regressor System Refactoring:**
    *   Introduced a new internal S3 class `Reg` for representing regressors.
    *   The main `regressor()` function now uses `Reg` internally but maintains backward compatibility (returns class `c("regressor", "Reg", "list")`).
    *   Deprecated `single_trial_regressor()` and `null_regressor()` in favour of using `regressor()` directly.
    *   Unified regressor evaluation under the `evaluate.Reg` S3 method, supporting different calculation methods ("fft", "conv", "loop", "Rconv").
    *   Evaluation methods now consistently use the refactored `evaluate.HRF` for HRF sampling.
    *   Refactored C++ evaluation code into a single wrapper (`evaluate_regressor_cpp`).
    *   Removed redundant internal helper functions (`fastevalreg`, `fastevalreg2`, `conform_len`, `dots`).
    *   Implemented `autoplot.Reg` (ggplot2) and `print.Reg` (cli) methods, deprecating older `plot.regressor` and `print.regressor`.
    *   Improved input validation and recycling for `regressor()` arguments using `vctrs`.
    *   Added optional sparse matrix output to `evaluate.Reg`.
    *   Added memoization for HRF sampling within evaluation.

---

# fmrireg 0.0.1

* Initial CRAN release. 
