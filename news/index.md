# Changelog

## fmrireg 0.2.0

### Datasets are now `fmri_frame` objects

- `fmrireg` has moved from the removed legacy dataset API of
  `fmridataset` (`matrix_dataset()`, `fmri_mem_dataset()`,
  `fmri_dataset()`, `latent_dataset()`, `get_data_matrix()`,
  `get_mask()`, `data_chunks()`, `read_fmri_config()`, …) to its
  canonical frame API. Every function that takes a `dataset`
  ([`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md),
  [`fmri_rlm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_rlm.md),
  [`fmri_model()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_model.md),
  [`create_fmri_model()`](https://bbuchsbaum.github.io/fmrireg/reference/create_fmri_model.md),
  [`estimate_betas()`](https://bbuchsbaum.github.io/fmrireg/reference/estimate_betas.md),
  [`estimate_hrf()`](https://bbuchsbaum.github.io/fmrireg/reference/estimate_hrf.md),
  [`fmri_latent_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_latent_lm.md),
  [`glm_ols()`](https://bbuchsbaum.github.io/fmrireg/reference/glm_ols.md),
  [`glm_lss()`](https://bbuchsbaum.github.io/fmrireg/reference/glm_lss.md),
  [`extract_nuisance_timeseries()`](https://bbuchsbaum.github.io/fmrireg/reference/extract_nuisance_timeseries.md),
  [`build_model()`](https://bbuchsbaum.github.io/fmrireg/reference/build_model.md),
  and the engine plugin API) now requires a
  [`fmridataset::fmri_frame`](https://bbuchsbaum.github.io/fmridataset/reference/fmri_frame.html)
  (or a lazy `fmri_view` of one). The old dataset classes are no longer
  accepted, and `fmridataset (>= 0.10.0.9000)` is required.

- Four thin constructors build a frame from the inputs `fmrireg` users
  typically hold; each stores the run structure as `run_id` / `TR`
  observation columns, an optional logical `censor` column, and the
  event table as a keyed
  [`fmridataset::event_table()`](https://bbuchsbaum.github.io/fmridataset/reference/event_table.html):

  - `matrix_frame(datamat, TR, run_length, event_table, censor)` for a
    time-by-feature matrix (an `index_space`; replaces
    `matrix_dataset()`).
  - `neurovec_frame(scans, mask, TR, ...)` for in-memory `NeuroVec` runs
    (a `volume_space` on the mask; replaces `fmri_mem_dataset()`).
  - `nifti_frame(scans, mask, TR, ...)` for NIfTI files, lazily read
    through
    [`fmridataset::nifti_array_source()`](https://bbuchsbaum.github.io/fmridataset/reference/nifti_array_source.html)
    (replaces `fmri_dataset()`).
  - `latent_frame(x, TR, run_length, ...)` for a
    [`fmristore::LatentNeuroVec`](https://rdrr.io/pkg/fmrilatent/man/LatentNeuroVec.html):
    component scores become the assay and the loadings a synthesis-only
    `basis_space` (replaces `latent_dataset()`). For fMRIPrep
    derivatives use
    [`fmridataset::read_bids_bold()`](https://bbuchsbaum.github.io/fmridataset/reference/read_bids_bold.html)
    directly.

- Reading a frame back uses the `fmridataset` API instead of list
  fields: `fmridataset::as_sampling_frame(x)` replaces
  `x$sampling_frame`, `fmridataset::temporal_schema(x)` exposes run
  lengths, TR, and censoring, `fmridataset::collect_assay(x)` replaces
  `get_data_matrix()`, the event table is
  `fmridataset::event_data(x$tables$events)`, and the spatial reference
  is recovered from `fmridataset::space(x)`. Censoring for
  `ar_options = list(censor = "auto")` is read from the frame’s `censor`
  column.

- [`estimate_betas()`](https://bbuchsbaum.github.io/fmrireg/reference/estimate_betas.md)
  has a single `fmri_frame` method. Frames on a `volume_space` return
  `NeuroVec` betas as the volumetric method did; frames on an
  `index_space` or `basis_space` return coefficient matrices as the
  matrix and latent methods did. The experimental `prewhiten` argument
  of the latent method, which never ran, has been dropped.

- [`dataset_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/dataset_spec.md)
  /
  [`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md)
  /
  [`realize_dataset()`](https://bbuchsbaum.github.io/fmrireg/reference/realize_dataset.md)
  name the new constructors: file-backed bindings produce
  `"nifti_frame"` specs and inline matrices produce `"matrix_frame"`
  specs; the constructor allowlist is `matrix_frame`, `nifti_frame`,
  `neurovec_frame`, and `latent_frame`.

- `simulate_fmri_matrix()$time_series` is now an `fmri_frame`.

- Engine capability `forbid_by_cluster_dataset_classes` is matched
  against the dataset’s class or its feature-space class
  (`"basis_space"` marks latent frames).

- `read_fmri_config()` is gone with the legacy API; nothing in `fmrireg`
  used it outside commented-out tests, so it was not ported.

- Fixed the `unused variable 'P'` compiler warning in
  [`ols_t_cpp()`](https://bbuchsbaum.github.io/fmrireg/reference/ols_t_cpp.md)
  (`src/ols_t.cpp`).

- Internally, chunkwise and runwise fitting iterate over lazy frame
  views (`x[rows, ]`, `x[, cols]` + `collect_assay()`) instead of
  `data_chunks()` / `exec_strategy()`, with the same chunk partition, so
  fitted numbers are unchanged.

### HRF Estimation

- [`estimate_hrf()`](https://bbuchsbaum.github.io/fmrireg/reference/estimate_hrf.md)
  is now a vectorized, condition-level smooth FIR estimator. It
  constructs an explicit event-aligned spline basis, removes baseline
  and fixed nuisance designs once, fits all voxels with one penalized
  multiresponse solve, and can choose a shared smoothing strength by
  scale-normalized GCV. The new `fmri_hrf_estimate` result preserves
  curve and voxel labels and provides standard errors, confidence
  intervals,
  [`tidy()`](https://bbuchsbaum.github.io/fmrireg/reference/tidy.md),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html), and
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) methods. This
  replaces a voxel-by-voxel GAM path that treated convolved design
  values as if they were post-stimulus time and failed before
  prediction.

### Statistical Corrections

- **AR degrees of freedom.**
  [`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
  no longer deflates the residual degrees of freedom when an AR
  structure is used. The previous adjustment multiplied `n - p` by
  `1 / (1 + 2 * sum(1 - k/n))`, which is the variance inflation one
  would obtain if *every* autocorrelation equalled 1. It never depended
  on the fitted AR coefficients — at `n = 200`, `p = 12`, AR(1) it
  returned 62.9 whether the true autocorrelation was 0.05 or 0.9 — and
  it was applied on top of prewhitening, double-counting a correction
  already made.

  **This changes reported AR p-values.** They become less conservative;
  degrees of freedom roughly triple for AR(1). The previous behaviour
  cost statistical power rather than inflating false positives, so
  results that were significant before remain significant, but p-values
  and any power-sensitive analyses will differ from earlier versions.

  `calculate_effective_df()` gains a working `method = "satterthwaite"`
  (previously a verbatim duplicate of `method = "simple"`) that computes
  `tr(RV)^2 / tr(RVRV)` from a supplied design and post-whitening
  covariance, reducing exactly to `n - p` when the errors are
  uncorrelated. Residual correlation surviving the filter is what
  legitimately costs degrees of freedom; AR order by itself does not.

### Bug Fixes

- The `fmridataset` requirement is now `(>= 0.11.0.9000)`. The pre-frame
  API was removed upstream without a version change, so both sides of
  the break reported `0.10.0.9000` and the previous constraint could not
  tell them apart: a stale install resolved cleanly and then failed at
  lazy loading with
  `object 'as.matrix_dataset' is not exported by 'namespace:fmridataset'`.
  fmridataset bumped to `0.11.0.9000` for this purpose
  (bbuchsbaum/fmridataset#92), so the mismatch is now caught at
  dependency resolution with a message that names the real problem.

- The test suite has been ported to the frame API. The frame migration
  and a parallel test-coverage branch were developed against the same
  base and merged independently, so the merged tree combined frame-only
  package code with ~50 test files still building legacy
  `matrix_dataset()`, `fmri_mem_dataset()`, and `latent_dataset()`
  fixtures. `R CMD check` failed on every platform with 106 test
  failures that neither branch saw on its own. Fixtures now use
  [`matrix_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md),
  [`neurovec_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/neurovec_frame.md),
  and
  [`latent_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/latent_frame.md);
  no assertion was relaxed to accommodate the port.

- `chunkwise_lm()` dispatches on `fmri_frame` from outside the package.
  [`chunkwise_lm.fmri_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/chunkwise_lm.fmri_frame.md)
  was never registered as an S3 method, so dispatch from outside the
  namespace failed with “no applicable method” even though internal
  calls resolved.

- `.rrr_extract_response_matrix()` validates its input again. A
  non-frame argument fell through to `collect_assay()` and failed inside
  `fmridataset` with a message naming neither the argument nor the
  caller.

- `extract_censor_from_dataset()` no longer errors on a non-frame
  dataset. It now reports “no censoring” instead, which matters because
  censoring is consulted on every fit under the default change below.

- A censor column carried by a dataset is now used by default.
  [`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
  previously consulted it only when `ar_options = list(censor = "auto")`
  was passed as well, so a `matrix_frame(censor = )` (and, before it,
  `fmri_dataset(censor = )`) was silently discarded: fits on a censored
  frame were bit-identical to fits on an uncensored one, with no
  warning. An unset `censor` now resolves against the dataset, which is
  what `"auto"` asked for explicitly; pass `censor = "none"` to ignore a
  censor column deliberately. **This changes results** for anyone who
  set censoring on a dataset and did not opt in, since those fits were
  not censored at all. An explicit `censor` vector in `ar_options` still
  overrides the dataset.

- [`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
  now warns when censoring cannot affect the fit. Censoring feeds AR
  estimation and whitening only and never drops volumes from the
  regression, so under the default iid noise model it does nothing at
  all. That was previously silent, and indistinguishable from censoring
  having been applied.

- Shared AR estimation now pools voxel residual autocovariances by
  default instead of fitting the cross-voxel mean residual series. The
  former targets a typical voxel covariance; the latter suppresses
  voxel-specific noise and targets a different, often more
  autocorrelated coherent component. The old behavior remains available
  explicitly as `noise_spec(..., shared_estimator = "mean_series")`.
  Fixed-order fits now use fmriAR’s pooled AR autocovariance path as
  well; the previous ARMA(p, 0) shortcut accidentally collapsed both
  choices to the mean series.

- `shared_estimator = "mean_series"` now retains the matching OLS
  design-bias correction. OLS projection is linear across response
  columns, so averaging residuals after projection is exactly the same
  as projecting the mean series.

- The OLS residual-bias correction now uses a separate covariance-tail
  lag budget rather than stopping at the fitted AR order. Projection
  mixes tail autocovariances into lags 0:p, so a p-lag correction left
  avoidable bias.

- Global AR estimation over separately fitted runs now supplies the
  exact block-diagonal residual-forming design. Row-binding the per-run
  designs incorrectly described a single coefficient vector shared
  across runs.

- Design-corrected AR coefficients are estimated once from initial OLS
  residuals and held fixed for the GLS solve. Later GLS residuals have a
  different residual-forming operator and no longer reuse the OLS
  `design=` correction. The low-rank engine now follows this contract
  too, and the low-rank and RRR whitening paths reset at run boundaries.

- GDS exports now release unreachable HDF5 handles before atomic
  finalization, allowing Windows to move and immediately reopen the
  completed file.

- `write_results(strategy = "by_stat")` no longer scatters contrast maps
  onto the wrong voxels. `.compute_statistical_volumes()` called
  [`as.logical()`](https://rdrr.io/r/base/logical.html) on the mask,
  which drops `dim` for plain arrays (the type
  `.fmri_dataset_mask_space()` always supplies), so
  `which(..., arr.ind = TRUE)` returned linear indices instead of 3D
  coordinates.

- Oversized residual-bootstrap blocks now remain one contiguous temporal
  block instead of silently becoming interleaved odd/even samples.

- Mixed-model solver failures now warn and return `NA` coefficients
  rather than plausible zeros; malformed response dimensions fail before
  solver dispatch.

- `robust_psi` now has an effect. It was previously unreachable:
  `robust` defaulted to `FALSE` rather than `NULL`, so
  `robust_options$type` was always claimed before `robust_psi` was
  consulted. `robust`’s default is now `NULL`, which also makes
  “unspecified” distinguishable from “explicitly off”.

- An explicit `robust = FALSE` is no longer silently overridden by
  `robust_options = list(type = "huber")`; the combination is now an
  error.

- Supplying `cfg` alongside other configuration arguments is now an
  error instead of silently discarding them. Previously
  `cor_struct = "ar2"` with a `cfg` naming `"iid"` ran IID without
  warning.

- AR shorthands (`cor_struct`, `cor_iter`, `cor_global`,
  `ar1_exact_first`, `ar_p`, `ar_voxelwise`) that disagree with the
  corresponding `ar_options` entry now error rather than being silently
  dropped. Shorthands that agree are still accepted.

- `engine` now warns about the execution arguments it ignores
  (`strategy`, `nchunks`, `use_fast_path`, `progress`,
  `parallel_voxels`, `parallel_chunks`), which were previously accepted
  and silently discarded.

- `compute_sandwich_variance()` computed its meat matrix as
  `X' diag(e^4) X` instead of `X' diag(e^2) X`, giving standard errors
  about 3.5x too large. Both sandwich helpers are now checked against
  [`sandwich::vcovHC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html).

### Performance

- Voxelwise AR fitting no longer recomputes the run-level design
  projection for every voxel. `.fast_preproject()` performs an `n x n`
  solve and was being called once per voxel on a design that does not
  vary by voxel; hoisting it gives roughly a 1.7x speedup on the
  voxelwise AR path (4.95s to 2.98s for 2000 voxels at 300 timepoints).

### Documentation

- `ar_voxelwise` was documented as overriding `ar_options$voxelwise`. It
  never did — the options list took precedence — and disagreement is now
  an error.
- `use_fast_path` was documented as defaulting to `FALSE` in the runwise
  and chunkwise fitters; it defaults to `TRUE`.
- Corrected claims that sandwich variance estimation is applied
  automatically. It is not used by any
  [`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
  fitting path; the robust path reports a model-based
  weighted-least-squares variance. The accompanying degrees-of- freedom
  formula was also wrong and has been rewritten.

### Changes

- `fmriAR` is now declared with a minimum version (`>= 0.3.3`) and
  listed in `Remotes:`; it previously had neither.

### New Features

- Reexported
  [`feature()`](https://bbuchsbaum.github.io/fmridesign/reference/feature.html)
  (from fmridesign) and
  [`feature_regressor()`](https://bbuchsbaum.github.io/fmrihrf/reference/feature_regressor.html)
  (from fmrihrf) so mixed event formulas such as
  `onset ~ hrf(condition) + feature(rms, dt = 0.1)` work after
  [`library(fmrireg)`](https://bbuchsbaum.github.io/fmrireg/). Feature
  terms are sampled series, not trials; helpers that assume `hrfspec`
  trial structure skip or unwrap them instead of erroring
  ([`fitted_hrf()`](https://bbuchsbaum.github.io/fmrireg/reference/fitted_hrf.md),
  [`shortnames()`](https://bbuchsbaum.github.io/fmrireg/reference/shortnames.md),
  [`design_plot()`](https://bbuchsbaum.github.io/fmrireg/reference/design_plot.md),
  [`preflight()`](https://bbuchsbaum.github.io/fmrireg/reference/preflight.md),
  pair-contrast resolution).

- [`write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/write_results.md)
  now exports statistical maps as NIfTI volumes via `format = "nifti"`
  (or `format = c("h5", "nifti")`), reusing the same BIDS entity and
  filename machinery as the HDF5 backend. This removes the need to
  hand-roll loops over
  [`coef_image()`](https://bbuchsbaum.github.io/fmrireg/reference/coef_image.md) +
  [`neuroim2::write_vol()`](https://bbuchsbaum.github.io/neuroim2/reference/write_vol-methods.html).

- Added
  [`coef_images()`](https://bbuchsbaum.github.io/fmrireg/reference/coef_images.md),
  a plural companion to
  [`coef_image()`](https://bbuchsbaum.github.io/fmrireg/reference/coef_image.md)
  that returns a named list of `NeuroVol`s for every coefficient of a
  given `type` (`"estimates"`, `"contrasts"`, or `"F"`) and `statistic`.

### Bug Fixes

- [`write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/write_results.md)
  now labels F-contrast maps as `fstat`/`fpval` instead of reusing the
  `tstat`/`pval` labels. F statistics follow an F distribution (not
  Student’s t), so mixing them under a `tstat` label produced mislabeled
  maps and, under `strategy = "by_stat"`, could collide t- and
  F-contrast outputs in a single file. t- and F-contrasts are now
  written to separate, correctly labeled files.

## fmrireg 0.1.2

### Changes

- Moved PDF report generation (`report.fmri_lm`) to the standalone
  `fmrireport` package. Use `fmrireport::report()` for report rendering.

### Bug Fixes

- Fixed low-rank/sketch residual degrees of freedom to use original
  timepoints (`T - p`) instead of sketch rows.
- Replaced placeholder
  [`fit_contrasts.fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fit_contrasts.fmri_lm.md)
  outputs with computed standard errors, test statistics, and p-values.
- Consolidated duplicate meta-analysis pathways to shared
  implementations and removed inconsistent behavior across
  runwise/chunkwise pooling.
- Fixed F-contrast meta pooling to combine evidence via p-values rather
  than averaging F statistics.
- Fixed thread-safety issue in C++ meta kernels by deferring
  Paule-Mandel non-convergence warnings until after OpenMP regions.
- Added bounds validation and safer inverse fallback in sketch kernels
  to prevent unsafe indexing and hard failures on ill-conditioned data.
- Improved AR pipeline consistency for whitened
  covariance/fitted/residual handling and unified AR effective-df
  calculations.
- Added regression coverage for
  [`fit_contrasts.fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fit_contrasts.fmri_lm.md)
  to prevent placeholder-stat regressions.

------------------------------------------------------------------------

## fmrireg 0.1.1

### Bug Fixes

- Fixed
  [`glm_lss()`](https://bbuchsbaum.github.io/fmrireg/reference/glm_lss.md)
  tests that incorrectly expected specific Cholesky decomposition
  errors.
- Fixed convolution test with non-strictly-increasing onsets within
  blocks.
- Fixed testthat API usage (`expect_lt`/`expect_gt` no longer use
  deprecated `info` argument).
- Fixed `latent_dataset` API usage: now uses `get_latent_scores()`
  instead of deprecated `get_data()`.
- Changed
  [`glm_lss()`](https://bbuchsbaum.github.io/fmrireg/reference/glm_lss.md)
  `use_cpp` parameter default from `TRUE` to `FALSE` (C++ implementation
  retired; fmrilss package now used).
- Added `fmrireg.suppress_deprecation` option check to all deprecated
  functions for cleaner test output.
- Added `tests/testthat/setup.R` to suppress expected deprecation
  warnings during testing.
- Suppressed expected kmeans convergence warnings in landmark SRHT
  tests.

### Internal

- Test warnings reduced from 129 to 31 (remaining warnings are from
  external packages).

------------------------------------------------------------------------

## fmrireg 0.1.0

### Breaking Changes

- **Design Matrix Column Naming:** The naming scheme for columns in
  design matrices generated by
  [`event_model()`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html)
  has been completely revised for consistency and clarity.
  - All column names now strictly follow the format: `term_tag` + `_` +
    `condition_tag` + \[`_b##` basis suffix\].
  - `term_tag`: Automatically generated from variable names (e.g.,
    `var1_var2`) or user-provided `id=` in
    [`hrf()`](https://bbuchsbaum.github.io/fmridesign/reference/hrf.html),
    sanitized (dots become underscores), and made unique with `#` suffix
    if needed (e.g., `cond`, `cond#1`).
  - `condition_tag`: Represents factor levels (e.g., `Factor.Level`),
    continuous basis columns (e.g., `poly_RT_01`, `z_RT`), or
    interactions joined by `_` (e.g., `Factor.Level_poly_RT_01`).
  - `_b##`: Optional suffix added *only* when the HRF has multiple basis
    functions (e.g., `_b01`, `_b02`).
  - The previous `style` argument (`"compact"`, `"qualified"`, `"uid"`)
    in
    [`design_matrix()`](https://bbuchsbaum.github.io/fmridesign/reference/design_matrix.html)
    is removed. Only the single canonical format is produced.
  - Scripts or analyses that relied on matching previous column name
    formats (e.g., using `Var[Level]`, `Var:Level`, `:basis[]`) **will
    need to be updated** to use the new `term_tag_Condition.Tag_b##`
    format.

### Major Changes

- **Regressor System Refactoring:**
  - Introduced a new internal S3 class `Reg` for representing
    regressors.
  - The main
    [`regressor()`](https://bbuchsbaum.github.io/fmrihrf/reference/regressor.html)
    function now uses `Reg` internally but maintains backward
    compatibility (returns class `c("regressor", "Reg", "list")`).
  - Deprecated
    [`single_trial_regressor()`](https://bbuchsbaum.github.io/fmrihrf/reference/single_trial_regressor.html)
    and `null_regressor()` in favour of using
    [`regressor()`](https://bbuchsbaum.github.io/fmrihrf/reference/regressor.html)
    directly.
  - Unified regressor evaluation under the `evaluate.Reg` S3 method,
    supporting different calculation methods (“fft”, “conv”, “loop”,
    “Rconv”).
  - Evaluation methods now consistently use the refactored
    `evaluate.HRF` for HRF sampling.
  - Refactored C++ evaluation code into a single wrapper
    (`evaluate_regressor_cpp`).
  - Removed redundant internal helper functions (`fastevalreg`,
    `fastevalreg2`, `conform_len`, `dots`).
  - Implemented `autoplot.Reg` (ggplot2) and `print.Reg` (cli) methods,
    deprecating older `plot.regressor` and `print.regressor`.
  - Improved input validation and recycling for
    [`regressor()`](https://bbuchsbaum.github.io/fmrihrf/reference/regressor.html)
    arguments using `vctrs`.
  - Added optional sparse matrix output to `evaluate.Reg`.
  - Added memoization for HRF sampling within evaluation.

------------------------------------------------------------------------

## fmrireg 0.0.1

- Initial CRAN release.
