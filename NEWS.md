# fmrireg 0.3.0

## Breaking: design generics now come from fmridesign (fmridesign >= 0.6.1.9000)

fmrireg no longer defines its own copies of `longnames()`, `shortnames()`,
`construct()` and `correlation_map()`, and no longer registers methods for
classes that fmridesign owns (`event_model`, `event_term`, `convolved_term`,
`feature_term`, `baseline_model`, ...). It re-exports fmridesign's generics,
so `fmrireg::longnames` and `fmridesign::longnames` are the same function and
a method registered by either package is reached from both. `.onLoad` no
longer reaches into fmridesign's namespace. Coefficient and design-matrix
column names do not change. What does change:

* **`blockids(<event_model>)` now returns one run id per event, not one per
  scan.** fmrireg used to override fmridesign's method with per-scan ids
  whenever it was loaded, so the result depended on whether fmrireg was
  attached. Code that indexes scans with it must change:

  ```r
  blockids(em)                  # per event (fmridesign's method)
  blockids(em$sampling_frame)   # per scan: what blockids(em) used to return
  blocklens(em)                 # scans per run (unchanged)
  ```

  fmridesign prints a one-time message the first time `blockids()` is called
  on an `event_model` in a session.

* **`longnames()` uses fmridesign's `.` format instead of `#`.**

  | design | before | now |
  |---|---|---|
  | `hrf(cond)` | `cond#A` | `cond.A` |
  | `hrf(rt)` | `rt#rt` | `rt` |
  | `hrf(cond, rt)` | `cond#A` | `cond.A_rt` |
  | `hrf(cond, attn)` | `cond#A:attn#x` | `cond.A_attn.x` |
  | `trialwise()` | `.trial_factor...#01` | `.trial_factor....01` |

  A long name is the design-matrix column name without its `<term>_` prefix
  (and without the `_bNN` basis suffix unless `expand_basis = TRUE`): column
  `cond_cond.A` has long name `cond.A`. See `?fmridesign::longnames`.
  `shortnames()` is unchanged for factor terms, but parametric modulators now
  keep the modulator (`hrf(cond, rt)` gives `A:rt`, formerly `A`), and the
  result for an `event_model` is an unnamed vector (formerly named).

* **Empty interaction cells** are still left out of `longnames()` and
  `shortnames()` (the default, `drop.empty = TRUE`), so names line up with
  columns. `drop.empty = FALSE` lists the full factor grid. `conditions()`
  always returns the full grid, including empty cells.

* **`correlation_map(<baseline_model>)`** is fmridesign's method. It keeps
  fmrireg's former behaviour: `within_run = TRUE` by default (run intercepts
  dropped, columns centred within runs, run-specific columns correlated on
  their own run), with identical cells and values. `label_values` is accepted
  as an alias for `annotate`, and cells are labelled by default at 12 or
  fewer columns, as before. The drawing (VIF diagonal, outlined high
  correlations) follows fmridesign's style. Arguments that are not
  `geom_tile()` arguments now raise an error instead of being ignored.

* `design_matrix(<convolved_term>, blockid = )` keeps a one-column result as a
  matrix.

## Other fixes in this change

* `standard_error(<fmri_latent_lm>)` (including `recon = TRUE`) and the
  internal `pull_stat_revised()` name their columns with the design-matrix
  column names, matching `coef()` and `stats(recon = TRUE)`. For simple
  designs the names change from the condition name to the column name (for
  example `a.1` becomes `a_a.1`). Multi-basis HRFs and designs with an empty
  interaction cell, where the number of conditions differs from the number of
  columns, used to fail with an error and now work.
* `glm_ols()` and `glm_lss()`, which accept an HRF basis by name
  (`"HRF_SPMG1"` etc.), now resolve it with `getExportedValue()` and accept
  any `HRF_*` object that fmrihrf exports. The names `"HRF_AFNI"`,
  `"HRF_GAM"`, `"HRF_IL"` and `"HRF_DD"` were listed as valid but never
  existed in fmrihrf; they now give "Unknown HRF basis name" instead of an
  "object not found" error.
* The fmrigds reducer check calls `fmrigds::get_reducer()` directly instead
  of looking it up in fmrigds's namespace.

## Breaking: corrected SPMG kernel (fmrihrf >= 0.4.0)

* fmrireg now builds against fmrihrf 0.4.0, which corrects the SPM canonical
  HRF used by `HRF_SPMG1`, `HRF_SPMG2`, and `HRF_SPMG3`. The kernel now has the
  SPM double-gamma shape (about 9% undershoot rather than 0.6%), and its raw
  peak is about ten times lower (0.175 rather than 1.75). For the same data,
  betas and standard errors for the canonical and temporal-derivative columns
  are therefore about ten times larger; t and F statistics change only through
  the corrected shape. The third `HRF_SPMG3` column is now a dispersion
  derivative, not a second time derivative, so its coefficients are not
  comparable across versions at all. Because `hrf()` defaults to an
  unnormalized SPMG1 basis, this affects default models. Refit existing
  analyses rather than comparing raw coefficients across versions.
* The bundled `fmri_benchmark_datasets` were regenerated with the corrected
  kernel (same seeds), and the data-raw generator was updated for the
  `fmri_frame` dataset API.
* NIfTI maps written by `write_results()` are now explicitly FLOAT32 (NIfTI
  datatype 16), as documented. This was already the neuroim2 default for
  in-memory volumes; it now also holds for volumes whose source file was
  stored as DOUBLE.

## Breaking change: `coef()` on `fmri_lm` fits is now voxels x terms

* `coef.fmri_lm()` now returns one orientation for every argument
  combination: **one row per voxel, one column per coefficient**, with term
  names on the column margin (#227). Previously the default call
  (`type = "betas"`, `include_baseline = FALSE`) returned the transpose,
  terms x voxels with term names as row names, while
  `include_baseline = TRUE` and `type = "contrasts"` were already
  voxels x terms. `coef()` now also agrees with `stats()`,
  `standard_error()` and `p_values()`, so `coef(fit) / standard_error(fit)`
  lines up element-wise.

  **Update code that indexed the old default by row:** `coef(fit)[term, ]`
  becomes `coef(fit)[, term]`, `rownames(coef(fit))` becomes
  `colnames(coef(fit))` (or `coef_names(fit)`), and a `t(coef(fit))` added
  to reach voxels x terms must be dropped. When the number of voxels equals
  the number of event coefficients the old and new results have the same
  `dim()`, so such code does not error: check it by hand. The orientation
  flip has no deprecation period, since an orientation cannot warn only
  the callers who depend on it. `coef()` now has its own help page,
  `?coef.fmri_lm`.

* The old orientation caused two silent transpositions, both now fixed:
  - `fit_contrasts()` on an `fmri_lm` fit guessed the orientation of `coef()`
    from its dimensions. When the voxel count equalled the number of event
    coefficients, it computed the contrast from the transposed betas.
  - The "known vs recovered" table in the Package Overview vignette (two
    voxels, two conditions) showed each condition's estimates across voxels
    next to each voxel's true values.

* `reduce_betas()` and latent-space reconstruction in
  `coef(<fmri_latent_lm>, recon = TRUE)` also compensated for the old
  orientation. They now use the single one, and their output is unchanged.

## Runwise fits: correct baseline columns and near-exact fits

* Multi-run runwise fits pooled betas by run-local column position. Each
  run's design has the event columns plus only that run's baseline columns,
  so the stored betas were too narrow (6 of 10 columns for two runs with
  default drift). They mixed different runs' drift and intercept terms, and
  `coef(fit, include_baseline = TRUE)` gave them the wrong global names, for
  example run 1's intercept labelled `base_bs1_block_2`. Betas are now pooled
  by global design column. A run-specific baseline column keeps its own
  run's estimate. A column shared by all runs (event regressors, a global
  intercept) is pooled by inverse-variance weighting. The stored betas now
  have one column per design-matrix column. Event-column estimates are
  unchanged.

* Runwise pooling returned `NaN` for every coefficient when a run's
  standard error was exactly zero. With inverse-variance weighting, a
  near-exact fit made 0/0. Pooling now takes the zero-variance limit.

* The memory-lean residual-sum-of-squares computation (`y'y - b'X'y`)
  cancelled to zero for near-exact fits, which reported `se = 0`. Voxels
  whose RSS falls below the precision of `y'y` are now recomputed from
  explicit residuals.

* The internal demo fit used in examples and tests now has 3 voxels instead
  of 2, so its voxels x terms results are not square.

## `stats(fit, "estimates")` renamed to `stats(fit, "betas")`

* `stats(fit, type = "estimates")` returned t-statistics, and the name
  invited reading them as estimates (#217). The parameter family is now
  called `"betas"`, matching `coef(type = "betas")`, and that is the default
  for `stats()`, `standard_error()` and `p_values()` on `fmri_lm` and
  `fmri_latent_lm` fits. `stats(fit, type = "estimates")` still works but
  raises a deprecation warning (class
  `fmrireg_deprecated_estimates_type`) that points to `coef()`; set
  `options(fmrireg.suppress_deprecation = TRUE)` to silence it.
  `standard_error()` and `p_values()` accept `"estimates"` as a silent
  synonym, because their function names already say what they return. The
  `type` documentation of `stats()` now says it returns t-statistics and
  links to `coef()`.

* The coefficient view of `autoplot()` / `plot()` for `fmri_lm` is named
  `type = "betas"` to match. `type = "estimates"` stays a silent synonym,
  because the plot is labelled and cannot be mistaken for the estimates
  themselves. The plots now read coefficient values from `coef()` instead of
  rebuilding them as t x SE, which gave 0 where the SE was numerically zero.

* Duplicate definitions of `coef()`, `stats()` and `standard_error()` for
  `fmri_lm` in `R/fmrilm.R` and `R/fmri_lm_methods.R` have been merged
  into one. Which copy ran used to depend on file collation order.

## Plotting

* New `autoplot()` / `plot()` methods for `fmri_lm` fits, with five views
  selected by `type`:
  - `"betas"` / `"contrasts"`: t-statistics across voxels as a sina plot
    with p < .001 reference lines and the share of voxels beyond each tail;
    with `voxel =`, a coefficient plot with confidence intervals, one panel
    per event term.
  - `"hrf"`: estimated responses with pointwise confidence bands from the
    coefficient covariance (one voxel) or a t-based interval of the mean
    (several voxels). Peak latency is labelled only when the basis has more
    than one function.
  - `"timecourse"`: observed and fitted signal for one voxel, onset lanes per
    condition, and residuals on the same scale, with the partial R-squared of
    the event regressors.
  - `"residuals"`: residual autocorrelation at lags 1-10, by default the
    median and 10th-90th percentile over up to 200 voxels. The white-noise
    range is centred on the (negative) autocorrelation that least-squares
    residuals have under the design and uses the exact null covariance, and
    the subtitle counts voxels failing a portmanteau test. For AR fits the
    prewhitened (GLS) residuals are shown alongside, which checks that the
    noise model removed the serial correlation.
* New `theme_fmrireg()`, used by every fmrireg plot. Set
  `options(fmrireg.plot_theme = <theme>)` to substitute your own.
* All plots share one colour-blind-safe condition palette and readable
  regressor labels (`condition_condition.faces` is shown as `faces`),
  built from term metadata.
* `design_map()` draws the run boundary at the correct scan (it was placed
  using per-event block ids), keeps columns in model order grouped by term,
  scales each column to its maximum absolute value by default
  (`scale_columns = TRUE`), draws run intercepts in grey, and reports the
  rank and condition number of the design. Extra arguments now go to
  `geom_raster()`, so `geom_tile()` arguments such as `colour =` are
  ignored with a warning.
* `correlation_map()` now defaults to `half_matrix = TRUE` and a new
  `within_run = TRUE`, which centres columns within runs and drops run
  intercepts before correlating, blanks pairs of columns from different runs,
  and correlates pairs involving a run-specific column over that run's scans.
  The diagonal is left blank, constant columns are shown in grey, event
  regressors' VIFs are given in the caption, and `absolute_limits = FALSE`
  now spans plus or minus the largest absolute correlation.
  `within_run = FALSE, half_matrix = FALSE` is closest to the previous
  behaviour.
* `plot.fmri_model()` is now implemented in fmrireg (it no longer delegates to
  fmridesign or needs cowplot): small multiples of every regressor on one time
  axis, with nuisance regressors included, run-specific columns drawn only in
  their run, and no lines joined across runs.
* `autoplot()` for regressors marks event onsets, draws single-event responses
  under the summed curve when they overlap, extends the time axis to cover
  event durations, and names SPMG2/SPMG3 basis functions.
* `cowplot` is no longer a suggested package.

## Datasets are now `fmri_frame` objects

* `fmrireg` has moved from the removed legacy dataset API of `fmridataset`
  (`matrix_dataset()`, `fmri_mem_dataset()`, `fmri_dataset()`,
  `latent_dataset()`, `get_data_matrix()`, `get_mask()`, `data_chunks()`,
  `read_fmri_config()`, ...) to its canonical frame API. Every function that
  takes a `dataset` (`fmri_lm()`, `fmri_rlm()`, `fmri_model()`,
  `create_fmri_model()`, `estimate_betas()`, `estimate_hrf()`,
  `fmri_latent_lm()`, `glm_ols()`, `glm_lss()`, `extract_nuisance_timeseries()`,
  `build_model()`, and the engine plugin API) now requires a
  `fmridataset::fmri_frame` (or a lazy `fmri_view` of one). The old dataset
  classes are no longer accepted, and `fmridataset (>= 0.10.0.9000)` is
  required.

* Four thin constructors build a frame from the inputs `fmrireg` users
  typically hold; each stores the run structure as `run_id` / `TR`
  observation columns, an optional logical `censor` column, and the event
  table as a keyed `fmridataset::event_table()`:
  - `matrix_frame(datamat, TR, run_length, event_table, censor)` for a
    time-by-feature matrix (an `index_space`; replaces `matrix_dataset()`).
  - `neurovec_frame(scans, mask, TR, ...)` for in-memory `NeuroVec` runs
    (a `volume_space` on the mask; replaces `fmri_mem_dataset()`).
  - `nifti_frame(scans, mask, TR, ...)` for NIfTI files, lazily read through
    `fmridataset::nifti_array_source()` (replaces `fmri_dataset()`).
  - `latent_frame(x, TR, run_length, ...)` for a `fmristore::LatentNeuroVec`:
    component scores become the assay and the loadings a synthesis-only
    `basis_space` (replaces `latent_dataset()`).
  For fMRIPrep derivatives use `fmridataset::read_bids_bold()` directly.

* Reading a frame back uses the `fmridataset` API instead of list fields:
  `fmridataset::as_sampling_frame(x)` replaces `x$sampling_frame`,
  `fmridataset::temporal_schema(x)` exposes run lengths, TR, and censoring,
  `fmridataset::collect_assay(x)` replaces `get_data_matrix()`, the event
  table is `fmridataset::event_data(x$tables$events)`, and the spatial
  reference is recovered from `fmridataset::space(x)`. Censoring for
  `ar_options = list(censor = "auto")` is read from the frame's `censor`
  column.

* `estimate_betas()` has a single `fmri_frame` method. Frames on a
  `volume_space` return `NeuroVec` betas as the volumetric method did; frames
  on an `index_space` or `basis_space` return coefficient matrices as the
  matrix and latent methods did. The experimental `prewhiten` argument of the
  latent method, which never ran, has been dropped.

* `dataset_spec()` / `instantiate()` / `realize_dataset()` name the new
  constructors: file-backed bindings produce `"nifti_frame"` specs and inline
  matrices produce `"matrix_frame"` specs; the constructor allowlist is
  `matrix_frame`, `nifti_frame`, `neurovec_frame`, and `latent_frame`.

* `simulate_fmri_matrix()$time_series` is now an `fmri_frame`.

* Engine capability `forbid_by_cluster_dataset_classes` is matched against the
  dataset's class or its feature-space class (`"basis_space"` marks latent
  frames).

* `read_fmri_config()` is gone with the legacy API; nothing in `fmrireg`
  used it outside commented-out tests, so it was not ported.

* Fixed the `unused variable 'P'` compiler warning in `ols_t_cpp()`
  (`src/ols_t.cpp`).

* Internally, chunkwise and runwise fitting iterate over lazy frame views
  (`x[rows, ]`, `x[, cols]` + `collect_assay()`) instead of
  `data_chunks()` / `exec_strategy()`, with the same chunk partition, so
  fitted numbers are unchanged.

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

* **Explicit volume weights.** `weights_spec(values = )` now enables weighting
  without requiring a construction method, and joint IID OLS applies the same
  weights to its shared design and every response chunk. Zero-weight volumes
  are excluded from residual degrees of freedom. Reference backends and joint
  AR/robust combinations that cannot apply requested volume weights now fail
  instead of fitting an unweighted model.

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

* **Time-sketched GLM inference (`engine = "latent_sketch"`).** Reported
  statistics now describe the estimator that is actually returned.

  - *Sketch-and-solve standard errors are now honest.* `"srht"`,
    `"gaussian"` and `"countsketch"` fit the model to `m` sketched rows, so
    their coefficients vary about `sqrt(T / m)` times more than full-data
    OLS, but they reported OLS-sized standard errors with `T - p` degrees of
    freedom. Under the null the type-I rate at nominal 0.05 was 0.42 (SRHT),
    0.45 (Gaussian) and 0.45 (CountSketch) at `m = 40`, `T = 200`, and
    0.16-0.23 at `m = 120`. The fits now report the conditional-on-sketch
    covariance `(Xs'Xs)^-1 Xs'SS'Xs (Xs'Xs)^-1`, an unbiased residual
    variance `||r_s||^2 / tr(P SS')` and Satterthwaite residual degrees of
    freedom (about `m - p`), which restores type-I rates of 0.047-0.051 for
    every sketch-and-solve method, with and without AR prewhitening (OLS on
    the same data: 0.052). **Standard errors,
    t-statistics, p-values and `df.residual` from these methods change**:
    standard errors grow by about `sqrt(T / m)`.
  - *Fields are on the data scale.* The SRHT used an unnormalised Hadamard
    transform, so for `"srht"` and `"ihs"` `sigma2`, `result$sigma`,
    `resvar` and `rss` were about `T` times too large and `cov.unscaled`
    about `T` times too small (the two cancelled in standard errors). All
    sketches are now normalised (`E||Sr||^2 = ||r||^2`), and `sigma2` is on
    the scale of the noise variance for every method. For sketch-and-solve
    fits `cov.unscaled` is the conditional covariance above; with landmarks,
    `sigma2` remains the interpolated variance `sum_l w_l^2 sigma2_l`.
  - *IHS is OLS-exact.* `"ihs"` now returns the exact `(X'X)^-1`, the exact
    full-data residual variance and `T - p` degrees of freedom; previously
    its covariance came from one random sketch (scaling the whole SE map by
    a random factor of 0.98-1.56 at `m = 40`) and its variance from a second,
    independent sketch. By default it iterates until no coefficient moves by
    more than `tol = 1e-3` OLS standard errors, up to `iters = 100`
    iterations (median 11, max 16 at `m = 8p` over 200 sketch seeds, leaving
    at most 6e-4 standard errors of error; the former fixed 3 iterations left
    a median of 0.5 and a maximum of 1.9). `time_sketch$tol = 0` runs exactly
    `iters` iterations, and a fit that stops at `iters` before reaching `tol`
    warns. Its standard errors, t-statistics and p-values now match OLS.
  - Fits carry `$sketch` (method, `m`, residual df, IHS iterations and
    convergence). Sketch-and-solve now requires `m > p` and errors
    otherwise.
  - `lowrank_control()` documents all four methods and the `iters` and
    `tol` controls, validates `time_sketch$method`, and no longer carries an
    `iters = 0L` default that was invalid for `"ihs"`.

## Bug Fixes

* Multi-run fits with `baseline_model(intercept = "global")` no longer
  duplicate the global intercept. `term_matrices()` and
  `design_matrix(<fmri_model>, blockid =)` now build the whole design without
  a run filter when every run is selected, so the design holds one
  `constant_global` column instead of one per run, `baseline_term_indices`
  stays within the design, and joint fits no longer alias the duplicate
  columns to `NA`. Selections of a subset of runs are unchanged. Runs are
  identified from the sampling frame, so a run without events still counts
  when every run is requested, and `term_matrices()` now errors if its
  baseline term matrices do not span exactly the baseline design.

* Parcel-pooled AR (`noise_spec(pooling = "parcel")`, `by_cluster`) in
  `engine = "latent_sketch"` summed the sketched Gram matrices of all
  parcels and solved each parcel's cross-products against that sum,
  shrinking every coefficient by about the number of parcels (with 10
  parcels the reported residual variance was 10^4-10^6 times too large and
  no null test ever rejected). Each parcel is now solved with its own whitened design and
  covariance, and `method = "ihs"`, which previously ran plain SRHT here, is
  now honoured. Because the coefficient covariance differs between parcels,
  post-hoc `fit_contrasts()` on such fits errors; declare contrasts in the
  model instead. Parcel labels containing `NA` now error instead of leaving voxels
  without a cluster.
* `engine = "latent_sketch"` is now rank safe. Its solver inverted the
  sketched Gram matrix with `chol()` plus a silent ridge fallback, and
  `chol()` often succeeds on a numerically singular matrix. With one
  aliased nuisance column (`T = 200`, `m = 40`), `"countsketch"` reported
  a mean `sigma2` of 221 against a true 1 and event standard errors of
  38-166 against about 3.6. `"srht"` and `"gaussian"` gave aliased
  coefficients small finite standard errors, and `"ihs"` reported aliased
  variances of about 2e13. Estimability is now judged on the full
  (whitened, per parcel for `by_cluster`) design with the pivoted-QR rule
  that exact fits use. The sketch is solved on the estimable columns only,
  through a QR of the sketched design. As in exact fits, aliased
  coefficients and their standard errors are `NA`, `cov.unscaled` and
  `covariance_by_cluster` carry the `aliased` attribute, and contrasts that
  load on an aliased column are `NA` with a warning naming it. The
  estimable coefficients, standard errors, `sigma2` and Satterthwaite df
  equal those of the same sketch applied to the reduced design. If the
  sketched design restricted to the estimable columns is rank deficient,
  `m` is too small for the design and the fit stops with an error naming
  `time_sketch$m`. The ridge fallback is gone. `"ihs"` likewise no longer
  falls back to a pseudo-inverse of a singular sketched Hessian (with `m`
  below the number of columns it returned coefficients 10-800 OLS
  standard errors off while reporting the exact OLS covariance); it
  requires `m` at least the number of estimable columns and errors
  otherwise. For `by_cluster` fits the non-estimable-contrast warning is
  issued once, not once per cluster.
* Reporting for `engine = "latent_sketch"` fits:
  - `write_results()` metadata and `result$df$method` now record
    `DegreesOfFreedomMethod = "satterthwaite"` for sketch-and-solve fits.
    `"ihs"` and exact fits keep `"residual"`.
  - `by_cluster` fits now report covariance scope `"cluster"` (with the
    per-cluster covariances) instead of `"summary"`.
  - `result$rss` is now the residual sum of squares of the fitted
    (sketched) rows. It was `sigma2 * rdf`, which is not a residual sum of
    squares because `rdf` is a Satterthwaite df. Its expectation is
    `sigma2 * kappa`, where `kappa = tr(P SS')` is stored in
    `result$sketch$kappa`. For `"ihs"`, `rss` is unchanged (the full-data
    RSS). For landmark fits it is interpolated from the landmarks, like
    `sigma2`.
* The single-voxel HRF plot (`autoplot(fit, type = "hrf", voxel = v)`)
  of a `by_cluster` sketch fit drew no confidence band and always marked the
  curve as significant, because it read the shared `cov.unscaled`, which
  such fits do not have. It now uses the covariance of the voxel's cluster.
* `time_sketch` lists without an `m` element (for example
  `list(method = "ihs")`) failed with "m <= Tlen is not TRUE", because
  `sk$m` partially matched `sk$method`. The default `m = min(8p, T)` now
  applies.

* The `"ihs"` time sketch in `engine = "latent_sketch"` now performs an
  actual iterative Hessian sketch. Each iteration previously used a sketched
  gradient as well as a sketched Hessian, so it re-solved an independent
  sketched problem and never converged to the least-squares solution; more
  iterations did not help. The gradient now uses the full data, and a step
  halving keeps the residual sum of squares non-increasing, so `iters`
  controls accuracy as documented. Iterations start from the sketch-and-solve
  solution, so the baseline is fitted from the first step. `iters < 1` and
  non-finite inputs are now errors.
* The `fmridataset` requirement is now `(>= 0.11.0.9000)`. The pre-frame API
  was removed upstream without a version change, so both sides of the break
  reported `0.10.0.9000` and the previous constraint could not tell them
  apart: a stale install resolved cleanly and then failed at lazy loading with
  `object 'as.matrix_dataset' is not exported by 'namespace:fmridataset'`.
  fmridataset bumped to `0.11.0.9000` for this purpose
  (bbuchsbaum/fmridataset#92), so the mismatch is now caught at dependency
  resolution with a message that names the real problem.

* The test suite has been ported to the frame API. The frame migration and a
  parallel test-coverage branch were developed against the same base and
  merged independently, so the merged tree combined frame-only package code
  with ~50 test files still building legacy `matrix_dataset()`,
  `fmri_mem_dataset()`, and `latent_dataset()` fixtures. `R CMD check` failed
  on every platform with 106 test failures that neither branch saw on its own.
  Fixtures now use `matrix_frame()`, `neurovec_frame()`, and `latent_frame()`;
  no assertion was relaxed to accommodate the port.

* `chunkwise_lm()` dispatches on `fmri_frame` from outside the package.
  `chunkwise_lm.fmri_frame()` was never registered as an S3 method, so
  dispatch from outside the namespace failed with "no applicable method"
  even though internal calls resolved.

* `.rrr_extract_response_matrix()` validates its input again. A non-frame
  argument fell through to `collect_assay()` and failed inside `fmridataset`
  with a message naming neither the argument nor the caller.

* `extract_censor_from_dataset()` no longer errors on a non-frame dataset.
  It now reports "no censoring" instead, which matters because censoring is
  consulted on every fit under the default change below.

* A censor column carried by a dataset is now used by default. `fmri_lm()`
  previously consulted it only when `ar_options = list(censor = "auto")` was
  passed as well, so a `matrix_frame(censor = )` (and, before it,
  `fmri_dataset(censor = )`) was silently discarded: fits on a censored frame
  were bit-identical to fits on an uncensored one, with no warning. An unset
  `censor` now resolves against the dataset, which is what `"auto"` asked for
  explicitly; pass `censor = "none"` to ignore a censor column deliberately.
  **This changes results** for anyone who set censoring on a dataset and did
  not opt in, since those fits were not censored at all. An explicit `censor`
  vector in `ar_options` still overrides the dataset.

* `fmri_lm()` now warns when censoring cannot affect the fit. Censoring feeds
  AR estimation and whitening only and never drops volumes from the
  regression, so under the default iid noise model it does nothing at all.
  That was previously silent, and indistinguishable from censoring having
  been applied.

* Shared AR estimation now pools voxel residual autocovariances by default
  instead of fitting the cross-voxel mean residual series. The former targets
  a typical voxel covariance; the latter suppresses voxel-specific noise and
  targets a different, often more autocorrelated coherent component. The old
  behavior remains available explicitly as
  `noise_spec(..., shared_estimator = "mean_series")`. Fixed-order fits now use
  fmriAR's pooled AR autocovariance path as well; the previous ARMA(p, 0)
  shortcut accidentally collapsed both choices to the mean series.
* `shared_estimator = "mean_series"` now retains the matching OLS design-bias
  correction. OLS projection is linear across response columns, so averaging
  residuals after projection is exactly the same as projecting the mean series.
* The OLS residual-bias correction now uses a separate covariance-tail lag
  budget rather than stopping at the fitted AR order. Projection mixes tail
  autocovariances into lags 0:p, so a p-lag correction left avoidable bias.
* Global AR estimation over separately fitted runs now supplies the exact
  block-diagonal residual-forming design. Row-binding the per-run designs
  incorrectly described a single coefficient vector shared across runs.
* Design-corrected AR coefficients are estimated once from initial OLS
  residuals and held fixed for the GLS solve. Later GLS residuals have a
  different residual-forming operator and no longer reuse the OLS `design=`
  correction. The low-rank engine now follows this contract too, and the
  low-rank and RRR whitening paths reset at run boundaries.
* GDS exports now release unreachable HDF5 handles before atomic finalization,
  allowing Windows to move and immediately reopen the completed file.

* `write_results(strategy = "by_stat")` no longer scatters contrast maps onto
  the wrong voxels. `.compute_statistical_volumes()` called `as.logical()` on
  the mask, which drops `dim` for plain arrays (the type
  `.fmri_dataset_mask_space()` always supplies), so `which(..., arr.ind = TRUE)`
  returned linear indices instead of 3D coordinates.
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

## Deprecated

* `lowrank_control(time_sketch = list(method = "ihs"))` is deprecated. In
  `engine = "latent_sketch"` the design is narrow (p of about 10-30) and the
  response wide (1e4-1e5 columns), so least squares costs one `X'Z` pass,
  and every iterative-Hessian-sketch iteration repeats that pass to form the
  exact gradient: IHS could never beat one exact solve (9.35 s against
  0.12-0.23 s for OLS). `"ihs"` now computes the exact OLS fit through the
  solver of the exact joint path, so coefficients, standard errors,
  contrasts, rank handling (aliased columns `NA`) and the `T - rank`
  residual df are those of exact OLS; without AR it is identical to
  `fmri_lm(..., strategy = "chunkwise")`. It emits a once-per-session
  message pointing to `"countsketch"`/`"gaussian"` for speed or exact OLS
  for accuracy. This supersedes the IHS iteration and `tol` entries under
  Statistical Corrections and Bug Fixes. The controls `time_sketch$iters`
  and `time_sketch$tol` are accepted and ignored for one release, `m` is
  ignored for `"ihs"` (it no longer has to reach the number of estimable
  columns), and `result$sketch` no longer carries `iters` or `converged`.
  The internal IHS kernel (`cpp_ihs_latent()`, `ihs_latent_solve()`) is
  removed.

## Performance

* The SRHT sketched solve in `engine = "latent_sketch"` is 17-83x faster
  (`T` from 200 to 1000, 1e4-5e4 response columns, 2 threads). Its
  Walsh-Hadamard transform walked the column-major data row by row with a
  temporary vector per butterfly; it now transforms each response column in
  place in a contiguous buffer and, when built with OpenMP, runs in parallel
  over columns (`options(fmrireg.num_threads)`). Output is bit-identical to
  the previous kernel. The sketched solve now forms coefficients and
  residuals with two BLAS products through the explicit `Q` factor instead
  of `qr.coef()`/`qr.resid()`, which applied Householder reflections one
  response column at a time, and the sketch Gram matrix `SS'` is computed
  once per fit instead of once per cluster. At `T = 400`, `p = 20`,
  `m = 8p` and 5e4 response columns the SRHT solve fell from 3.5 s to
  0.07 s, and the Gaussian and CountSketch solves from about 0.3 s to
  0.04 s (exact OLS: 0.05 s). SRHT still trails CountSketch and Gaussian
  at realistic sizes, so `?lowrank_control` now recommends those two for
  speed.
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
