# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Build and Development Commands

### Package Building and Checking

- `devtools::document()` - Generate/update documentation from roxygen2
  comments
- `devtools::check(cran = TRUE)` - Run comprehensive CRAN checks (this
  is the gold standard)
- `devtools::spell_check()` - Check spelling in documentation and
  DESCRIPTION
- `urlchecker::url_check()` - Validate all URLs in documentation
- `devtools::check_win_devel()` - Check on Windows development version
- `rhub::check_for_cran()` - Run multi-platform CRAN checks

### Testing

- `devtools::test()` - Run all tests
- [`testthat::test_local()`](https://testthat.r-lib.org/reference/test_package.html) -
  Run tests locally
- `testthat::test_file("tests/testthat/test_specific.R")` - Run a single
  test file
- Tests are located in `tests/testthat/`

### Documentation

- [`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html) -
  Build the package website
- Documentation follows CRAN guidelines strictly (see CRAN_guidance.md)

## High-Level Architecture

fmrireg is an R package for fMRI time series regression analysis with a
layered architecture:

### Core Abstractions

1.  **Event System** (`R/event-classes.R`, `R/event_vector.R`)
    - `event_factor`: Categorical experimental events
    - `event_variable`: Continuous experimental variables
    - `event_matrix`: Matrix-based events
    - Events are convolved with HRFs to create regressors
2.  **HRF System** (`R/hrf.R`, `R/hrf-*.R`)
    - Base `HRF` class with various implementations (gamma, gaussian,
      spline, etc.)
    - HRFs can be modified via decorators (lag, block, normalize)
    - Supports custom basis sets and parametric variations
3.  **Model Hierarchy**
    - `baseline_model`: Nuisance regressors (drift, motion, etc.)
    - `event_model`: Experimental design (events + HRFs)
    - `fmri_model`: Complete model (baseline + events)
4.  **Dataset Abstractions** (`R/fmri_dataset.R`)
    - `fmri_dataset`: Base class for fMRI data
    - `matrix_dataset`: In-memory matrix data
    - `fmri_mem_dataset`: Memory-mapped fMRI data
    - `latent_dataset`: Reduced dimensionality data
5.  **Model Fitting** (`R/fmri_model.R`, `R/fmrilm.R`)
    - `fmri_lm`: Standard GLM fitting
    - `fmri_rlm`: Robust GLM fitting
    - Supports different strategies: runwise, chunkwise, trial-wise
    - C++ implementations for performance (mixed_solve, AR whitening)
6.  **Contrast System** (`R/contrast.R`)
    - Flexible contrast specification via formulas
    - Support for F-contrasts, pairwise, polynomial contrasts

### Key Design Patterns

- **S3 Object System**: All generics defined in `R/all_generic.R`
- **Builder Pattern**: Models built incrementally (events → event_model
  → fmri_model)
- **Strategy Pattern**: Different fitting algorithms (OLS, robust,
  regularized)
- **Decorator Pattern**: HRF modifications (lag_hrf, block_hrf)

### Performance Considerations

- C++ implementations via Rcpp for computationally intensive operations
- Parallelization has two layers: R-level loops go through
  `.fmrireg_lapply()` (R/parallel.R), which uses future.apply under the
  user’s
  [`future::plan()`](https://future.futureverse.org/reference/plan.html);
  C++ kernels use OpenMP threads from `.fmrireg_num_threads()`
  (`options(fmrireg.num_threads)` / `FMRIREG_NUM_THREADS`). Do not add
  another parallel backend (no furrr, RcppParallel, parallel::mclapply).
- Chunked processing for large datasets to manage memory

### CRAN Compliance

The package strictly follows CRAN guidelines: - S3 generics are the
primary documentation site - Methods use `@rdname` to link to generic
documentation - All examples must run quickly (\< 5 seconds) - Avoid
modifying user options or par() without restoration - See
`CRAN_guidance.md` for detailed guidelines

### Development Principles

From `data-raw/principles.md`: - Single source of truth for data
representations - Encapsulation over reconstruction - Clear separation
of concerns - Functional composition preferred - Object-oriented design
with S3 - One way to do one thing - Fail fast and locally with clear
error messages

## AR MODELING MIGRATION TO fmriAR (2025-09-18)

### Migration Complete

All AR (autoregressive) modeling functionality has been migrated to the
specialized fmriAR package:

1.  **Dependencies Updated**: ✓
    - Added fmriAR to Imports and Remotes in DESCRIPTION
    - fmriAR provides enhanced AR/ARMA modeling capabilities
2.  **Integration Adapter**: ✓
    - Created `R/fmriAR_adapter.R` for seamless integration
    - Maps fmrireg configurations to fmriAR parameters
    - Maintains backward compatibility
3.  **Legacy Code Removed**: ✓
    - Removed `R/fmri_ar_modeling.R`, `R/ar_utils.R`
    - Removed `src/ar_whiten.cpp` (C++ implementation)
    - Created compatibility layer for legacy function calls
4.  **Benefits**:
    - Access to advanced features: ARMA models, multiscale pooling,
      better diagnostics
    - Improved performance via optimized C++ implementation
    - Reduced maintenance burden (~600 lines of code removed)
    - Single source of truth for AR functionality

### API Compatibility

All existing
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
calls continue to work unchanged: - `cor_struct = "ar1"` → delegates to
fmriAR with p=1 - `cor_struct = "ar2"` → delegates to fmriAR with p=2 -
`cor_struct = "arp"` → delegates to fmriAR with user-specified p - All
options preserved: `exact_first`, `iter_gls`, `global`, etc.

### Testing

- Integration tests in `test_fmriAR_integration.R`
- All existing AR tests continue to pass
- Numerical results match within tolerance

## CRITICAL TEST ISSUES RESOLVED (2025-05-28)

All critical functionality tests are now passing after the following
fixes:

1.  **Config structure mismatch**: FIXED ✓
    - Updated `solve_integrated_glm` to handle both config structures
      (`config$ar$struct` and `config$ar_options$cor_struct`)
    - Added proper config extraction in all pipeline functions
    - Files modified: `fmri_lm_integrated_solver.R`
2.  **Effective df for robust models**: FIXED ✓
    - Corrected the formula to use sum of weights approach:
      `df = sum(weights) - p`
    - File modified: `fmri_lm_effective_df.R`
3.  **Missing XtXinv in results**: FIXED ✓
    - Added XtXinv to all pipeline results for contrast computation
    - Added dfres to core solver results
    - Files modified: `fmri_lm_integrated_solver.R`, `fmri_lm_solver.R`

All 6 critical tests now pass: - ✓ AR whitening integration (now via
fmriAR) - ✓ Robust fitting - ✓ Contrast computation - ✓ Effective df
calculations - ✓ Bootstrap functionality - ✓ Sandwich variance estimator

The refactored code is now ready for integration with the existing
codebase.

### Correction (2026-08-07)

The “sandwich variance estimator” item above overstated what shipped.
Both sandwich helpers (`compute_sandwich_variance`,
`calculate_sandwich_variance`) are reachable only by direct call — no
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
path invokes them, and the robust path reports a model-based WLS
variance instead. The test that covered them asserted only dimensions
and positivity, which is why an arithmetic error (`X'diag(e^4)X` in
place of `X'diag(e^2)X`, giving standard errors 3.5x too large) survived
in `compute_sandwich_variance` until 2026-08-07. That error is now fixed
and both helpers are checked against
[`sandwich::vcovHC`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html)
in `tests/testthat/test-variance-closed-form.R`; they remain unwired.
(2026-10: both helpers were removed as unreachable.
`variance_spec("sandwich")` is the live sandwich path; see
`data-raw/simplification-2026-10-plan.md`.)

## REVIEW 2026-10: CORRECTNESS AND PERFORMANCE PASS

Plan, checklist and results: `data-raw/review-2026-10-plan.md`. Every
fix has a regression test that fails on the pre-review code; the tests
live in `tests/testthat/test-review-{A,B,C,D,E,F}-*.R` and compare
against `lm`, `anova`, `t.test`, `sandwich`,
[`metafor::rma`](https://wviechtb.github.io/metafor/reference/rma.uni.html),
closed forms, or pure-noise calibration simulations. Keep that standard:
a test that only checks shapes or `is.finite` is how most of these bugs
survived.

Conventions established by this pass (rely on them, keep them true): -
**Residual scale is per voxel everywhere.** `result$sigma`, `rss` and
`resvar` have one value per voxel for every strategy (multi-run runwise
fits stored one pooled scalar until this pass). - **Robust fits are
reproducible WLS.** Every robust path reports the weights used for the
final betas, sigma^2 = RSS_w / (sum(w \> 0) - rank), and the effective
df `.robust_effective_df()`: sum(w) - rank. There is one effective-df
implementation (`R/effective_df.R`). - **Stored covariance matches
stored scale.** `cov.unscaled` / `covariance_by_voxel` are in the fitted
(whitened/weighted) domain, so post-hoc
[`fit_contrasts()`](https://bbuchsbaum.github.io/fmrireg/reference/fit_contrasts.md)
equals in-fit contrasts. Multi-run runwise fits that whiten or weight
per run have no such covariance;
[`fit_contrasts()`](https://bbuchsbaum.github.io/fmrireg/reference/fit_contrasts.md)
refuses them (`.posthoc_check_supported`). - **Use `[[` for result
fields.** `result$sigma` partially matches `sigma2`; that bug inflated
reference-backend AR standard errors ~3x. - **Group-analysis sign.**
`fmri_ttest` contrasts are computed once, in the requested `sign`;
`.fmri_ttest_apply_group_sign` flips coefficient rows only. -
**Meta-analysis weights.** `weights = "equal"/"custom"` are relative
weights with an estimated residual scale (WLS); `"ivw"` treats V as
known sampling variances. `method = "reml"` runs Paule-Mandel and says
so. - **Spatial FDR** estimates its empirical null by Efron’s central
matching and numbers groups by sorted label / factor level. Since P1
(below) the computation is
[`fmrigds::spatial_fdr_sabha()`](https://rdrr.io/pkg/fmrigds/man/spatial_fdr_sabha.html);
[`spatial_fdr()`](https://bbuchsbaum.github.io/fmrireg/reference/spatial_fdr.md)
wraps it. - **Chunks are lazy.** `.dset_chunks()` returns descriptors;
data are read in the worker via `.chunk_data()`; run rows come from
`.dset_run_rows()`. - **`.fast_preproject()` never forms an n x n
matrix**; `proj$Pinv` is computed on demand.

Corrections to the 2025 notes above: the robust effective df is sum(w) -
rank (the “sum(weights) - p” fix lived in an unreachable file until this
pass); “Bootstrap functionality” referred to
`multiresponse_bootstrap_lm`, which was unreachable and has been removed
(`bootstrap_glm_inference`, also never wired into
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md),
was removed in the 2026-10 simplification pass).

## SIMPLIFICATION 2026-10

Plan, bugs found and remaining phases:
`data-raw/simplification-2026-10-plan.md`. Phase 1 removed code with no
caller reachable from an export, S3 method or `.onLoad`. Before adding a
helper, check it will have a live caller; before keeping one “for
later”, note why in the plan.

Group analysis is moving to fmrigds (plan Phase 4 item 1). Since P1:
fmrireg registers exactly one fmrigds reducer, `"fmrireg:meta_reg"` (the
`.meta_fit_core()` estimator), plus the post-hoc alias
`"fmrireg:fdr:spatial"` (= fmrigds `"fdr:sabha"`), at load and never
over an fmrigds name. The gds path of
[`fmri_meta()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md)
is reduce -\>
[`fmrigds::coef_contrast()`](https://rdrr.io/pkg/fmrigds/man/coef_contrast.html)
-\> wrap; since P2,
[`fmri_ttest()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_ttest.md)
normalises every input with `.fmri_ttest_input()` and its classic /
Welch engines call fmrigds’s `"ols:voxelwise"` / `"ttest:welch"` kernels
(`.fmri_ttest_ols()`, `.fmri_ttest_welch()`); reducer options are read
with `[[` (fmrigds’s `options$weights` partially matches a
`weights_custom` element, hence the option name `custom_weights`).
