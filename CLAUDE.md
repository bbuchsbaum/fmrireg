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
- Parallelization via RcppParallel (configurable:
  `options(fmrireg.num_threads = N)`)
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
