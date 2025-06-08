# fmrireg to fmrihrf Migration Summary

## Date: 2025-06-07

## Overview
Successfully migrated fmrireg to use the fmrihrf package for core HRF functionality, removing code duplication while maintaining backward compatibility.

## Changes Made

### 1. Added Dependencies
- Added `fmrihrf (>= 0.1.0)` to Imports in DESCRIPTION
- Added `bbuchsbaum/fmrihrf` to Remotes (already present)

### 2. Removed Redundant Files
The following files were removed as their functionality is now provided by fmrihrf:
- `R/hrf.R` - Core HRF class
- `R/hrf-functions.R` - HRF implementations  
- `R/hrf-afni.R` - AFNI HRF functions
- `R/hrf_decorators.R` - HRF decorators
- `R/hrf_from_coefficients.R` - HRF from coefficients
- `R/reg-constructor.R` - Regressor construction
- `R/reg-methods.R` - Regressor methods
- `R/regressor.R` - Consolidated regressor file
- `R/sampling_frame.R` - Temporal structure
- `R/penalty_matrix.R` - Penalty matrices
- `R/reconstruction_matrix.R` - Reconstruction matrix
- `R/evaluate-helpers.R` - Evaluation helpers

### 3. Kept Regression-Specific Files
- `R/hrf-formula.R` - Contains critical `hrf()` formula function for regression models
- `R/hrf_smoothing_kernel.R` - Creates temporal similarity matrices for regression

### 4. Updated Import Statements
Updated 42 source files and 26 test files to use `fmrihrf::` namespace for HRF functions:
- Heavy users: event_model.R, event_model_helpers.R, afni.R, simulate.R
- Medium users: benchmark_datasets.R, fmri_model.R, baseline_model.R
- Light users: Various files with 1-5 HRF references
- Test files: Comprehensive test suite updates

### 5. Created Import/Export File
Created `R/fmrihrf-imports.R` to:
- Import all necessary functions from fmrihrf
- Re-export commonly used functions for backward compatibility

### 6. Fixed Issues
- Removed references to non-existent `getHRF()` function
- Updated hrf-formula.R to use HRF constants directly
- Updated DESCRIPTION Collate section to remove deleted files
- Fixed generic function definitions in all_generic.R

## Benefits Achieved

1. **Reduced Code Duplication**: Removed ~12 files worth of duplicated HRF code
2. **Cleaner Architecture**: Clear separation between HRF functionality and regression
3. **Easier Maintenance**: Single source of truth for HRF implementations
4. **Backward Compatible**: Re-exported functions ensure existing code continues to work
5. **Better Modularity**: HRF functionality can be used independently via fmrihrf

## Testing Status

- Package loads successfully
- HRF constants accessible (e.g., `HRF_SPMG1`)
- Basic HRF functionality working through fmrihrf
- Full test suite needs to be run to verify complete functionality

## Next Steps

1. Run full test suite: `devtools::test()`
2. Run CRAN checks: `devtools::check(cran = TRUE)`
3. Update any vignettes that reference HRF internals
4. Consider adding migration notes to NEWS.md
5. Update CI/CD pipelines to install fmrihrf dependency

## Rollback Instructions

If needed, the migration can be rolled back:
```bash
git checkout fmrihrf-migration-backup
```

This will restore the pre-migration state with all original HRF files.