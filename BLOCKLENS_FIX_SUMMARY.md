# Summary of Fixes for fmrihrf Integration Issues

## Problem
After migrating to use fmrihrf, tests were failing with errors like:
- `no applicable method for 'blocklens' applied to an object of class "sampling_frame"`
- `no applicable method for 'blockids' applied to an object of class "sampling_frame"`
- `no applicable method for 'nbasis' applied to an object of class "c('HRF', 'function')"`
- `could not find function "Reg"`

## Root Cause
Both fmrireg and fmrihrf export the same generic functions (blocklens, blockids, nbasis, evaluate, etc.), causing namespace conflicts. When fmrireg code called these functions without explicit namespace qualification, R couldn't determine which version to use.

## Solution
Updated all unqualified function calls to use explicit `fmrihrf::` namespace qualification.

### Files Modified

1. **Direct field access replaced with function calls:**
   - `afni.R`: `$blocklens` → `fmrihrf::blocklens()`
   - `covariate.R`: `$blocklens` → `fmrihrf::blocklens()`
   - `event_model_helpers.R`: `$blocklens` → `fmrihrf::blocklens()`
   - `event_vector.R`: `$blocklens` → `fmrihrf::blocklens()`
   - `fmri_dataset.R`: `$blocklens` → `fmrihrf::blocklens()`
   - `fmri_model.R`: `$blocklens` → `fmrihrf::blocklens()`
   - `fmrilm.R`: `$blocklens` → `fmrihrf::blocklens()`

2. **Unqualified function calls updated:**
   - `blocklens()` → `fmrihrf::blocklens()` (7 locations)
   - `blockids()` → `fmrihrf::blockids()` (13 locations)
   - `samples()` → `fmrihrf::samples()` (2 locations)
   - `global_onsets()` → `fmrihrf::global_onsets()` (3 locations)
   - `nbasis()` → `fmrihrf::nbasis()` (4 locations)
   - `evaluate()` → `fmrihrf::evaluate()` (2 locations)
   - `Reg()` → `fmrihrf::regressor()` (2 locations)

3. **Import statements updated:**
   - Added imports for `blocklens`, `blockids`, `samples`, `global_onsets` in `fmrihrf-imports.R`

## Result
All AR GLM integration tests now pass successfully (8 tests passed, 0 failed).

## Lessons Learned
1. When two packages export the same generics, explicit namespace qualification is essential
2. Direct field access (e.g., `obj$field`) bypasses S3 method dispatch and should be replaced with accessor functions
3. Migration scripts should check for both direct field access and unqualified function calls

## Recommendation
Consider adding a linter rule or check to ensure all calls to fmrihrf functions use explicit namespace qualification to prevent future issues.