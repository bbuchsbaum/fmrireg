# Summary of Namespace Fixes After fmrihrf Migration

## Latest Update: Fixed `evaluate` and `nbasis` namespace conflicts

## Problem
After migrating HRF functionality to the fmrihrf package, tests were failing with errors like:
- `no applicable method for 'blocklens' applied to an object of class "sampling_frame"`
- `no applicable method for 'samples' applied to an object of class "sampling_frame"`
- `no applicable method for 'nbasis' applied to an object of class "c('HRF', 'function')"`

## Root Cause
Both fmrireg and fmrihrf were defining the same generic functions, causing namespace conflicts. When fmrireg code called these functions without explicit namespace qualification, R couldn't determine which version to use.

## Solution
1. **Commented out generic definitions in fmrireg** that are now provided by fmrihrf:
   - `blocklens()`, `blockids()`, `samples()`, `global_onsets()` in `R/all_generic.R`

2. **Updated imports and re-exports** in `R/fmrihrf-imports.R`:
   - Already importing these generics from fmrihrf
   - Added re-exports so users can use them without qualification

3. **Fixed unqualified function calls in tests**:
   - `test_baseline.R`: Fixed `blocklens()` and `samples()` calls
   - `test_event_model.R`: Fixed `blocklens()` and `nbasis()` calls
   - `test_condition_basis_list.R`: Fixed `blocklens()` call
   - `test_sampling_frame.R`: Fixed all `sampling_frame()`, `samples()`, and `global_onsets()` calls

## Results
- ✅ `test_baseline.R`: All 22 tests passing
- ✅ `test_sampling_frame.R`: All 31 tests passing
- ✅ `test_event_model.R`: Most tests passing (113 pass, 4 fail due to other issues)
- ✅ `test_benchmark_datasets.R`: All 67 tests passing

## Additional Fixes (Latest)
1. **Commented out `evaluate` generic** in `R/all_generic.R` (line 909-911)
2. **Added re-exports** in `R/fmrihrf-imports.R`:
   - `fmrihrf::evaluate`
   - `fmrihrf::nbasis`
3. **Fixed unqualified function calls**:
   - `R/benchmark_datasets.R`: Fixed `evaluate()` call
   - `R/simulate.R`: Fixed `evaluate()` and `nbasis()` calls
   - `R/fmri_lm_methods.R`: Fixed `evaluate()` call
   - `R/hrf-formula.R`: Fixed `evaluate()` call in `evaluate.hrfspec`
   - `R/autoplot-methods.R`: Fixed `evaluate()` and `nbasis()` calls
   - `R/event_vector.R`: Fixed `nbasis()` calls
   - Test files: Fixed `nbasis()` calls in `test_event_model.R`

## Remaining Issues
1. Some tests in `test-regressor.R` are testing old internal functions that have been replaced
2. A few tests in `test_event_model.R` fail due to issues with Fourier basis (not namespace related)

## Recommendations
1. Review and update tests that use old internal functions (like `Reg()` instead of `regressor()`)
2. Consider removing test files that duplicate fmrihrf functionality testing
3. Ensure all new code uses `fmrihrf::` qualification for clarity when calling fmrihrf functions