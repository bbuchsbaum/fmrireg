# Critical Test Coverage Gaps in fmrireg

## Summary of Test Coverage Analysis

### 1. **Bootstrap Functionality** (CRITICAL GAP)
- **Status**: No existing tests
- **Issues Found**: 
  - `multiresponse_bootstrap_lm` has structural issues with contrast handling
  - Missing column 'estimate' in results
  - Covariance calculation fails
- **Impact**: Bootstrap confidence intervals are untested and likely broken

### 2. **Effective Degrees of Freedom** (CRITICAL GAP)
- **Status**: No existing tests
- **Issues Found**:
  - Functions `calculate_effective_df` and `sandwich_var_adjustment` don't exist
  - Core statistical inference functionality is missing
- **Impact**: P-values and statistical inference may be incorrect

### 3. **Solver Edge Cases** (CRITICAL BUGS)
- **Status**: Limited testing
- **Issues Found**:
  - Cannot handle rank-deficient matrices (crashes with Cholesky decomposition)
  - Robust fitting missing required argument `X_orig_for_resid`
  - `mixed_solve_cpp` has index out of bounds errors
- **Impact**: Common real-world scenarios will crash the solver

### 4. **AR Modeling** (CRITICAL BUG)
- **Status**: Some tests exist but incomplete
- **Issues Found**:
  - `ar_whiten_inplace` not applying correct transformation
  - First row scaling incorrect
  - Subsequent row differencing incorrect
- **Impact**: AR correction is broken, leading to incorrect standard errors

### 5. **Latent Models** (UNTESTED)
- **Status**: Minimal coverage
- **Issues Found**:
  - API mismatches with `matrix_dataset`
  - `chunkwise_lm.latent_dataset` untested
- **Impact**: Dimensionality reduction workflows are untested

### 6. **Integration Scenarios** (UNTESTED)
- **Status**: No integration tests
- **Critical Combinations**:
  - AR + Robust fitting
  - Bootstrap + Robust fitting
  - Multi-run with different strategies
  - Missing data handling
- **Impact**: Complex workflows may fail in production

## Most Critical Bugs to Fix

1. **Rank-deficient matrix handling** - Common in real fMRI data
2. **AR whitening implementation** - Core functionality is broken
3. **Robust fitting argument mismatch** - Prevents robust regression
4. **Mixed solve indexing** - Causes crashes

## Recommendations

1. **Fix critical bugs first** before adding more tests
2. **Add integration tests** for common workflows
3. **Test numerical edge cases** systematically
4. **Ensure bootstrap functionality** works for confidence intervals
5. **Validate statistical calculations** against known results

## Test Files Created

1. `test_bootstrap.R` - Reveals structural issues in bootstrap implementation
2. `test_effective_df.R` - Shows missing statistical functions
3. `test_solver_edge_cases.R` - Exposes numerical instability
4. `test_fmri_latent_lm.R` - Has syntax errors, needs fixing
5. `test_integration_scenarios.R` - Has syntax errors, needs fixing
6. `test_critical_gaps.R` - Working tests that expose real bugs