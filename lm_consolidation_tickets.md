# LM Consolidation Sprint Tickets

## Critical Path (Must Complete in Order)

### CONS-001: Backup and Disable Conflicting Implementation ✓
**Priority**: CRITICAL
**Effort**: 15 min
**Status**: COMPLETE
**Tasks**:
- [x] Rename `fmrilm_new.R` to `fmrilm_new.R.bak`
- [x] Remove `fmrilm_new.R` from DESCRIPTION Collate field  
- [x] Run `devtools::document()` to update NAMESPACE
- [x] Verify `fmri_lm` from `fmrilm.R` is now active
- [x] Run tests to establish new baseline failure count
**Result**: Test failures reduced from 132 to 12!

### CONS-002: Inventory Missing Functions ✓
**Priority**: CRITICAL  
**Effort**: 30 min
**Status**: COMPLETE
**Tasks**:
- [x] List all functions in `fmrilm.R` 
- [x] List all functions in modular files
- [x] Identify functions used by tests but missing from modular implementation
- [x] Create mapping of functions to their intended modular homes
- [x] Document in `lm_consolidation_mapping.md`
**Result**: Identified key gaps - fmri_lm_fit doesn't use modular solvers, AR parameter mapping missing

### CONS-003: Restore AR Parameter Compatibility ✓
**Priority**: HIGH
**Effort**: 1 hour
**Status**: COMPLETE
**Tasks**:
- [x] Add AR parameter handling to `fmrilm.R`'s `fmri_lm` function
- [x] Map old parameters (cor_struct, cor_iter, etc.) to new ar_options structure
- [x] Ensure backward compatibility for existing test code
- [x] Test with `devtools::test(filter='ar_')`
**Result**: AR parameter mapping already existed in fmri_lm. Fixed config attributes and glm_context issues.

### CONS-004: Connect fmri_lm to Integrated Solver ✓
**Priority**: HIGH
**Effort**: 2 hours
**Status**: COMPLETE
**Tasks**:
- [x] Update `fmri_lm_fit` in `fmrilm.R` to detect when AR+Robust is requested
- [x] Route AR+Robust requests to `solve_integrated_glm` 
- [x] Ensure config object is properly created and passed
- [x] Handle the "whiten then robustly weight" pipeline
- [x] Test with AR+Robust combination examples
**Result**: Created `process_run_integrated` wrapper and connected to runwise_lm in fast path

### CONS-005: Implement runwise_lm Using Modular Components
**Priority**: HIGH
**Effort**: 2 hours  
**Blocked by**: CONS-004
**Tasks**:
- [ ] Update `runwise_lm` to use `solve_integrated_glm` for each run
- [ ] Ensure proper handling of run indices
- [ ] Aggregate results across runs appropriately
- [ ] Maintain backward compatibility with existing return structure
- [ ] Test with `devtools::test(filter='runwise')`

### CONS-006: Implement chunkwise_lm Using Modular Components
**Priority**: HIGH
**Effort**: 2 hours
**Blocked by**: CONS-004  
**Tasks**:
- [ ] Update `chunkwise_lm` to use modular solver
- [ ] Ensure proper chunking of voxels
- [ ] Handle contrast computation per chunk
- [ ] Aggregate results appropriately
- [ ] Test with `devtools::test(filter='chunkwise')`

### CONS-007: Fix fast_lm Path
**Priority**: MEDIUM
**Effort**: 1 hour
**Blocked by**: CONS-004
**Tasks**:
- [ ] Identify what `use_fast_path=TRUE` is supposed to do
- [ ] Update fast path to use modular components
- [ ] Ensure it works with AR options
- [ ] Test with `use_fast_path=TRUE` examples

### CONS-008: Migrate Helper Functions
**Priority**: MEDIUM
**Effort**: 2 hours
**Blocked by**: CONS-005, CONS-006
**Tasks**:
- [ ] Move `.fast_preproject` to `fmri_lm_internal.R` (if not already there)
- [ ] Move `.fast_lm_matrix` to appropriate module
- [ ] Move `fast_rlm_run` to `fmri_robust_fitting.R`
- [ ] Update all references to use new locations
- [ ] Remove duplicates

### CONS-009: Consolidate S3 Methods
**Priority**: MEDIUM
**Effort**: 1 hour
**Blocked by**: CONS-008
**Tasks**:
- [ ] Ensure all S3 methods (coef, fitted_hrf, stats, etc.) work with new structure
- [ ] Update methods to use modular result structure
- [ ] Test each S3 method individually
- [ ] Update documentation

### CONS-010: Remove Duplicate Implementations
**Priority**: LOW
**Effort**: 1 hour
**Blocked by**: CONS-009
**Tasks**:
- [ ] Remove duplicate functions from `fmrilm.R`
- [ ] Clean up commented code
- [ ] Ensure `fmrilm.R` only contains public API
- [ ] Update internal function calls to use modular versions

## Testing & Validation

### CONS-011: Create Integration Tests
**Priority**: HIGH
**Effort**: 2 hours
**Blocked by**: CONS-006
**Tasks**:
- [ ] Create test for basic OLS fitting
- [ ] Create test for AR(1) fitting
- [ ] Create test for Robust fitting
- [ ] Create test for AR+Robust "whiten then robustly weight"
- [ ] Create test for contrasts with all methods
- [ ] Ensure tests cover both runwise and chunkwise strategies

### CONS-012: Fix Remaining Test Failures
**Priority**: HIGH
**Effort**: 4 hours
**Blocked by**: CONS-011
**Tasks**:
- [ ] Run full test suite
- [ ] Categorize remaining failures
- [ ] Fix failures related to consolidation
- [ ] Document any tests that need updating
- [ ] Achieve 0 test failures

## Documentation

### CONS-013: Update Documentation
**Priority**: MEDIUM
**Effort**: 2 hours
**Blocked by**: CONS-012
**Tasks**:
- [ ] Update `fmri_lm` documentation with all parameters
- [ ] Document the AR+Robust pipeline
- [ ] Add examples showing different use cases
- [ ] Update any vignettes that use `fmri_lm`
- [ ] Create developer documentation for modular architecture

### CONS-014: Clean Up and Archive
**Priority**: LOW
**Effort**: 30 min
**Blocked by**: CONS-013
**Tasks**:
- [ ] Archive `fmrilm_new.R.bak` to a deprecated folder
- [ ] Remove temporary consolidation files
- [ ] Update DESCRIPTION if needed
- [ ] Create final consolidation report

## Verification Checklist

### CONS-015: Final Validation
**Priority**: CRITICAL
**Effort**: 1 hour
**Blocked by**: CONS-014
**Tasks**:
- [ ] All tests pass (0 failures)
- [ ] AR+Robust integration works as specified
- [ ] Performance benchmarks are acceptable
- [ ] Code follows modular architecture principles
- [ ] No duplicate implementations remain
- [ ] Documentation is complete and accurate

---

## Summary

**Total Tickets**: 15
**Estimated Total Effort**: ~25 hours
**Critical Path**: CONS-001 → CONS-002 → CONS-003 → CONS-004 → (CONS-005 & CONS-006 in parallel) → CONS-011 → CONS-012 → CONS-015

## Key Risks

1. **Test Failures**: Consolidation may reveal additional issues
2. **Performance**: Modular approach might be slower
3. **API Changes**: Some internal APIs might need adjustment
4. **Hidden Dependencies**: Tests might depend on implementation details

## Success Metrics

1. **Test Pass Rate**: 100% (0 failures)
2. **Code Reduction**: Remove ~500+ lines of duplicate code
3. **Modular Coverage**: All major components in separate files
4. **AR+Robust Working**: Seamless "whiten then robustly weight" pipeline