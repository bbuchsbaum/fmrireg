# fMRI Linear Model Consolidation Plan

## Current State Analysis

We have a fragmented implementation with three overlapping codebases:

1. **Original `fmrilm.R`** (2076 lines)
   - Monolithic implementation with 22 functions
   - Contains the complete original functionality
   - Exports `fmri_lm` and other methods

2. **Partial refactor `fmrilm_new.R`** (452 lines)
   - Simplified implementation with 5 functions
   - Exports `fmri_lm` (overrides the original!)
   - Missing many functions that tests depend on

3. **Modular components** (created during refactoring)
   - `fmri_lm_config.R` - Configuration management
   - `fmri_lm_context.R` - GLM context objects
   - `fmri_lm_solver.R` - Core solver
   - `fmri_ar_modeling.R` - AR utilities
   - `fmri_robust_fitting.R` - Robust fitting engine
   - `fmri_lm_integrated_solver.R` - Integrated AR+Robust solver
   - Various other modular files

## Problems

1. **Function Override**: `fmrilm_new.R` loads after `fmrilm.R` and overrides the main `fmri_lm` function
2. **Missing Functionality**: The new implementation lacks many functions, causing 132 test failures
3. **Incomplete Integration**: The modular components aren't fully integrated with the main API
4. **Duplicated Code**: Same functions implemented in multiple places
5. **AR/Robust Integration**: The "whiten then robustly weight" pipeline isn't properly exposed

## Consolidation Goals

Following the principles in `LM_Refactor.md`:

1. **Single Source of Truth**: One implementation of each function
2. **Modular Architecture**: Clear separation of concerns
3. **Seamless AR+Robust**: "Whiten then robustly weight" as the primary pipeline
4. **Backward Compatibility**: Existing tests should pass
5. **Clean API**: User-facing functions in `fmrilm.R`, internals in modular files

## Consolidation Strategy

### Phase 1: Immediate Stabilization
1. Rename `fmrilm_new.R` to `fmrilm_new.R.bak` to stop function override
2. Ensure tests can run with original implementation
3. Document what functionality is missing from modular components

### Phase 2: Integration
1. Update `fmrilm.R` to use modular components where they exist
2. Implement missing pieces in modular files
3. Ensure AR+Robust pipeline works end-to-end

### Phase 3: Cleanup
1. Remove duplicated code
2. Move all internal functions to appropriate modular files
3. Keep only public API in `fmrilm.R`

### Phase 4: Validation
1. All tests pass
2. AR+Robust integration works seamlessly
3. Performance is maintained or improved

## Risk Mitigation

1. **Backup Current State**: Tag current commit before changes
2. **Incremental Changes**: Small, testable commits
3. **Test Coverage**: Run tests after each change
4. **Documentation**: Update as we go

## Success Criteria

1. Single, coherent implementation
2. All 1450+ tests passing
3. AR+Robust "whiten then robustly weight" working
4. Clean modular architecture per `LM_Refactor.md`
5. No duplicated code
6. Clear separation of public API and internals