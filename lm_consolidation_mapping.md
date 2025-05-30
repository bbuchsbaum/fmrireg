# LM Consolidation Function Mapping

## Current State Analysis

### Functions in fmrilm.R (Main Implementation)
1. **Public API Functions:**
   - `fmri_lm` - Main entry point for GLM fitting
   - `chunkwise_lm.fmri_dataset_old` - Old chunkwise implementation
   - `runwise_lm` - Run-by-run fitting
   - `multiresponse_lm` - Multi-response fitting
   - `create_fmri_model` - Model creation

2. **S3 Methods:**
   - `coef.fmri_lm`
   - `fitted_hrf.fmri_lm`
   - `print.fmri_lm`
   - `standard_error.fmri_lm`
   - `stats.fmri_lm`

3. **Internal Functions:**
   - `fmri_lm_fit` - Core fitting logic
   - `fit_lm_contrasts` - Contrast computation
   - `.fast_lm_matrix` - Fast matrix operations
   - `.fast_preproject` - Projection matrix computation
   - `fast_rlm_run` - Fast robust fitting per run
   - `pull_stat` / `pull_stat_revised` - Extract statistics
   - `reshape_coef` - Reshape coefficients
   - `unpack_chunkwise` - Unpack chunk results

4. **Helper Functions:**
   - `get_formula.fmri_model`
   - `term_matrices.fmri_model`
   - `is.formula`

### Functions in Modular Files

1. **fmri_lm_config.R / fmri_lm_config_wrapper.R:**
   - `fmri_lm_config` - Configuration object creation
   - `fmri_lm_control` - Control parameters

2. **fmri_lm_context.R:**
   - `glm_context` - Context object for data passing

3. **fmri_lm_solver.R:**
   - `solve_glm_core` - Core GLM solver

4. **fmri_lm_integrated_solver.R:**
   - `solve_integrated_glm` - Unified solver with AR+Robust
   - `solve_ols_pipeline`, `solve_ar_pipeline`, `solve_robust_pipeline`, `solve_ar_robust_pipeline`

5. **fmri_lm_ar_integration.R:**
   - `whiten_glm_context` - AR whitening
   - `iterative_ar_solve` - Iterative AR estimation
   - `compute_ar_effective_df` - AR effective df

6. **fmri_robust_fitting.R:**
   - `iterative_robust_solve` - IRLS implementation
   - `compute_robust_weights` - Weight computation

7. **fmri_lm_chunkwise.R:**
   - `chunkwise_lm.fmri_dataset` - New chunkwise implementation
   - `chunkwise_lm_fast`, `chunkwise_lm_slow`

8. **fmri_lm_runwise.R:**
   - `runwise_lm` - Duplicated from fmrilm.R
   - `runwise_lm_fast`, `runwise_lm_slow`, `runwise_lm_voxelwise`
   - `pool_runwise_results`

9. **fmri_lm_methods.R:**
   - All S3 methods (duplicated from fmrilm.R)

10. **fmri_lm_internal.R:**
    - `.fast_lm_matrix`, `.fast_preproject` (duplicated)
    - Voxelwise contrast functions
    - Meta analysis functions

## Gap Analysis

### Missing from Modular Implementation:
1. **Main Entry Point:**
   - `fmri_lm` function not in modular files - still in fmrilm.R
   - `fmri_lm_fit` core logic not fully migrated

2. **Model Creation:**
   - `create_fmri_model` not in modular files

3. **Multi-response:**
   - `multiresponse_lm` not in modular files

4. **Fast Robust:**
   - `fast_rlm_run` not in fmri_robust_fitting.R

5. **Unpacking:**
   - `unpack_chunkwise` not in modular files

### Duplicated Across Files:
1. **S3 Methods** - in both fmrilm.R and fmri_lm_methods.R
2. **Internal functions** - `.fast_lm_matrix`, `.fast_preproject` in both
3. **runwise_lm** - in both fmrilm.R and fmri_lm_runwise.R

### Integration Issues:
1. **fmri_lm_fit in fmrilm.R doesn't call modular solvers**
2. **AR parameter mapping not implemented**
3. **Fast path (`use_fast_path=TRUE`) bypasses modular code**

## Target Function Mapping

### fmrilm.R (Public API Only)
- `fmri_lm` - Calls modular implementation
- `create_fmri_model` - Keep here or move to fmri_model.R

### fmri_lm_core.R (New File)
- `fmri_lm_fit` - Refactored to use modular components
- `multiresponse_lm` - If keeping this functionality

### fmri_lm_internal.R (Already Exists)
- `.fast_lm_matrix` (already there)
- `.fast_preproject` (already there)
- `unpack_chunkwise` (move here)

### fmri_robust_fitting.R
- `fast_rlm_run` (move here)

### fmri_lm_methods.R (Already Exists)
- All S3 methods (already there)

### Delete/Remove
- Duplicate functions in fmrilm.R after migration
- Old implementations once new ones are working

## Critical Path Forward

1. **Add AR parameter mapping to fmri_lm** (CONS-003)
2. **Connect fmri_lm_fit to modular solvers** (CONS-004)
3. **Update runwise_lm to use modular components** (CONS-005)
4. **Update chunkwise_lm to use modular components** (CONS-006)
5. **Migrate remaining functions** (CONS-008)
6. **Remove duplicates** (CONS-010)

## Current Test Failures (12 total)

### Categories:
1. **AFNI HRF Aliases (2 failures)**
   - Missing `name` attribute on HRF objects
   
2. **AR Parameter Estimation (7 failures)**
   - AR parameters not being estimated correctly
   - Config structure issues (ar$struct vs ar_options)

3. **AR + Robust Combined (3 failures)**  
   - Config structure mismatch preventing AR+Robust combination
   - Need to ensure proper parameter routing