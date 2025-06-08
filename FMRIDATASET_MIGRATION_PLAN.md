# fmridataset Migration Plan

## Overview
Migrate fmrireg to use the external fmridataset package instead of internal dataset code in `R/fmri_dataset.R`.

## Phase 1: Analysis (Current State)

### Files to Migrate
1. **R/fmri_dataset.R** - Contains all dataset classes and methods
   - Classes: `matrix_dataset`, `fmri_mem_dataset`, `latent_dataset`, `fmri_file_dataset`
   - Functions: `read_fmri_config`, `data_chunks`, chunking strategies
   - Methods: `get_data`, `get_data_matrix`, `get_mask`, `blocklens`, `as.matrix_dataset`

### Dependencies to Check
Files that import or use dataset functionality:
- `R/fmri_lm*.R` - Linear model fitting functions
- `R/fmri_model.R` - Model construction
- Tests in `tests/testthat/test_dataset.R`
- Any other files using `matrix_dataset()`, `fmri_mem_dataset()`, etc.

## Phase 2: Preparation

### 1. Add fmridataset Dependency
```r
# In DESCRIPTION:
Imports:
  fmridataset (>= 0.1.0)
```

### 2. Create Import File
Create `R/fmridataset-imports.R`:
```r
#' @import fmridataset
#' @importFrom fmridataset matrix_dataset fmri_mem_dataset latent_dataset
#' @importFrom fmridataset get_data get_data_matrix get_mask blocklens
#' @importFrom fmridataset data_chunks as.matrix_dataset
#' @importFrom fmridataset read_fmri_config sampling_frame
NULL

# Re-export commonly used functions
#' @export
fmridataset::matrix_dataset

#' @export
fmridataset::fmri_mem_dataset

#' @export
fmridataset::latent_dataset

#' @export
fmridataset::read_fmri_config
```

## Phase 3: Migration Steps

### 1. Namespace Updates
- Remove generics from `R/all_generic.R` that conflict:
  - `get_data`
  - `get_data_matrix`  
  - `get_mask`
  - `as.matrix_dataset`
  - `data_chunks`

### 2. Update Function Calls
Replace unqualified calls with namespace-qualified versions:
- `matrix_dataset()` → `fmridataset::matrix_dataset()`
- `get_data()` → `fmridataset::get_data()`
- etc.

### 3. Handle API Differences
Check for any API differences between internal and external versions:
- Parameter names
- Return values
- Class structures

### 4. Update Tests
- Modify tests to use fmridataset functions
- Add explicit namespace qualification where needed

## Phase 4: File Removal

After confirming everything works:
1. Remove `R/fmri_dataset.R`
2. Remove any dataset-specific methods from other files
3. Update NAMESPACE and documentation

## Phase 5: Validation

### Test Strategy
1. Run all existing tests
2. Check that all dataset creation works
3. Verify data access methods
4. Test chunking functionality
5. Ensure backward compatibility

### Key Areas to Test
- Model fitting with different dataset types
- Chunked processing
- Memory-mapped datasets
- Config file reading

## Potential Issues

### 1. Namespace Conflicts
- Both packages may export the same generics
- Solution: Comment out conflicting generics in fmrireg

### 2. API Differences  
- External package may have slightly different signatures
- Solution: Create wrapper functions if needed

### 3. Class Inheritance
- Check that class hierarchies are compatible
- Ensure S3 dispatch works correctly

### 4. Backend Support
- fmridataset has pluggable backends (nifti, h5, matrix)
- Ensure fmrireg code doesn't assume specific backend

## Implementation Order

1. **Add dependency** - Update DESCRIPTION
2. **Create imports** - Add fmridataset-imports.R
3. **Fix namespaces** - Update conflicting generics
4. **Update code** - Add namespace qualifications
5. **Test** - Run all tests
6. **Remove old code** - Delete fmri_dataset.R
7. **Document** - Update any documentation

## Success Criteria

- [ ] All tests pass
- [ ] No functions from fmri_dataset.R remain
- [ ] Dataset creation works for all types
- [ ] Chunked processing works
- [ ] No namespace conflicts
- [ ] Backward compatibility maintained