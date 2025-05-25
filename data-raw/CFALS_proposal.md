# CF-ALS Integration Proposal for fmrireg

## Overview

**Goal:** Integrate the CF-ALS (Confound-Free Alternating Least Squares) rank-1 HRF estimation method into the `fmrireg` R package to provide users with a fast, accurate, and robust alternative for data-driven HRF estimation.

## Algorithm Core: CF-ALS Method

### Mathematical Foundation
CF-ALS estimates hemodynamic response functions (HRFs) and activation amplitudes simultaneously using a rank-1 decomposition model: **Y ≈ D(h)β^T**, where:
- **Y**: n × v BOLD data matrix (n timepoints, v voxels)
- **h**: d × 1 HRF coefficient vector (d basis functions)
- **β**: k × v amplitude matrix (k conditions, v voxels)
- **D(h)**: n × k design matrix formed by convolving stimulus onsets with HRF h

### Key Algorithmic Features
1. **Confound Projection**: QR-based orthogonal projection to remove nuisance regressors
2. **Efficient Precomputation**: X^T X and X^T Y matrices computed once
3. **SVD Initialization**: Regularized least squares followed by SVD for robust starting values
4. **Single Alternation**: Empirically, one β-h alternation post-initialization is sufficient
5. **Separate Regularization**: Independent penalties λ_β and λ_h for fine-grained control
6. **Vectorized Identifiability**: Efficient scaling and sign alignment across all voxels

### Algorithm Steps

**Input:**
- Y: n × v BOLD data matrix
- X.list: List of k design matrices (n × d each), one per condition
- Z: Optional n × q confound matrix
- λ_β, λ_h: Regularization parameters
- fullXtX: Boolean for cross-condition terms in h-update

**Steps:**
1. **Confound Projection**: Project Y and X.list orthogonally to Z using QR decomposition
2. **Precomputation**: 
   - X^T X[[c]] = X_c^T X_c for each condition c
   - X^T Y[[c]] = X_c^T Y for each condition c
   - If fullXtX: X^T X.full[[l,m]] = X_l^T X_m for all pairs l,m
3. **SVD Initialization**: 
   - Form X_big = [X_1, ..., X_k]
   - Solve (X_big^T X_big + λ_init I)γ = X_big^T Y
   - For each voxel: SVD of reshaped γ to initialize h and β
4. **CF-ALS Alternation** (typically 1 iteration):
   - β-update: Solve (G_mat + λ_β I)β_v = D(h)^T y_v
   - h-update: Solve (LHS + λ_h I)h_v = RHS
5. **Identifiability**: Normalize ||h||_∞ = 1 and align sign with canonical HRF

## Integration with fmrireg HRF System

### HRF Basis Compatibility
CF-ALS will work with **any HRF basis from the fmrireg library** (provided nbasis > 1), including:

**From `hrf.R` registry:**
- **bspline**: B-spline basis functions (flexible HRF modeling)
- **tent**: Tent/triangular basis functions  
- **fourier**: Fourier basis functions
- **daguerre**: Daguerre basis functions
- **spmg2/spmg3**: SPM canonical with derivatives
- **Custom generators**: Any user-defined HRF basis

### Technical Implementation

**Input Interface:**
```r
fmrireg_hrf_cfals(
    fmri_data_obj, 
    event_data_obj, 
    hrf_basis,        # fmrireg HRF basis object
    confound_obj = NULL,
    lam_beta = 10, 
    lam_h = 1,
    R_mat = NULL,     # Penalty matrix for h coefficients
    fullXtX = FALSE,
    max_alt = 1,
    ...
)
```

**HRF Basis Processing:**
1. **Dimensionality**: Extract d = nbasis(hrf_basis)
2. **Design Matrices**: Use fmrireg's existing convolution system to generate X.list
3. **Reconstruction Matrix Φ**: Extract from hrf_basis for converting coefficients h to HRF shapes
4. **Penalty Matrix R**: Default to diag(d), or use basis-specific smoothness penalties
5. **Identifiability**: Scale and align reconstructed HRF shapes (Φh) rather than raw coefficients

**Key Modifications for Generic Basis Support:**
- Replace FIR-specific logic with basis-agnostic coefficient handling
- Use Φ matrix for HRF shape reconstruction in identifiability constraints
- Support basis-specific penalty matrices R for h regularization
- Generate canonical reference HRF in same basis space for sign alignment

## Implementation Plan

### Core Development (MVP)

**CFALS-CORE-01**: CF-ALS Engine Implementation
- Adapt lean R algorithm for generic HRF basis support
- Implement numerical stability checks and error handling
- Support fullXtX toggle and configurable regularization

**CFALS-DESIGN-02**: HRF Basis Integration  
- Interface with fmrireg HRF system for design matrix generation
- Extract reconstruction matrices and penalty matrices from HRF basis objects
- Handle confound projection with QR decomposition

**CFALS-WRAP-03**: User Interface
- Main `fmrireg_hrf_cfals()` function accepting standard fmrireg objects
- Parameter validation and sensible defaults
- Integration with fmrireg data structures

**CFALS-OUTPUT-04**: Output Object & Methods
- `fmrireg_cfals_fit` S3 class with comprehensive results
- Methods: print, summary, plot, residuals
- Store both coefficients h and reconstructed HRF shapes

### Enhanced Output Object

```r
fmrireg_cfals_fit:
  h                 # d × v matrix of HRF coefficients  
  beta              # k × v matrix of condition amplitudes
  reconstructed_hrfs # p × v matrix of HRF shapes (Φh)
  residuals         # n × v matrix of model residuals
  hrf_basis_used    # Original HRF basis object
  lambda_used       # Regularization parameters applied
  design_info       # Dimensions and flags (d, k, n, fullXtX)
  gof_per_voxel     # Optional: R² and other fit statistics
```

### Validation & Documentation

**CFALS-TEST-05**: Comprehensive Testing
- Unit tests for algorithmic correctness across different HRF bases
- Integration tests with real fMRI datasets
- Comparison with canonical GLM approaches
- Edge case handling (low signal, different basis types)

**CFALS-DOC-06**: User Documentation
- Complete parameter documentation with usage guidance
- Examples demonstrating different HRF basis types
- Interpretation guidelines for outputs
- Performance and memory considerations

## Future Enhancements

### Advanced Features (Post-MVP)
- **Parameter Selection**: Cross-validation for λ_β and λ_h optimization
- **Performance**: Rcpp optimization for voxel loops
- **Scalability**: Chunked processing for memory-limited scenarios  
- **Diagnostics**: Advanced visualization and comparison utilities
- **Basis Extensions**: Support for custom penalty matrices and constraints

### Integration Benefits
- **Consistency**: Leverages existing fmrireg HRF infrastructure
- **Flexibility**: Works with any basis type (nbasis > 1) from the registry
- **Efficiency**: Single alternation provides speed advantage over iterative methods
- **Robustness**: Separate regularization and numerical stability features
- **Usability**: Familiar interface for fmrireg users with comprehensive output objects

This integration will provide fmrireg users with a powerful, flexible tool for data-driven HRF estimation that seamlessly works with the package's existing HRF basis system while offering significant computational and methodological advantages.
