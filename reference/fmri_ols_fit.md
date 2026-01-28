# OLS Fit with Optional Voxelwise Covariate

Wrapper for OLS t-tests that supports an optional voxelwise covariate.

## Usage

``` r
fmri_ols_fit(
  Y,
  X,
  voxelwise = NULL,
  center_voxelwise = TRUE,
  voxel_name = "voxel_cov"
)
```

## Arguments

- Y:

  Outcome matrix (S x P)

- X:

  Design matrix (S x K)

- voxelwise:

  Optional voxelwise covariate matrix (S x P)

- center_voxelwise:

  Logical; center voxelwise covariate per feature

- voxel_name:

  Name for voxelwise coefficient

## Value

List with OLS results
