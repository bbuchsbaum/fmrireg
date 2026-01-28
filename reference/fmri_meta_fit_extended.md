# Extended Meta-Analysis Fit with Voxelwise Covariate

Wrapper for meta-analysis that supports an optional voxelwise covariate.
This extends the basic fmri_meta_fit to handle per-voxel covariates.

## Usage

``` r
fmri_meta_fit_extended(
  Y,
  V,
  X,
  method = c("pm", "dl", "fe", "reml"),
  robust = c("none", "huber"),
  huber_c = 1.345,
  robust_iter = 2,
  voxelwise = NULL,
  center_voxelwise = TRUE,
  voxel_name = "voxel_cov",
  n_threads = getOption("fmrireg.num_threads", 0)
)
```

## Arguments

- Y:

  Matrix of effect sizes (S x P)

- V:

  Matrix of variances (S x P)

- X:

  Design matrix (S x K)

- method:

  Meta-analysis method

- robust:

  Robust estimation method

- huber_c:

  Huber tuning constant

- robust_iter:

  Number of IRLS iterations

- voxelwise:

  Optional voxelwise covariate matrix (S x P)

- center_voxelwise:

  Logical; center voxelwise covariate per feature

- voxel_name:

  Name for voxelwise coefficient

- n_threads:

  Number of threads

## Value

List with meta-analysis results
