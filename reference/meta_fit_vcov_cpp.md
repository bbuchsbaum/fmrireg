# Meta-regression with ONE voxelwise covariate

Meta-regression with ONE voxelwise covariate

## Usage

``` r
meta_fit_vcov_cpp(
  Y,
  V,
  X,
  C,
  method,
  robust,
  huber_c = 1.345,
  robust_iter = 2L,
  n_threads = 0L
)
```

## Arguments

- Y:

  S x P matrix of effect sizes

- V:

  S x P matrix of variances

- X:

  S x K design matrix

- C:

  S x P matrix of voxelwise covariates

- method:

  Meta-analysis method

- robust:

  Robust estimation method

- huber_c:

  Huber tuning constant

- robust_iter:

  Number of IRLS iterations

- n_threads:

  Number of OpenMP threads

## Value

List with results
