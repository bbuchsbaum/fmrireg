# OLS with ONE voxelwise covariate

OLS with ONE voxelwise covariate

## Usage

``` r
ols_t_vcov_cpp(Y, X, C)
```

## Arguments

- Y:

  S x P matrix of outcomes

- X:

  S x K design matrix

- C:

  S x P matrix of voxelwise covariates

## Value

List with beta ((K+1) x P), se, t, df
