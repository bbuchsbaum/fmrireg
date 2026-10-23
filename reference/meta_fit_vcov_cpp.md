# Meta-regression with one voxelwise covariate

Fits the meta-regression of
[`fmri_meta_fit()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta_fit.md)
with the design augmented, per feature, by that feature's column of `C`
as a last predictor.

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

  Meta-analysis method ("fe", "dl", "pm"; "reml" runs PM)

- robust:

  Robust estimation method ("none" or "huber")

- huber_c:

  Huber tuning constant

- robust_iter:

  Number of IRLS iterations

- n_threads:

  Number of OpenMP threads

## Value

List with beta, se, z ((K + 1) x P), tau2, Q_fe, I2_fe, df, ok, method
and robust
