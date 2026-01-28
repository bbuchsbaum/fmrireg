# Fit Meta-Analysis Models

Low-level function that calls the C++ meta-analysis implementation. This
is typically called internally by higher-level functions like
fmri_meta().

## Usage

``` r
fmri_meta_fit(
  Y,
  V,
  X,
  method = c("pm", "dl", "fe", "reml"),
  robust = c("none", "huber"),
  huber_c = 1.345,
  robust_iter = 2,
  n_threads = getOption("fmrireg.num_threads", 0)
)
```

## Arguments

- Y:

  Numeric matrix of effect sizes (subjects x features)

- V:

  Numeric matrix of variances (subjects x features)

- X:

  Numeric matrix; design matrix (subjects x predictors), including
  intercept

- method:

  Character scalar; meta-analysis method: "pm" (Paule-Mandel), "dl"
  (DerSimonian-Laird), "fe" (fixed-effects), or "reml" (REML, uses PM
  solver)

- robust:

  Character scalar; robust estimation method: "none" or "huber"

- huber_c:

  Numeric scalar; tuning constant for Huber M-estimator (default:
  1.345). Smaller values provide more robust estimates but may reduce
  efficiency.

- robust_iter:

  Integer scalar; number of IRLS iterations for robust estimation
  (default: 2)

- n_threads:

  Integer scalar; number of OpenMP threads (0 = use all available)

## Value

List with components:

- beta:

  Numeric matrix of coefficients (predictors x features)

- se:

  Numeric matrix of standard errors (predictors x features)

- z:

  Numeric matrix of z-scores (predictors x features)

- tau2:

  Numeric vector of between-study variance estimates

- Q_fe:

  Numeric vector of Q statistics from fixed-effects model

- I2_fe:

  Numeric vector of I² statistics from fixed-effects model

- df:

  Numeric vector of degrees of freedom

- ok:

  Logical vector indicating successful fits

## See also

[`fmri_meta`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md)
