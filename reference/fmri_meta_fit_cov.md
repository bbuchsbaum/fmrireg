# Fit Meta-Analysis and return packed covariance per voxel

Fit Meta-Analysis and return packed covariance per voxel

## Usage

``` r
fmri_meta_fit_cov(
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

List with base outputs and cov_tri (tsize x P) where tsize = K\*(K+1)/2
