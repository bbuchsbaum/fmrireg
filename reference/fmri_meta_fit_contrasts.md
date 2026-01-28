# Fit Meta-Analysis Models with Exact Contrasts

Low-level function that calls the C++ meta-analysis implementation and
returns exact contrast statistics c' (X' W X)^(-1) c for provided
contrasts.

## Usage

``` r
fmri_meta_fit_contrasts(
  Y,
  V,
  X,
  Cmat,
  method = c("pm", "dl", "fe", "reml"),
  robust = c("none", "huber"),
  huber_c = 1.345,
  robust_iter = 2,
  n_threads = getOption("fmrireg.num_threads", 0)
)
```

## Arguments

- Y, V, X:

  See fmri_meta_fit

- Cmat:

  Numeric matrix K x J (columns are contrasts over predictors)

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

List with base meta outputs plus c_beta, c_se, c_z
