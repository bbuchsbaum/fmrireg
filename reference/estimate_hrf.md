# Estimate smooth condition-level hemodynamic response functions

`estimate_hrf()` estimates condition-level HRF curves directly from an
fMRI dataset. It builds an event-aligned spline basis, removes baseline
and fixed nuisance effects once, and solves every voxel together with
one penalized multiresponse linear system. This replaces the former
voxel-by-voxel GAM implementation, whose predictor did not represent
post-stimulus time.

## Usage

``` r
estimate_hrf(
  form,
  fixed = NULL,
  block,
  dataset,
  bs = NULL,
  rsam = seq(0, 20, by = 1),
  basemod = NULL,
  k = 8L,
  fx = NULL,
  progress = FALSE,
  basis = c("bspline", "tent"),
  lambda = "gcv",
  lambda_grid = c(0, 10^seq(-4, 4, length.out = 25L)),
  ci_level = 0.95
)
```

## Arguments

- form:

  A two-sided event-model formula containing one or more
  [`hrf()`](https://bbuchsbaum.github.io/fmridesign/reference/hrf.html)
  terms.
  [`trialwise()`](https://bbuchsbaum.github.io/fmridesign/reference/trialwise.html)
  terms are rejected because this function estimates condition-level
  curves.

- fixed:

  Optional event-model formula whose design is treated as nuisance.

- block:

  Formula identifying acquisition runs or blocks.

- dataset:

  An `fmri_dataset`.

- bs:

  Deprecated legacy basis selector. Values from the former GAM API are
  mapped to `basis = "bspline"`.

- rsam:

  Strictly increasing, finite post-stimulus times beginning at zero.

- basemod:

  Optional baseline model. The default is a constant baseline.

- k:

  Number of free HRF basis functions per curve. At least four for the
  default cubic B-spline basis.

- fx:

  Deprecated legacy selector. `fx = TRUE` maps to `lambda = 0` when
  `lambda` is omitted.

- progress:

  Show progress while scanning the GCV grid.

- basis:

  Either `"bspline"` (cubic) or `"tent"` (piecewise linear).

- lambda:

  Non-negative smoothing strength or `"gcv"` for shared automatic
  selection.

- lambda_grid:

  Candidate smoothing strengths used for GCV.

- ci_level:

  Confidence level in `(0, 1)`, or `NULL` to omit intervals.

## Value

An object of class `fmri_hrf_estimate`. Its `estimate` and `std.error`
arrays have dimensions time by curve by voxel. The object also contains
labeled curve metadata, basis coefficients, smoothing diagnostics, and
the designs used for estimation. Use
[`tidy()`](https://bbuchsbaum.github.io/fmrireg/reference/tidy.md),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`coef()`](https://rdrr.io/r/stats/coef.html), or
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) for common
downstream representations.

## Details

A single smoothing parameter is shared across voxels. With
`lambda = "gcv"`, it is selected by scale-normalized generalized
cross-validation so high- variance voxels do not dominate the choice.
Curves are constrained to zero at the beginning and end of `rsam`.
Standard errors and confidence intervals use the fitted
penalized-linear-model covariance under independent, homoscedastic
time-point errors; they are not autocorrelation-robust.

## See also

[`fitted_hrf()`](https://bbuchsbaum.github.io/fmrireg/reference/fitted_hrf.md),
[`event_model()`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html),
[`baseline_model()`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html)

## Examples

``` r
set.seed(18)
n <- 80L
events <- data.frame(
  onset = seq(6, 62, by = 8),
  condition = factor(rep(c("A", "B"), 4)),
  run = 1L
)
dataset <- fmridataset::matrix_dataset(
  matrix(rnorm(n * 2), nrow = n),
  TR = 1,
  run_length = n,
  event_table = events
)
fit <- estimate_hrf(
  onset ~ hrf(condition),
  block = ~run,
  dataset = dataset,
  rsam = 0:12,
  k = 6,
  lambda = 1
)
tidy(fit, voxel = 1)
#> # A tibble: 26 × 9
#>     time curve       term      condition   voxel estimate std.error  lower upper
#>    <dbl> <chr>       <chr>     <chr>       <chr>    <dbl>     <dbl>  <dbl> <dbl>
#>  1     0 condition.A condition condition.A voxe…   0          0      0     0    
#>  2     1 condition.A condition condition.A voxe…   0.407      0.454 -0.499 1.31 
#>  3     2 condition.A condition condition.A voxe…   0.487      0.445 -0.401 1.37 
#>  4     3 condition.A condition condition.A voxe…   0.409      0.371 -0.331 1.15 
#>  5     4 condition.A condition condition.A voxe…   0.265      0.351 -0.434 0.964
#>  6     5 condition.A condition condition.A voxe…   0.0847     0.345 -0.603 0.772
#>  7     6 condition.A condition condition.A voxe…  -0.101      0.347 -0.793 0.590
#>  8     7 condition.A condition condition.A voxe…  -0.258      0.369 -0.993 0.477
#>  9     8 condition.A condition condition.A voxe…  -0.353      0.415 -1.18  0.475
#> 10     9 condition.A condition condition.A voxe…  -0.378      0.506 -1.39  0.630
#> # ℹ 16 more rows
```
