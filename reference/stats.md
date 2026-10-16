# Extract Test Statistics from a Model Fit

Extract test statistics (e.g., t-statistics, F-statistics) from a fitted
model object. This is part of a family of functions for extracting
statistical measures.

## Usage

``` r
stats(x, ...)

# S3 method for class 'fmri_lm'
stats(x, type = c("betas", "contrasts", "F"), ...)
```

## Arguments

- x:

  The fitted model object

- ...:

  Additional arguments passed to methods

- type:

  For `fmri_lm` fits, the parameter family whose test statistics are
  returned: `"betas"` (default) gives the **t-statistics** of the event
  coefficients, `"contrasts"` the t-statistics of the simple contrasts,
  and `"F"` the F-statistics of the F-contrasts. `stats()` never returns
  the estimates themselves; use
  [`coef()`](https://rdrr.io/r/stats/coef.html) for those. `"estimates"`
  is a deprecated synonym of `"betas"` (it reads as "the estimates" but
  returns t-statistics) and warns; set
  `options(fmrireg.suppress_deprecation = TRUE)` to silence it.

## Value

A tibble or matrix containing test statistics. For `fmri_lm` fits, a
tibble with one row per voxel and one column per term, in the same
orientation as [`coef()`](https://rdrr.io/r/stats/coef.html).

## See also

[`coef()`](https://bbuchsbaum.github.io/fmrireg/reference/coef.fmri_lm.md)
for the estimates.

Other statistical_measures:
[`coef_names()`](https://bbuchsbaum.github.io/fmrireg/reference/coef_names.md),
[`p_values()`](https://bbuchsbaum.github.io/fmrireg/reference/p_values.md),
[`standard_error()`](https://bbuchsbaum.github.io/fmrireg/reference/standard_error.md)

## Examples

``` r
# Create example data
event_data <- data.frame(
  condition = factor(c("A", "B", "A", "B")),
  onsets = c(1, 10, 20, 30),
  run = c(1, 1, 1, 1)
)

# Create sampling frame and dataset
sframe <- sampling_frame(blocklens = 50, TR = 2)
dset <- matrix_frame(
  matrix(rnorm(50 * 2), 50, 2),
  TR = 2,
  run_length = 50,
  event_table = event_data
)

# Fit model
fit <- fmri_lm(
  onsets ~ hrf(condition),
  block = ~run,
  dataset = dset
)

# Extract test statistics
tstats <- stats(fit)
```
