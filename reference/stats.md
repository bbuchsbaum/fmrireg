# Extract Test Statistics from a Model Fit

Extract test statistics (e.g., t-statistics, F-statistics) from a fitted
model object. This is part of a family of functions for extracting
statistical measures.

## Usage

``` r
stats(x, ...)

# S3 method for class 'fmri_lm'
stats(x, type = c("estimates", "contrasts", "F"), ...)

# S3 method for class 'fmri_lm'
stats(x, type = c("estimates", "contrasts", "F"), ...)
```

## Arguments

- x:

  The fitted model object

- ...:

  Additional arguments passed to methods

- type:

  The type of statistics to extract: "estimates", "contrasts", or "F"
  (default: "estimates")

## Value

A tibble or matrix containing test statistics

## See also

Other statistical_measures:
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
dset <- fmridataset::matrix_dataset(
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
