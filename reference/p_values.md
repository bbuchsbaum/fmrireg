# Extract P-values from a Model Fit

Extract p-values associated with parameter estimates or test statistics
from a fitted model object. This is part of a family of functions for
extracting statistical measures.

## Usage

``` r
p_values(x, ...)

# S3 method for class 'fmri_lm'
p_values(x, type = c("estimates", "contrasts"), ...)
```

## Arguments

- x:

  The fitted model object

- ...:

  Additional arguments passed to methods

- type:

  Character string specifying the type of p-values to extract. Options
  typically include "estimates" for parameter estimates and "contrasts"
  for contrast tests. Defaults to "estimates" in most methods.

## Value

A tibble or matrix containing p-values

## See also

Other statistical_measures:
[`standard_error()`](https://bbuchsbaum.github.io/fmrireg/reference/standard_error.md),
[`stats()`](https://bbuchsbaum.github.io/fmrireg/reference/stats.md)

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

# Extract p-values
pvals <- p_values(fit)
```
