# Extract Standard Errors from a Model Fit

Extract standard errors of parameter estimates from a fitted model
object. This is part of a family of functions for extracting statistical
measures.

## Usage

``` r
standard_error(x, ...)

# S3 method for class 'fmri_latent_lm'
standard_error(x, type = c("estimates", "contrasts"), recon = FALSE, ...)

# S3 method for class 'fmri_lm'
standard_error(x, type = c("estimates", "contrasts"), ...)

# S3 method for class 'fmri_lm'
standard_error(x, type = c("estimates", "contrasts"), ...)
```

## Arguments

- x:

  The fitted model object

- ...:

  Additional arguments passed to methods

- type:

  The type of standard errors to extract: "estimates" or "contrasts"
  (default: "estimates")

- recon:

  Logical; whether to reconstruct the full matrix representation
  (default: FALSE)

## Value

A tibble or matrix containing standard errors of parameter estimates

## See also

Other statistical_measures:
[`p_values()`](https://bbuchsbaum.github.io/fmrireg/reference/p_values.md),
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

# Extract standard errors
se <- standard_error(fit)
```
