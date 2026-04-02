# Extract Image/Volume for Coefficient

Extract Image/Volume for Coefficient

## Usage

``` r
# S3 method for class 'fmri_lm'
coef_image(
  object,
  coef = 1,
  statistic = c("estimate", "se", "tstat", "prob"),
  type = c("estimates", "contrasts", "F"),
  ...
)

coef_image(object, coef = 1, statistic = c("estimate", "se", "z", "p"), ...)
```

## Arguments

- object:

  An fmri_meta object

- coef:

  Coefficient name or index

- statistic:

  Type of statistic to extract ("estimate", "se", "z", "p")

- type:

  For `fmri_lm` objects: which coefficient set to index into:
  `"estimates"` (default), `"contrasts"`, or `"F"`.

- ...:

  Additional arguments (currently unused).

## Value

NeuroVol object or matrix

## Examples

``` r
# Create a small example
X <- matrix(rnorm(50 * 4), 50, 4)
edata <- data.frame(
  condition = factor(c("A", "B", "A", "B")),
  onsets = c(1, 12, 25, 38),
  run = c(1, 1, 1, 1)
)
dset <- fmridataset::matrix_dataset(X, TR = 2, run_length = 50,
                                    event_table = edata)
fit <- fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset)
# Get coefficient estimates as a numeric vector
coef_image(fit, coef = 1)
#> [1]  0.08603959 -0.05707743 -0.59749063 -0.11075571
toy_meta <- structure(
  list(
    coefficients = matrix(c(0.3, 0.1), nrow = 1,
      dimnames = list(NULL, c("A", "B"))),
    se = matrix(c(0.05, 0.06), nrow = 1)
  ),
  class = "fmri_meta"
)
coef_image(toy_meta, coef = "A")
#>   A 
#> 0.3 
```
