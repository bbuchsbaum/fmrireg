# Extract Image/Volume for Coefficient

Extract one coefficient's statistic map from a fitted model, as a
`NeuroVol` when the underlying data are volumetric and as a numeric
vector otherwise.

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

  A fitted model: an `fmri_lm` object (from
  [`fmri_lm`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md))
  or an `fmri_meta` object (from
  [`fmri_meta`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md)).
  Group t-test fits have their own method, documented in
  [`coef_image.fmri_ttest_fit`](https://bbuchsbaum.github.io/fmrireg/reference/coef_image.fmri_ttest_fit.md).

- coef:

  Coefficient name or index.

- statistic:

  Statistic to extract. The vocabulary depends on the class of `object`:

  `fmri_lm`

  :   One of `"estimate"` (default), `"se"`, `"tstat"`, or `"prob"`.

  `fmri_meta`

  :   One of `"estimate"` (default), `"se"`, `"z"` (estimate / se), or
      `"p"` (two-sided normal p-value).

- type:

  For `fmri_lm` objects only: which coefficient set to index into:
  `"estimates"` (default), `"contrasts"`, or `"F"`.

- ...:

  Additional arguments passed to methods (currently unused).

## Value

A `NeuroVol` for volumetric data, otherwise a numeric vector with one
value per voxel or feature.

## Examples

``` r
# Create a small example
X <- matrix(rnorm(50 * 4), 50, 4)
edata <- data.frame(
  condition = factor(c("A", "B", "A", "B")),
  onsets = c(1, 12, 25, 38),
  run = c(1, 1, 1, 1)
)
dset <- matrix_frame(X, TR = 2, run_length = 50,
                                    event_table = edata)
fit <- fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset)
# Get coefficient estimates as a numeric vector
coef_image(fit, coef = 1)
#> [1] -2.674999 -5.722917  3.523104  7.586770
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
