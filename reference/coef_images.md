# Extract Images/Volumes for Multiple Coefficients

A plural companion to
[`coef_image`](https://bbuchsbaum.github.io/fmrireg/reference/coef_image.md).
Returns a named list of volumes (or numeric vectors for non-spatial
datasets), one per coefficient, so callers can iterate over or write out
every coefficient without reimplementing the loop over
[`coef_names`](https://bbuchsbaum.github.io/fmrireg/reference/coef_names.md).

## Usage

``` r
# S3 method for class 'fmri_lm'
coef_images(
  object,
  statistic = c("estimate", "se", "tstat", "prob"),
  type = c("estimates", "contrasts", "F"),
  coefs = NULL,
  ...
)

coef_images(object, ...)
```

## Arguments

- object:

  A fitted model object.

- statistic:

  For `fmri_lm` objects: one of `"estimate"`, `"se"`, `"tstat"`, or
  `"prob"`.

- type:

  For `fmri_lm` objects: which coefficient set to extract: `"estimates"`
  (default), `"contrasts"`, or `"F"`.

- coefs:

  Optional character vector of coefficient names (a subset of
  `coef_names(object, type = type)`) to extract. Defaults to all
  coefficients of the requested `type`.

- ...:

  Additional arguments passed to methods (and on to
  [`coef_image`](https://bbuchsbaum.github.io/fmrireg/reference/coef_image.md)).

## Value

A named list of `NeuroVol` objects (or numeric vectors when the dataset
is non-spatial), named by coefficient.

## See also

[`coef_image`](https://bbuchsbaum.github.io/fmrireg/reference/coef_image.md),
[`coef_names`](https://bbuchsbaum.github.io/fmrireg/reference/coef_names.md)

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
# Named list of one volume per event regressor
imgs <- coef_images(fit, statistic = "estimate", type = "estimates")
names(imgs)
#> [1] "condition_condition.A" "condition_condition.B"
```
