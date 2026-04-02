# Get Available Coefficient Names

Return the names of available coefficients from a fitted model object.
This helps users discover which coefficient names can be passed to
[`coef_image`](https://bbuchsbaum.github.io/fmrireg/reference/coef_image.md)
or other extraction functions.

## Usage

``` r
coef_names(x, ...)

# S3 method for class 'fmri_lm'
coef_names(x, type = c("estimates", "contrasts", "F", "all"), ...)
```

## Arguments

- x:

  A fitted model object

- ...:

  Additional arguments passed to methods

- type:

  Which set of names to return: `"estimates"` (default) for event
  regressor names, `"contrasts"` for simple contrast names, `"F"` for
  F-contrast names, or `"all"` for a named list of all three.

## Value

A character vector of coefficient names

## See also

[`coef_image`](https://bbuchsbaum.github.io/fmrireg/reference/coef_image.md),
[`coef`](https://rdrr.io/r/stats/coef.html)

Other statistical_measures:
[`p_values()`](https://bbuchsbaum.github.io/fmrireg/reference/p_values.md),
[`standard_error()`](https://bbuchsbaum.github.io/fmrireg/reference/standard_error.md),
[`stats()`](https://bbuchsbaum.github.io/fmrireg/reference/stats.md)

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
coef_names(fit)
#> [1] "condition_condition.A" "condition_condition.B"
```
