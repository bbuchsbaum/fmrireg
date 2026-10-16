# Extract Coefficients from an fmri_lm Fit

Returns estimated coefficients from a fitted `fmri_lm` model. Every form
of the result has one orientation: **voxels x terms**, with one row per
voxel (or latent component) and one column per coefficient or contrast,
and the term names on the column margin. This matches
[`stats()`](https://bbuchsbaum.github.io/fmrireg/reference/stats.md),
[`standard_error()`](https://bbuchsbaum.github.io/fmrireg/reference/standard_error.md)
and
[`p_values()`](https://bbuchsbaum.github.io/fmrireg/reference/p_values.md),
so results from the accessor family line up element-wise (for example
`coef(fit) / standard_error(fit)`).

## Usage

``` r
# S3 method for class 'fmri_lm'
coef(
  object,
  type = c("betas", "contrasts"),
  include_baseline = FALSE,
  recon = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted `fmri_lm` object.

- type:

  `"betas"` (default) for regression coefficients, or `"contrasts"` for
  the estimates of the simple contrasts defined in the model.

- include_baseline:

  Logical. For `type = "betas"`, also return the baseline (drift, block
  intercept and nuisance) coefficients, in design matrix column order.
  The default, `FALSE`, returns event coefficients only.

- recon:

  Logical; used by latent-space fits to reconstruct voxel-space
  coefficients. Ignored for ordinary fits.

- ...:

  Unused.

## Value

For `type = "betas"`, a numeric matrix with one row per voxel and one
column per coefficient, columns named after the design matrix (see
[`coef_names()`](https://bbuchsbaum.github.io/fmrireg/reference/coef_names.md)).
For `type = "contrasts"`, a tibble with one row per voxel and one column
per contrast.

## Orientation change in fmrireg 0.3.0

Before 0.3.0 the default call (`type = "betas"`,
`include_baseline = FALSE`) returned the transpose, terms x voxels,
while the other forms were voxels x terms. Code written for the old
default that indexed `coef(fit)[term, ]` must now use
`coef(fit)[, term]`.

## See also

[`coef_names()`](https://bbuchsbaum.github.io/fmrireg/reference/coef_names.md),
[`stats()`](https://bbuchsbaum.github.io/fmrireg/reference/stats.md)
(t-statistics),
[`standard_error()`](https://bbuchsbaum.github.io/fmrireg/reference/standard_error.md),
[`p_values()`](https://bbuchsbaum.github.io/fmrireg/reference/p_values.md),
[`coef_image()`](https://bbuchsbaum.github.io/fmrireg/reference/coef_image.md).

## Examples

``` r
X <- matrix(rnorm(50 * 3), 50, 3)
edata <- data.frame(
  condition = factor(c("A", "B", "A", "B")),
  onsets = c(1, 12, 25, 38),
  run = 1
)
dset <- matrix_frame(X, TR = 2, run_length = 50, event_table = edata)
fit <- fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset)
b <- coef(fit) # 3 voxels x 2 event coefficients
dim(b)
#> [1] 3 2
colnames(b)
#> [1] "condition_condition.A" "condition_condition.B"
dim(coef(fit, include_baseline = TRUE))
#> [1] 3 6
```
