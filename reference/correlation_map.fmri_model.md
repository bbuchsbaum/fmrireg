# correlation_map.fmri_model

Generates a correlation heatmap of the columns in an `fmri_model`'s
combined event+baseline design matrix.

## Usage

``` r
# S3 method for class 'fmri_model'
correlation_map(
  x,
  method = c("pearson", "spearman"),
  half_matrix = FALSE,
  absolute_limits = TRUE,
  ...
)
```

## Arguments

- x:

  An `fmri_model`.

- method:

  Correlation method (e.g., "pearson", "spearman").

- half_matrix:

  Logical; if TRUE, display only the lower triangle of the matrix.

- absolute_limits:

  Logical; if TRUE, set color scale limits from -1 to 1.

- ...:

  Additional arguments passed to internal plotting functions.
