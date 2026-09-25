# correlation_map.fmri_model

Heatmap of pairwise correlations between the columns of an
`fmri_model`'s combined event + baseline design matrix, grouped into
panels by model term (event terms, drift, nuisance).

By default (`within_run = TRUE`) each column is centred within each run
and the run-intercept columns are dropped before correlating. This shows
the collinearity that matters for estimation once run means are
modelled, and removes spurious correlations that zero-padding creates
between columns belonging to different runs. Cells for pairs of columns
that live in different runs are left blank.

## Usage

``` r
# S3 method for class 'fmri_model'
correlation_map(
  x,
  method = c("pearson", "spearman"),
  half_matrix = TRUE,
  absolute_limits = TRUE,
  within_run = TRUE,
  label_values = NULL,
  ...
)
```

## Arguments

- x:

  An `fmri_model`.

- method:

  Correlation method (e.g., "pearson", "spearman").

- half_matrix:

  Logical; if TRUE (default), display only the lower triangle.

- absolute_limits:

  Logical; if TRUE, the colour scale spans -1 to 1; otherwise it spans
  plus or minus the largest absolute correlation shown.

- within_run:

  Logical; centre columns within runs and drop run intercepts before
  correlating (default `TRUE`). `FALSE` correlates the raw columns.

- label_values:

  Logical or `NULL`; print correlations in the cells. `NULL` (default)
  labels every cell when there are at most 12 columns, and otherwise
  only cells with \|r\| \>= 0.3 that involve an event regressor.

- ...:

  Additional arguments passed to
  [`geom_tile`](https://ggplot2.tidyverse.org/reference/geom_tile.html).

## Value

A ggplot2 object.
