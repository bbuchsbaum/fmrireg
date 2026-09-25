# Heatmap visualization of the combined fmri_model design matrix

Produces a single heatmap of all columns in the design matrix of an
`fmri_model`, which merges the event_model and baseline_model
regressors. Rows are scans (time runs downward); columns are regressors,
kept in model order and grouped into panels by term (event terms, drift,
intercept, nuisance). Run boundaries are marked with horizontal rules,
and the subtitle reports the rank and condition number of the design.

By default each column is divided by its maximum absolute value, so
regressors on very different scales (an intercept, a motion trace, a
convolved event train) are all legible on one colour scale. Set
`scale_columns = FALSE` to show raw values. Run-intercept columns are
drawn in grey, since they carry no information beyond run membership.

The condition number is computed after scaling every column to unit
length, so it reflects collinearity rather than the units of the
columns.

## Usage

``` r
# S3 method for class 'fmri_model'
design_map(
  x,
  block_separators = TRUE,
  rotate_x_text = TRUE,
  fill_midpoint = NULL,
  fill_limits = NULL,
  scale_columns = TRUE,
  ...
)
```

## Arguments

- x:

  An `fmri_model` object.

- block_separators:

  Logical; if `TRUE`, draw rules between runs.

- rotate_x_text:

  Logical; if `TRUE`, rotate x-axis labels to vertical.

- fill_midpoint:

  Numeric or `NULL`; centre of the diverging colour scale. Defaults to
  0.

- fill_limits:

  Numeric vector of length 2 or `NULL`; passed to the fill scale
  `limits=` argument. This can clip or expand the color range.

- scale_columns:

  Logical; if `TRUE` (default), scale each column to a maximum absolute
  value of 1.

- ...:

  Additional arguments passed to
  [`geom_raster`](https://ggplot2.tidyverse.org/reference/geom_tile.html).

## Value

A ggplot2 plot object.
