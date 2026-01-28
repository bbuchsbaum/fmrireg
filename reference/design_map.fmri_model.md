# Heatmap visualization of the combined fmri_model design matrix

Produces a single heatmap of *all* columns in the design matrix from an
`fmri_model` object, which merges both the event_model and
baseline_model regressors. Rows are scans; columns are regressors.
Optionally draws horizontal lines between blocks (runs), and rotates
x-axis labels diagonally for readability.

## Usage

``` r
# S3 method for class 'fmri_model'
design_map(
  x,
  block_separators = TRUE,
  rotate_x_text = TRUE,
  fill_midpoint = NULL,
  fill_limits = NULL,
  ...
)
```

## Arguments

- x:

  An `fmri_model` object.

- block_separators:

  Logical; if `TRUE`, draw white horizontal lines between blocks.

- rotate_x_text:

  Logical; if `TRUE`, rotate x-axis labels by 45 degrees.

- fill_midpoint:

  Numeric or `NULL`; if not `NULL`, passed to
  [`scale_fill_gradient2`](https://ggplot2.tidyverse.org/reference/scale_gradient.html)
  to center the color scale (e.g. `fill_midpoint=0`).

- fill_limits:

  Numeric vector of length 2 or `NULL`; passed to the fill scale
  `limits=` argument. This can clip or expand the color range.

- ...:

  Additional arguments passed to
  [`geom_tile`](https://ggplot2.tidyverse.org/reference/geom_tile.html).

## Value

A ggplot2 plot object.
