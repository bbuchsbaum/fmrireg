# Visualize the entire design matrix as a heatmap

Generate a heatmap visualization of a design matrix, showing regressor
values over time. This is useful for inspecting the temporal structure
of fMRI design matrices.

## Usage

``` r
# S3 method for class 'convolved_term'
conditions(x, ...)
```

## Arguments

- x:

  The model object (event_model, baseline_model, or fmri_model)

- ...:

  Additional arguments passed to methods. Common arguments include:

  rescale_cols

  :   Logical; if TRUE, columns are rescaled to (-1,1)

  block_separators

  :   Logical; if TRUE, draw white lines between blocks

  rotate_x_text

  :   Logical; if TRUE, rotate x-axis labels by 45 degrees

## Value

A ggplot2 object containing the design matrix heatmap
