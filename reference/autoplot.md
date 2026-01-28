# Autoplot method for Reg objects

Creates a ggplot visualization of an fMRI regressor object.

## Usage

``` r
# S3 method for class 'Reg'
autoplot(object, grid = NULL, precision = 0.1, method = "conv", ...)
```

## Arguments

- object:

  A `Reg` object (or one inheriting from it, like `regressor`).

- grid:

  Optional numeric vector specifying time points (seconds) for
  evaluation. If NULL, a default grid is generated based on the object's
  onsets and span.

- precision:

  Numeric precision for HRF evaluation if `grid` needs generation or if
  internal evaluation requires it (passed to `evaluate`).

- method:

  Evaluation method passed to `evaluate`.

- ...:

  Additional arguments (currently unused).

## Value

A ggplot object.
