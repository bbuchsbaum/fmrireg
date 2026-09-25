# Autoplot method for Reg objects

Plots the predicted BOLD time course of an fMRI regressor: the event
train convolved with its HRF. Event onsets are marked as ticks beneath
the curve, so the hemodynamic lag is visible. When responses to nearby
events overlap, each event's own response is drawn faintly under the
summed curve. Multi-basis HRFs are drawn as stacked panels that share
the time axis, each with its own vertical scale.

## Usage

``` r
# S3 method for class 'Reg'
autoplot(
  object,
  grid = NULL,
  precision = 0.1,
  method = "conv",
  title = NULL,
  ...
)
```

## Arguments

- object:

  A `Reg` object (or one inheriting from it, like `regressor`).

- grid:

  Optional numeric vector specifying time points (seconds) for
  evaluation. If NULL, a default grid spans the first onset to the end
  of the last event's response (including its duration).

- precision:

  Numeric precision for HRF evaluation if `grid` needs generation or if
  internal evaluation requires it (passed to `evaluate`).

- method:

  Evaluation method passed to `evaluate`.

- title:

  Optional plot title, e.g. the condition this regressor models.

- ...:

  Additional arguments (currently unused).

## Value

A ggplot object.
