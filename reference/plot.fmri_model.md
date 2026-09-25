# Plot the regressors of an fmri_model over time

Draws the design-matrix columns as time series in stacked small
multiples that share one time axis. Each event regressor gets its own
row, coloured by condition (designs with more than 12 event columns get
one row per event term instead). Baseline regressors are drawn in grey,
one row per baseline term (drift, nuisance). Run-specific columns are
drawn only within their own run, and no line is joined across a run
boundary. Run intercepts are omitted, since a constant carries no
temporal information.

## Usage

``` r
# S3 method for class 'fmri_model'
plot(x, baseline = TRUE, ...)
```

## Arguments

- x:

  An `fmri_model`.

- baseline:

  Logical; include baseline rows (default `TRUE`).

- ...:

  Unused.

## Value

A ggplot2 object.
