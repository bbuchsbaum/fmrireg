# Extract Estimated MA Parameters from fmri_lm Fit

Retrieves the estimated moving-average parameters from a fitted fMRI
linear model that used ARMA or MA temporal-noise modeling.

## Usage

``` r
ma_parameters(object, ...)

# S3 method for class 'fmri_lm'
ma_parameters(object, scope = c("average", "per_run", "raw"), ...)
```

## Arguments

- object:

  An object of class `fmri_lm`.

- ...:

  Additional arguments (currently unused).

- scope:

  Character; `"average"` (default) returns coefficients averaged across
  runs, `"per_run"` returns the run-level estimates, and `"raw"` returns
  the stored structure without post-processing.

## Value

Depending on `scope`, either a numeric vector of averaged MA
coefficients, a list of per-run coefficient vectors, or the raw stored
structure. Returns `NULL` when no MA modeling was performed.
