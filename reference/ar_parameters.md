# Extract Estimated AR Parameters from fmri_lm Fit

Retrieves the estimated autoregressive parameters from a fitted fMRI
linear model that used AR error modeling.

## Usage

``` r
ar_parameters(object, ...)

# S3 method for class 'fmri_lm'
ar_parameters(object, scope = c("average", "per_run", "raw"), ...)
```

## Arguments

- object:

  An object of class `fmri_lm`

- ...:

  Additional arguments (currently unused)

- scope:

  Character; `"average"` (default) returns the pooled average AR
  coefficients, `"per_run"` returns a list of the run-level estimates,
  and `"raw"` returns the stored structure without post-processing.

## Value

Depending on `scope`, either a numeric vector of averaged AR
coefficients, a list of per-run coefficient vectors, or the raw stored
structure. Returns `NULL` when no AR modeling was performed.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit model with AR(1) errors
fit <- fmri_lm(onset ~ hrf(cond), dataset = dset, cor_struct = "ar1")
ar_parameters(fit)  # Extract estimated AR(1) coefficient
} # }
```
