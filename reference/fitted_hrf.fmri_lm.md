# Compute Fitted Hemodynamic Response Functions for an fmri_lm Object

This method computes the fitted hemodynamic response functions (HRFs)
for an `fmri_lm` object.

## Usage

``` r
# S3 method for class 'fmri_lm'
fitted_hrf(x, sample_at = seq(0, 24, by = 1), ...)
```

## Arguments

- x:

  An `fmri_lm` object for which the fitted HRFs should be computed.

- sample_at:

  A numeric vector of time points at which the HRFs should be sampled.
  Default is `seq(0, 24, by = 1)`.

- ...:

  Additional arguments (currently unused).

## Value

A list where each element corresponds to an event term in the `fmri_lm`
object. Each element contains:

- `pred`:

  A matrix of predicted HRF values.

- `design`:

  A tibble containing the design matrix for the HRFs.
