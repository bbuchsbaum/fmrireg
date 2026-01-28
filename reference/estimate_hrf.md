# Estimate hemodynamic response function (HRF) using Generalized Additive Models (GAMs)

This function estimates the HRF using GAMs from the `mgcv` package. The
HRF can be estimated with or without fixed effects.

## Usage

``` r
estimate_hrf(
  form,
  fixed = NULL,
  block,
  dataset,
  bs = c("tp", "ts", "cr", "ps"),
  rsam = seq(0, 20, by = 1),
  basemod = NULL,
  k = 8,
  fx = TRUE,
  progress = TRUE
)
```

## Arguments

- form:

  A formula specifying the event model for the conditions of interest

- fixed:

  A formula specifying the fixed regressors that model constant effects
  (i.e., non-varying over trials); default is NULL

- block:

  A formula specifying the block factor

- dataset:

  An object representing the fMRI dataset

- bs:

  Basis function for the smooth term in the GAM; one of "tp" (default),
  "ts", "cr", or "ps"

- rsam:

  A sequence of time points at which the HRF is estimated (default:
  seq(0, 20, by = 1))

- basemod:

  A `baseline_model` instance to regress out of data before HRF
  estimation (default: NULL)

- k:

  the dimension of the basis, default is 8

- fx:

  indicates whether the term is a fixed d.f. regression spline (TRUE) or
  a penalized regression spline (FALSE); default is TRUE.

- progress:

  Logical; display progress during estimation.

## Value

A matrix with the estimated HRF values for each voxel

## See also

[`baseline_model`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html),
[`event_model`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html),
[`design_matrix`](https://bbuchsbaum.github.io/fmridesign/reference/design_matrix.html)

## Examples

``` r
# To be added
```
