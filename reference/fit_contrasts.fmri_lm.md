# Fit Contrasts for fMRI Linear Model Objects

S3 method for computing contrasts on fitted fmri_lm objects.

## Usage

``` r
# S3 method for class 'fmri_lm'
fit_contrasts(object, contrasts, ...)
```

## Arguments

- object:

  An fmri_lm object

- contrasts:

  A list of contrast specifications

- ...:

  Additional arguments (unused)

## Value

A list of contrast results

## Details

Multi-run fits estimated run by run (`runwise_meta`) with temporal noise
models, robust weighting, volume weights or projection whiten or weight
each run differently, so their pooled betas have no single stored
covariance; post-hoc contrasts are refused for them. Specify such
contrasts in the model formula, or fit with
`estimation_spec(scope = "joint")`.
