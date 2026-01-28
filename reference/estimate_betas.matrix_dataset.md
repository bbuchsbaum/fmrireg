# Estimate betas for a matrix dataset

This function estimates betas (regression coefficients) for fixed and
random effects in a matrix dataset using various methods.

## Usage

``` r
# S3 method for class 'matrix_dataset'
estimate_betas(
  x,
  fixed = NULL,
  ran,
  block,
  method = c("mixed", "lss", "ols"),
  basemod = NULL,
  fracs = 0.5,
  progress = TRUE,
  ...
)
```

## Arguments

- x:

  An object of class `matrix_dataset` representing the matrix dataset

- fixed:

  A formula specifying the fixed regressors that model constant effects
  (i.e., non-varying over trials)

- ran:

  A formula specifying the random (trialwise) regressors that model
  single trial effects

- block:

  A formula specifying the block factor

- method:

  The regression method for estimating trialwise betas; one of "mixed",
  "lss", or "ols" (default: "mixed")

- basemod:

  A `baseline_model` instance to regress out of data before beta
  estimation (default: NULL)

- fracs:

  Fraction of voxels used for prewhitening.

- progress:

  Logical; show progress bar.

- ...:

  Additional arguments passed to the estimation method

## Value

A list of class "fmri_betas" containing the following components:

- betas_fixed: Matrix representing the fixed effect betas

- betas_ran: Matrix representing the random effect betas

- design_ran: Design matrix for random effects

- design_fixed: Design matrix for fixed effects

- design_base: Design matrix for baseline model

## See also

[`matrix_dataset`](https://bbuchsbaum.github.io/fmridataset/reference/matrix_dataset.html),
[`baseline_model`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html)

Other estimate_betas:
[`estimate_betas()`](https://bbuchsbaum.github.io/fmrireg/reference/estimate_betas.md)
