# Estimate Beta Coefficients for fMRI Data

Estimate beta coefficients (regression parameters) from fMRI data using
various methods. This function supports different estimation approaches
for:

- single:

  Single-trial beta estimation

- effects:

  Fixed and random effects

- regularization:

  Various regularization techniques

- hrf:

  Optional HRF estimation

This function estimates betas (regression coefficients) for fixed and
random effects in a matrix dataset using various methods.

## Usage

``` r
estimate_betas(x, ...)

# S3 method for class 'latent_dataset'
estimate_betas(
  x,
  fixed = NULL,
  ran,
  block,
  method = c("mixed", "lss", "ols"),
  basemod = NULL,
  prewhiten = FALSE,
  progress = TRUE,
  ...
)
```

## Arguments

- x:

  An object of class `matrix_dataset` representing the matrix dataset

- ...:

  Additional arguments passed to the estimation method

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

- prewhiten:

  currently experimental, default to `FALSE`.

- progress:

  Logical; show progress bar.

## Value

A list of class "fmri_betas" containing:

- betas_fixed:

  Fixed effect coefficients

- betas_ran:

  Random (trial-wise) coefficients

- design_ran:

  Design matrix for random effects

- design_fixed:

  Design matrix for fixed effects

- design_base:

  Design matrix for baseline model

- method_specific:

  Additional components specific to the estimation method used

A list of class "fmri_betas" containing the following components:

- betas_fixed: Matrix representing the fixed effect betas

- betas_ran: Matrix representing the random effect betas

- design_ran: Design matrix for random effects

- design_fixed: Design matrix for fixed effects

- design_base: Design matrix for baseline model

## Details

This is a generic function with methods for different dataset types:

- fmri_dataset:

  For volumetric fMRI data

- matrix_dataset:

  For matrix-format data

- latent_dataset:

  For dimensionality-reduced data

Available estimation methods include:

- mixed:

  Mixed-effects model using rrBLUP

- r1:

  Rank-1 GLM with joint HRF estimation

- lss:

  Least-squares separate estimation

- pls:

  Partial least squares regression

- ols:

  Ordinary least squares

## References

Mumford, J. A., et al. (2012). Deconvolving BOLD activation in
event-related designs for multivoxel pattern classification analyses.
NeuroImage, 59(3), 2636-2643.

Pedregosa, F., et al. (2015). Data-driven HRF estimation for encoding
and decoding models. NeuroImage, 104, 209-220.

## See also

[`fmri_dataset`](https://bbuchsbaum.github.io/fmridataset/reference/fmri_dataset.html),
[`matrix_dataset`](https://bbuchsbaum.github.io/fmridataset/reference/matrix_dataset.html),
[`latent_dataset`](https://bbuchsbaum.github.io/fmridataset/reference/latent_dataset.html)

[`matrix_dataset`](https://bbuchsbaum.github.io/fmridataset/reference/matrix_dataset.html),
[`baseline_model`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html)

Other estimate_betas:
[`estimate_betas.matrix_dataset()`](https://bbuchsbaum.github.io/fmrireg/reference/estimate_betas.matrix_dataset.md)

## Examples

``` r
# Create example data
event_data <- data.frame(
  condition = factor(c("A", "B", "A", "B")),
  onset = c(1, 10, 20, 30),
  run = c(1, 1, 1, 1)
)

# Create sampling frame and dataset
sframe <- sampling_frame(blocklens = 100, TR = 2)
dset <- fmridataset::matrix_dataset(
  matrix(rnorm(100 * 2), 100, 2),
  TR = 2,
  run_length = 100,
  event_table = event_data
)

# Estimate betas using mixed-effects model
betas <- estimate_betas(
  dset,
  fixed = onset ~ hrf(condition),
  ran = onset ~ trialwise(),
  block = ~run,
  method = "mixed"
)
```
