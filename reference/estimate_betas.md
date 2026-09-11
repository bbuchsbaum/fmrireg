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
random effects using various regression methods including mixed models,
least squares, and PLS.

## Usage

``` r
estimate_betas(x, ...)

# S3 method for class 'fmri_frame'
estimate_betas(
  x,
  fixed = NULL,
  ran,
  block,
  method = c("mixed", "lss", "ols"),
  basemod = NULL,
  maxit = 1000,
  fracs = 0.5,
  progress = TRUE,
  ...
)
```

## Arguments

- x:

  An `fmri_frame` (see
  [`matrix_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md),
  [`neurovec_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/neurovec_frame.md),
  [`nifti_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/nifti_frame.md),
  and
  [`latent_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/latent_frame.md)).
  When the frame's feature space is a `volume_space`, the fixed and
  random betas are returned as `NeuroVec` objects on that grid;
  otherwise (index or basis spaces) they are returned as
  coefficient-by-feature matrices.

- ...:

  Additional arguments passed to the estimation method.

- fixed:

  A formula specifying the fixed regressors that model constant effects
  (i.e., non-varying over trials).

- ran:

  A formula specifying the random (trialwise) regressors that model
  single trial effects.

- block:

  A formula specifying the block factor.

- method:

  The regression method for estimating trialwise betas; one of "mixed",
  "lss", or "ols".

- basemod:

  A `baseline_model` instance to regress out of data before beta
  estimation (default: NULL).

- maxit:

  Maximum number of iterations for optimization methods (default: 1000).

- fracs:

  Fraction of voxels used for prewhitening.

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

- betas_fixed: fixed effect betas (NeuroVec for volumetric frames,
  matrix otherwise).

- betas_ran: random effect betas (NeuroVec for volumetric frames, matrix
  otherwise).

- design_ran: Design matrix for random effects.

- design_fixed: Design matrix for fixed effects.

- design_base: Design matrix for baseline model.

- basemod: Baseline model object.

- fixed_model: Fixed effect model object.

- ran_model: Random effect model object.

- estimated_hrf: The estimated HRF vector (NULL for most methods).

## Details

This is a generic function whose `fmri_frame` method adapts to the
frame's feature space: volumetric frames (`volume_space`) return
`NeuroVec` betas, while matrix-format (`index_space`) and latent
(`basis_space`) frames return coefficient matrices.

Available estimation methods include:

- mixed:

  Mixed-effects model using ridge/BLUP estimation

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

[`matrix_frame`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md),
[`neurovec_frame`](https://bbuchsbaum.github.io/fmrireg/reference/neurovec_frame.md),
[`latent_frame`](https://bbuchsbaum.github.io/fmrireg/reference/latent_frame.md)

[`matrix_frame`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md),
[`baseline_model`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html),
[`event_model`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html)

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
dset <- matrix_frame(
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

if (FALSE) { # \dontrun{
facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
facedes$frun <- factor(facedes$run)
scans <- paste0("rscan0", 1:6, ".nii")

dset <- nifti_frame(scans=scans, mask="mask.nii", TR=1.5,
        run_length=rep(436,6), event_table=facedes)
fixed = onset ~ hrf(run)
ran = onset ~ trialwise()
block = ~ run

betas <- estimate_betas(dset, fixed=fixed, ran=ran, block=block, method="mixed")
} # }
```
