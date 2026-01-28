# Estimate betas using various regression methods

This function estimates betas (regression coefficients) for fixed and
random effects using various regression methods including mixed models,
least squares, and PLS.

## Usage

``` r
# S3 method for class 'fmri_dataset'
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

  An object of class `fmri_dataset` representing the fMRI dataset.

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

- ...:

  Additional arguments passed to the estimation method.

## Value

A list of class "fmri_betas" containing the following components:

- betas_fixed: NeuroVec object representing the fixed effect betas.

- betas_ran: NeuroVec object representing the random effect betas.

- design_ran: Design matrix for random effects.

- design_fixed: Design matrix for fixed effects.

- design_base: Design matrix for baseline model.

- basemod: Baseline model object.

- fixed_model: Fixed effect model object.

- ran_model: Random effect model object.

- estimated_hrf: The estimated HRF vector (NULL for most methods).

## See also

[`fmri_dataset`](https://bbuchsbaum.github.io/fmridataset/reference/fmri_dataset.html),
[`baseline_model`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html),
[`event_model`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html)

## Examples

``` r
if (FALSE) { # \dontrun{
facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
facedes$frun <- factor(facedes$run)
scans <- paste0("rscan0", 1:6, ".nii")

dset <- fmri_dataset(scans=scans, mask="mask.nii", TR=1.5, 
        run_length=rep(436,6), event_table=facedes)
fixed = onset ~ hrf(run)
ran = onset ~ trialwise()
block = ~ run

betas <- estimate_betas(dset, fixed=fixed, ran=ran, block=block, method="mixed")
} # }
```
