# Fast fMRI Regression Model Estimation from a Latent Component Dataset

This function estimates a regression model for fMRI data using a latent
component dataset. The dataset must be an `fmri_frame` whose feature
space is a
[`fmridataset::basis_space`](https://bbuchsbaum.github.io/fmridataset/reference/basis_space.html)
(see
[`latent_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/latent_frame.md)),
so that the assay holds component scores and the space carries the
spatial loadings.

## Usage

``` r
fmri_latent_lm(
  formula,
  block,
  baseline_model = NULL,
  dataset,
  durations = 0,
  drop_empty = TRUE,
  robust = FALSE,
  autocor = c("none", "auto", "ar1", "ar2", "arma"),
  bootstrap = FALSE,
  nboot = 1000,
  control = NULL,
  ...
)
```

## Arguments

- formula:

  A formula specifying the regression model.

- block:

  A factor indicating the block structure of the data.

- baseline_model:

  An optional baseline model.

- dataset:

  A latent `fmri_frame` (see
  [`latent_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/latent_frame.md)).

- durations:

  The duration of events in the dataset. Default is 0 (impulse events),
  as in
  [`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md).

- drop_empty:

  Whether to drop empty events from the model. Default is TRUE.

- robust:

  Must be FALSE. Robust weights are voxel-specific and cannot be
  computed on component scores, so any robust request is an error.

- autocor:

  The autocorrelation model used when `control` is `NULL`. 'none' (OLS),
  'ar1' or 'ar2'. AR models use one set of AR coefficients per run
  shared across all components, estimated so that the fit equals a dense
  fit of the reconstructed voxel data. 'auto' and 'arma' are not
  implemented and are an error. When `control` is supplied, `autocor`
  may be omitted; if it is also given it must name the same structure as
  `control`.

- bootstrap:

  Must be FALSE; bootstrapping is not implemented for latent datasets
  and TRUE is an error.

- nboot:

  Unused (see `bootstrap`).

- control:

  Optional
  [`fmri_lm_control()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md)
  object. `NULL` (the default) fits `estimation_spec("joint")` with the
  noise model named by `autocor`, shared across components. Supply a
  control to choose the estimation scope and noise model explicitly,
  e.g.
  `fmri_lm_control(estimation = estimation_spec("runwise_meta"), noise = noise_spec("ar1", voxelwise = TRUE))`
  for a separate AR model per component (see Details).

- ...:

  Not used. Any argument passed here is an error, so a misspelled or
  unsupported argument is never silently dropped.

## Value

An object of class 'fmri_latent_lm' containing the regression model and
dataset.

## Details

**Default (shared noise) model.** With `control = NULL` the model is
fitted jointly over runs and, for AR noise, uses one AR model per run
shared by all components. Because voxel data are `Y = Z %*% t(L)`
(scores `Z`, loadings `L`) this is exactly the fit a dense GLS of the
reconstructed voxel data would give, so voxel-level standard errors and
statistics are available with `recon = TRUE`.

**Per-component noise.** Latent components can differ strongly in
temporal autocorrelation. A `control` with
`estimation_spec("runwise_meta")` (and, optionally,
`noise_spec(..., voxelwise = TRUE)`) treats each component as a separate
response with its own noise model and run weights, i.e. the same fit
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
gives for the score matrix as a
[`matrix_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md).
That model has no exact voxel-space equivalent: `coef(recon = TRUE)`
returns the loadings applied to the component estimates, but
[`standard_error()`](https://bbuchsbaum.github.io/fmrireg/reference/standard_error.md)
and [`stats()`](https://bbuchsbaum.github.io/fmrireg/reference/stats.md)
with `recon = TRUE` are an error because the component noise models do
not combine into one voxel covariance. Under joint estimation only
shared (run- or global-pooled) AR models are supported.

## Note

This method is currently experimental.

## See also

[`latent_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/latent_frame.md),
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md),
[`fmri_lm_control()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md)

## Examples

``` r


# Estimate the fMRI regression model using the latent dataset
#result <- fmri_latent_lm(formula = formula, block = block, dataset = dset,
#                          durations = NULL, drop_empty = TRUE, robust = FALSE)

# Per-component AR(1) noise, run-wise estimation
#ctrl <- fmri_lm_control(estimation = estimation_spec("runwise_meta"),
#                        noise = noise_spec("ar1", voxelwise = TRUE))
#result <- fmri_latent_lm(formula, block = block, dataset = dset,
#                         control = ctrl)

# Print the result
#print(result)
```
