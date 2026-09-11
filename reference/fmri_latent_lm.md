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
  durations,
  drop_empty = TRUE,
  robust = FALSE,
  autocor = c("none", "auto", "ar1", "ar2", "arma"),
  bootstrap = FALSE,
  nboot = 1000,
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

  The duration of events in the dataset.

- drop_empty:

  Whether to drop empty events from the model. Default is TRUE.

- robust:

  Whether to use robust regression methods. Default is FALSE.

- autocor:

  The autocorrelation correction method to use on components. One of
  'none', 'auto', 'ar1', 'ar2', or 'arma'. Default is 'none'.

- bootstrap:

  Whether to compute bootstrapped parameter estimates. Default is FALSE.

- nboot:

  The number of bootstrap iterations. Default is 1000.

- ...:

  Additional arguments.

## Value

An object of class 'fmri_latent_lm' containing the regression model and
dataset.

## Note

This method is currently experimental.

## Examples

``` r


# Estimate the fMRI regression model using the latent dataset
#result <- fmri_latent_lm(formula = formula, block = block, dataset = dset,
#                          durations = NULL, drop_empty = TRUE, robust = FALSE)

# Print the result
#print(result)
```
