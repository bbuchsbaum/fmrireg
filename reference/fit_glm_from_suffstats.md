# Fit GLM from sufficient statistics

Computes OLS/GLS-equivalent estimates from cross-products without
materializing the transformed series. Engines can stream `XtX`, `XtS`,
and `StS` and call this to obtain a standard `fmri_lm` object. AR/robust
are not estimated here; this is a low-level helper for OLS-equivalent
inference from suffstats.

## Usage

``` r
fit_glm_from_suffstats(
  model,
  XtX,
  XtS,
  StS,
  df,
  cfg = NULL,
  dataset = NULL,
  strategy = "external",
  engine = "external"
)
```

## Arguments

- model:

  An `fmri_model` describing the design.

- XtX:

  p×p cross-product of the design matrix.

- XtS:

  p×V cross-product of design with data.

- StS:

  length-V vector of sum of squares per voxel.

- df:

  Residual degrees of freedom.

- cfg:

  Optional `fmri_lm_config`; used for metadata only.

- dataset:

  Optional dataset backing the model.

- strategy:

  Character label for the returned object.

- engine:

  Character label for the returned object.

## Value

An object of class `fmri_lm`.
