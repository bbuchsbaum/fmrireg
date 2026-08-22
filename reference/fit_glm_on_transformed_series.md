# Fit the core GLM on a transformed time-series matrix

This helper performs IID, non-robust OLS inference. It rejects controls
that request AR or robust fitting; use
[`fit_glm_with_config()`](https://bbuchsbaum.github.io/fmrireg/reference/fit_glm_with_config.md)
for those models.

## Usage

``` r
fit_glm_on_transformed_series(
  model,
  Y,
  cfg = NULL,
  dataset = NULL,
  strategy = "external",
  engine = "external"
)
```

## Arguments

- model:

  An `fmri_model` describing the design.

- Y:

  Numeric matrix with `nrow(Y)` time points and columns matching voxels.

- cfg:

  Optional `fmri_lm_control`; defaults to
  [`fmri_lm_control()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md)
  and must request IID noise and no robust fitting.

- dataset:

  Optional dataset whose feature axis describes the columns of `Y`. It
  is never inferred from `model$dataset`, because a transformed response
  can have a different feature space.

- strategy:

  Character label recorded on the returned object. Defaults to
  "external".

- engine:

  Character label indicating the engine that produced the fit. Defaults
  to "external".

## Value

An object of class `fmri_lm`.
