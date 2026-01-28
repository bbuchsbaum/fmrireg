# Fit the core GLM on a transformed time-series matrix

Fit the core GLM on a transformed time-series matrix

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

  Optional `fmri_lm_config`; defaults to
  [`fmri_lm_control()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md).

- dataset:

  Optional dataset backing the model. Defaults to `model$dataset` when
  available.

- strategy:

  Character label recorded on the returned object. Defaults to
  "external".

- engine:

  Character label indicating the engine that produced the fit. Defaults
  to "external".

## Value

An object of class `fmri_lm`.
