# Fit a Linear Regression Model for fMRI Data Analysis

`fmri_lm` is a generic for fitting fMRI regression models. The default
interface accepts a model formula and dataset. An alternative method can
be used with a preconstructed `fmri_model` object that already contains
the design and data.

The public fitting boundary separates statistical choices in
[`fmri_lm_control()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md)
from mechanical execution choices in
[`compute_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/compute_spec.md).
Optional engines use the same control object and receive their private
arguments through `engine_args`.

## Usage

``` r
fmri_lm(formula, ...)

# S3 method for class 'formula'
fmri_lm(
  formula,
  block,
  baseline_model = NULL,
  dataset,
  durations = 0,
  drop_empty = TRUE,
  control = fmri_lm_control(estimation = estimation_spec("runwise_meta")),
  compute = compute_spec(),
  engine = NULL,
  engine_args = list(),
  lowrank = NULL,
  ...
)

# S3 method for class 'fmri_model'
fmri_lm(
  formula,
  dataset = NULL,
  control = fmri_lm_control(estimation = estimation_spec("runwise_meta")),
  compute = compute_spec(),
  engine = NULL,
  engine_args = list(),
  lowrank = NULL,
  ...
)
```

## Arguments

- formula:

  A model formula describing the event structure or an `fmri_model`
  object.

- ...:

  Transitional legacy arguments retained for one compatibility window.
  New code must express statistical choices in `control` and execution
  choices in `compute`; unknown arguments are rejected.

- block:

  Formula describing run/block structure.

- baseline_model:

  Optional baseline/nuisance model.

- dataset:

  An `fmri_dataset`. For an `fmri_model` method this may be omitted when
  the model already owns its dataset.

- durations:

  Event durations passed to model construction.

- drop_empty:

  Remove empty factor levels during model construction.

- control:

  A validated
  [`fmri_lm_control()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md).

- compute:

  A validated
  [`compute_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/compute_spec.md).

- engine:

  Optional registered engine name.

- engine_args:

  Named list passed only to `engine`.

- lowrank:

  Optional
  [`lowrank_control()`](https://bbuchsbaum.github.io/fmrireg/reference/lowrank_control.md)
  for the latent-sketch engine.

## Value

An object of class `fmri_lm`.
