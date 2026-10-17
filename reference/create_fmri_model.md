# Create an fMRI Model

This function creates an `fmri_model` by combining an event model and a
baseline model. If a baseline model is not provided, a default one is
created based on the dataset.

## Usage

``` r
create_fmri_model(
  formula,
  block,
  baseline_model = NULL,
  dataset,
  drop_empty = TRUE,
  durations = 0
)
```

## Arguments

- formula:

  The model formula for experimental events.

- block:

  The model formula for block structure.

- baseline_model:

  (Optional) A `baseline_model` object. If `NULL`, a default baseline
  model is created.

- dataset:

  An `fmri_frame` containing the event table and sampling frame.

- drop_empty:

  Logical. Whether to remove factor levels with zero size. Default is
  `TRUE`.

- durations:

  A vector of event durations. Default is `0`.

## Value

An `fmri_model` object.
