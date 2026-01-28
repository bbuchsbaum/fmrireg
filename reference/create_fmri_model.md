# Create an fMRI Model

This function creates an fMRI model consisting of an event model and a
baseline model.

This function creates an `fmri_model` object from a formula, block
specification, and dataset. It's a convenience function that combines
event and baseline models.

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

create_fmri_model(
  formula,
  block,
  baseline_model = NULL,
  dataset,
  drop_empty = TRUE,
  durations = 0
)

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

  An `fmri_dataset` containing the event table and sampling frame.

- drop_empty:

  Logical. Whether to remove factor levels with zero size. Default is
  `TRUE`.

- durations:

  A vector of event durations. Default is `0`.

## Value

An `fmri_model` object.

An fmri_model object that also stores the `dataset`

An `fmri_model` object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have an fmri_dataset object named ds and a formula for events:
fmri_mod <- create_fmri_model(formula = onset ~ hrf(x) + hrf(y),
                              block = ~ run,
                              dataset = ds,
                              drop_empty = TRUE,
                              durations = rep(0, nrow(ds$event_table)))
} # }
```
