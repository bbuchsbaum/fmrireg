# Construct an fMRI Regression Model

This function constructs an fMRI regression model consisting of an event
model and a baseline model. The resulting model can be used for the
analysis of fMRI data.

## Usage

``` r
fmri_model(event_model, baseline_model, dataset)
```

## Arguments

- event_model:

  An object of class "event_model" representing the event-related part
  of the fMRI regression model.

- baseline_model:

  An object of class "baseline_model" representing the baseline-related
  part of the fMRI regression model.

- dataset:

  An `fmri_dataset` used to build the model.

## Value

An object of class `fmri_model` containing the event and baseline models
along with the dataset.

## See also

event_model, baseline_model
