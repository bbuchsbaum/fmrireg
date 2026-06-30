# Assemble the full fMRI model for a job

Builds the per-subject `baseline_model` (from the template's
[baseline_spec](https://bbuchsbaum.github.io/fmrireg/reference/baseline_spec.md),
the dataset's sampling frame, and any nuisance regressors) and combines
it with the event model into an `fmri_model`.

## Usage

``` r
build_model(job, dataset = realize_dataset(job))
```

## Arguments

- job:

  An
  [fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md).

- dataset:

  The realized dataset (defaults to `realize_dataset(job)`).

## Value

An `fmri_model`.

## See also

[`realize_dataset()`](https://bbuchsbaum.github.io/fmrireg/reference/realize_dataset.md),
`run()`
