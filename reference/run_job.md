# Run a single job

Realizes the dataset, assembles the model, fits it with the template's
control options (and engine, if any), and applies the template's
reducer.

## Usage

``` r
run_job(job, progress = FALSE)
```

## Arguments

- job:

  An
  [fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md).

- progress:

  Logical; show a fitting progress bar.

## Value

The reducer's output, or the fitted `fmri_lm` if the template has no
reducer.

## Details

Template contrasts given as contrast specifications are computed in the
fit (see
[`build_model()`](https://bbuchsbaum.github.io/fmrireg/reference/build_model.md)).
Template contrasts given as weight vectors or matrices are computed
afterwards with
[`fit_contrasts()`](https://bbuchsbaum.github.io/fmrireg/reference/fit_contrasts.md)
and added to the fit's contrast table; post-hoc contrasts are not
available for multi-run `runwise_meta` fits with a temporal noise model
or robust weighting.

## See also

[`run_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/run_jobs.md),
[`build_model()`](https://bbuchsbaum.github.io/fmrireg/reference/build_model.md)
