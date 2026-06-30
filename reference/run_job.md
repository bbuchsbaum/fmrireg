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

## See also

[`run_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/run_jobs.md),
[`build_model()`](https://bbuchsbaum.github.io/fmrireg/reference/build_model.md)
