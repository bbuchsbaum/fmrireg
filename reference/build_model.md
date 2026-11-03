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

## Details

Each contrast specification in the template's `contrasts` (e.g. the
elements of a
[`contrast_set()`](https://bbuchsbaum.github.io/fmrireg/reference/contrast_set.md))
is attached to the event term it resolves against, exactly as if it had
been written in that term (`hrf(..., contrasts = )`). The fit therefore
computes it like any formula contrast, under every estimation scope, and
[`write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/write_results.md)
/
[`reduce_write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_write_results.md)
write it.

## See also

[`realize_dataset()`](https://bbuchsbaum.github.io/fmrireg/reference/realize_dataset.md),
[`run_job()`](https://bbuchsbaum.github.io/fmrireg/reference/run_job.md)
