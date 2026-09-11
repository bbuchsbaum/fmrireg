# Realize the dataset described by a job

Reconstructs the `fmri_frame` from the job's
[dataset_spec](https://bbuchsbaum.github.io/fmrireg/reference/dataset_spec.md).
For file-backed specs this is where data first becomes addressable
(still lazily, through the frame's NIfTI array source).

## Usage

``` r
realize_dataset(job)
```

## Arguments

- job:

  An
  [fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md).

## Value

An `fmri_frame`.

## See also

[`build_model()`](https://bbuchsbaum.github.io/fmrireg/reference/build_model.md),
[`run_job()`](https://bbuchsbaum.github.io/fmrireg/reference/run_job.md)
