# Realize the dataset described by a job

Reconstructs the `fmri_dataset` from the job's
[dataset_spec](https://bbuchsbaum.github.io/fmrireg/reference/dataset_spec.md).
For file-backed specs this is where data first becomes addressable
(still lazily, per the dataset backend).

## Usage

``` r
realize_dataset(job)
```

## Arguments

- job:

  An
  [fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md).

## Value

An `fmri_dataset`.

## See also

[`build_model()`](https://bbuchsbaum.github.io/fmrireg/reference/build_model.md),
`run()`
