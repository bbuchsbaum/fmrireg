# Describe how to (re)construct a subject's dataset

A small, serializable recipe for a dataset: the name of a dataset
constructor plus the arguments to call it with. For file-backed data the
arguments are paths and run lengths (no voxel data), which keeps the
enclosing
[fmri_job](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md)
tiny and portable. The dataset is realized lazily on the worker.

## Usage

``` r
dataset_spec(constructor, args = list(), source = c("file", "inline"))
```

## Arguments

- constructor:

  Name of a dataset constructor (a string), e.g. `"fmri_dataset"` or
  `"matrix_dataset"`. Resolved at run time, so the data is not loaded
  when the spec is built.

- args:

  A named list of arguments passed to `constructor` (for
  `"fmri_dataset"`: `scans`, `TR`, `run_length`, `event_table`, `mask`,
  `base_path`, ...).

- source:

  Either `"file"` (paths; nothing loaded until run) or `"inline"` (data
  already in `args`, e.g. a `matrix_dataset`).

## Value

An object of class `dataset_spec`.

## See also

[`fmri_job()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md),
[`fmri_dataset()`](https://bbuchsbaum.github.io/fmridataset/reference/fmri_dataset.html),
[`matrix_dataset()`](https://bbuchsbaum.github.io/fmridataset/reference/matrix_dataset.html)

## Examples

``` r
dataset_spec("fmri_dataset",
             args = list(scans = c("run-1_bold.nii.gz", "run-2_bold.nii.gz"),
                         TR = 2, run_length = c(200, 200)),
             source = "file")
#> <dataset_spec>
#>   constructor: fmri_dataset()  source: file
#>   args: scans, TR, run_length
```
