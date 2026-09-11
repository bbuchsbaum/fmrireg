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

  Name of a frame constructor (a string): one of `"nifti_frame"`,
  `"matrix_frame"`, `"neurovec_frame"`, or `"latent_frame"`. Resolved at
  run time, so the data is not loaded when the spec is built.

- args:

  A named list of arguments passed to `constructor` (for
  `"nifti_frame"`: `scans`, `TR`, `run_length`, `event_table`, `mask`,
  `base_path`, ...).

- source:

  Either `"file"` (paths; nothing loaded until run) or `"inline"` (data
  already in `args`, e.g. a matrix for `matrix_frame`).

## Value

An object of class `dataset_spec`.

## See also

[`fmri_job()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_job.md),
[`nifti_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/nifti_frame.md),
[`matrix_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md)

## Examples

``` r
dataset_spec("nifti_frame",
             args = list(scans = c("run-1_bold.nii.gz", "run-2_bold.nii.gz"),
                         TR = 2, run_length = c(200, 200)),
             source = "file")
#> <dataset_spec>
#>   constructor: nifti_frame()  source: file
#>   args: scans, TR, run_length
```
