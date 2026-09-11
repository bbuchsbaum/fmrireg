# Build a lazy fmri_frame from NIfTI files

Wraps
[`fmridataset::nifti_array_source()`](https://bbuchsbaum.github.io/fmridataset/reference/nifti_array_source.html),
so headers and the mask are read at construction but no BOLD volumes are
read until the fit needs them. For fMRIPrep derivatives see
[`fmridataset::read_bids_bold()`](https://bbuchsbaum.github.io/fmridataset/reference/read_bids_bold.html).

## Usage

``` r
nifti_frame(
  scans,
  mask,
  TR,
  run_length = NULL,
  event_table = NULL,
  censor = NULL,
  base_path = ".",
  assay = "bold"
)
```

## Arguments

- scans:

  Character vector of NIfTI paths, one per run, all on the same grid.

- mask:

  A NIfTI mask path, a
  [`neuroim2::LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html),
  or a
  [`fmridataset::volume_space()`](https://bbuchsbaum.github.io/fmridataset/reference/volume_space.html).

- TR:

  Repetition time in seconds; one value, or one per run.

- run_length:

  Volumes per run; defaults to the fourth header dimension of each file.

- event_table, censor:

  As in
  [`matrix_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md).

- base_path:

  Directory that relative `scans` and `mask` paths are resolved against.

- assay:

  Name of the assay.

## Value

An `fmri_frame` over the NIfTI source's `volume_space`.

## See also

[`matrix_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md),
[`neurovec_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/neurovec_frame.md),
[`fmridataset::nifti_array_source()`](https://bbuchsbaum.github.io/fmridataset/reference/nifti_array_source.html)
