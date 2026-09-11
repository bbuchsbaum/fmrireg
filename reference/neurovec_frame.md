# Build an fmri_frame from in-memory NeuroVec objects

Extracts the masked time series of one or more
[`neuroim2::NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
runs and builds an
[`fmridataset::fmri_frame()`](https://bbuchsbaum.github.io/fmridataset/reference/fmri_frame.html)
whose feature space is the
[`fmridataset::volume_space()`](https://bbuchsbaum.github.io/fmridataset/reference/volume_space.html)
described by `mask`. Voxel time series are read with
[`neuroim2::series()`](https://bbuchsbaum.github.io/neuroim2/reference/series-methods.html),
so a `SparseNeuroVec` does not need to be densified.

## Usage

``` r
neurovec_frame(
  scans,
  mask,
  TR,
  run_length = NULL,
  event_table = NULL,
  censor = NULL,
  assay = "bold"
)
```

## Arguments

- scans:

  A `NeuroVec`, or a list of them, one per run.

- mask:

  A
  [`neuroim2::LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html)
  (or a `NeuroVol` whose non-zero voxels form the mask) on the same 3D
  grid as `scans`.

- TR:

  Repetition time in seconds; one value, or one per run.

- run_length:

  Volumes per run; defaults to the fourth dimension of each scan.

- event_table, censor:

  As in
  [`matrix_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md).

- assay:

  Name of the assay.

## Value

An `fmri_frame` over a `volume_space`.

## See also

[`matrix_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md),
[`nifti_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/nifti_frame.md)
