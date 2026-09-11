# Build an fmri_frame from a latent (basis) decomposition

A
[`fmristore::LatentNeuroVec`](https://rdrr.io/pkg/fmrilatent/man/LatentNeuroVec.html)
stores temporal scores (time by component) and spatial loadings (voxel
by component). This helper places the scores in the assay and the
loadings in a
[`fmridataset::basis_space()`](https://bbuchsbaum.github.io/fmridataset/reference/basis_space.html)
whose parent is the mask's `volume_space`, so fmrireg's latent engines
can fit in component space and reconstruct voxel maps through
[`fmridataset::basis_synthesis()`](https://bbuchsbaum.github.io/fmridataset/reference/basis-operators.html).

## Usage

``` r
latent_frame(
  x,
  TR,
  run_length = NULL,
  event_table = NULL,
  censor = NULL,
  assay = "scores"
)
```

## Arguments

- x:

  A
  [`fmristore::LatentNeuroVec`](https://rdrr.io/pkg/fmrilatent/man/LatentNeuroVec.html),
  or a list with elements `basis` (time by component scores), `loadings`
  (masked-voxel by component), and `mask` (a
  [`neuroim2::LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html)).

- TR:

  Repetition time in seconds; one value, or one per run.

- run_length:

  Volumes per run; defaults to a single run spanning all score rows.

- event_table, censor:

  As in
  [`matrix_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md).

- assay:

  Name of the assay.

## Value

An `fmri_frame` over a synthesis-only `basis_space`.

## See also

[`matrix_frame()`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md),
[`fmridataset::basis_space_from_decoder()`](https://bbuchsbaum.github.io/fmridataset/reference/basis_space_from_decoder.html)
