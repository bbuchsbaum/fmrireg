# Tidy fitted HRF curves from an fmri_lm fit

Convert
[`fitted_hrf()`](https://bbuchsbaum.github.io/fmrireg/reference/fitted_hrf.md)
output into a long tibble suitable for plotting.

## Usage

``` r
tidy_fitted_hrf(
  x,
  sample_at = seq(0, 24, by = 1),
  term = NULL,
  term_match = c("contains", "exact", "regex"),
  voxel = 1L,
  average_voxels = FALSE,
  ...
)
```

## Arguments

- x:

  An `fmri_lm` object.

- sample_at:

  Numeric vector of time points passed to
  [`fitted_hrf()`](https://bbuchsbaum.github.io/fmrireg/reference/fitted_hrf.md).

- term:

  Optional term selector. When `NULL`, all event terms are returned.

- term_match:

  Matching mode for `term`: `"contains"` (default), `"exact"`, or
  `"regex"`.

- voxel:

  Integer voxel index to extract from each term's prediction matrix.

- average_voxels:

  Logical; if `TRUE`, average across voxels instead of selecting one
  voxel.

- ...:

  Additional arguments passed to
  [`fitted_hrf()`](https://bbuchsbaum.github.io/fmrireg/reference/fitted_hrf.md).

## Value

A tibble with columns `term`, `time`, `condition`, `estimate`, and
`voxel`.
