# Internal handler for spatial FDR on fmri_meta objects

Internal handler for spatial FDR on fmri_meta objects

## Usage

``` r
spatial_fdr_fmri_meta(
  object,
  coef = 1,
  group = NULL,
  alpha = 0.05,
  tau = 0.5,
  lambda = 1,
  neighbors = NULL,
  min_pi0 = 0.05,
  empirical_null = TRUE,
  verbose = FALSE
)
```

## Arguments

- object:

  An fmri_meta object

- coef:

  Character or integer; coefficient name or index to test

- group:

  Grouping variable for spatial FDR. Can be:

  - NULL: Auto-create blocks (voxelwise) or use ROIs (ROI-wise)

  - Integer/factor vector: Custom grouping

  - "blocks": Create 3D blocks (voxelwise only)

  - "parcels": Use existing parcellation

- ...:

  Additional arguments passed to `spatial_fdr`
