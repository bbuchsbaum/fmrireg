# Temporal-noise specification

Temporal-noise specification

## Usage

``` r
noise_spec(
  struct = c("iid", "ar1", "ar2", "arp"),
  p = NULL,
  q = 0L,
  iter_gls = 1L,
  pooling = c("run", "global", "parcel"),
  parcels = NULL,
  voxelwise = FALSE,
  exact_first = FALSE,
  censor = NULL,
  shrink_c0 = 100L
)
```

## Arguments

- struct:

  One of iid, ar1, ar2, or arp.

- p:

  AR order for `struct = "arp"`.

- q:

  Reserved MA order. Nonzero MA order is rejected until the fitter
  exposes an MA-capable inference path.

- iter_gls:

  Number of GLS refinement iterations.

- pooling:

  Temporal covariance pooling scope.

- parcels:

  Optional parcel labels for parcel pooling.

- voxelwise:

  Estimate temporal covariance separately by voxel.

- exact_first:

  Use exact first-observation AR scaling.

- censor:

  Optional censor indices, logical mask, or `"auto"`.

- shrink_c0:

  Parcel shrinkage constant used by supporting engines.
