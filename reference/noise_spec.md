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

  Nonnegative moving-average order. Positive values request an
  ARMA(p, q) model, with the AR order determined by `struct` and `p`.
  The current MA-capable path supports runwise meta-estimation with run
  pooling, without censoring, parcel pooling, or voxelwise covariance
  estimation.

- iter_gls:

  Number of GLS refinement iterations.

- pooling:

  Temporal covariance pooling scope. With built-in shared AR estimation,
  `"run"` estimates one coefficient vector per run from the cross-voxel
  mean residual series; `"global"` concatenates those runwise mean
  series before estimation.

- parcels:

  Optional parcel labels for parcel pooling.

- voxelwise:

  Estimate temporal covariance separately by voxel. In the built-in
  fitter this currently requires AR-only runwise meta-estimation,
  `pooling = "run"`, `iter_gls = 1`, no censoring, and no volume
  weighting or soft-subspace projection. Robust fitting is supported,
  but robust AR re-estimation is not. Registered engines may define
  broader capabilities.

- exact_first:

  Use exact first-observation AR scaling.

- censor:

  Optional censor indices, logical mask, or `"auto"`.

- shrink_c0:

  Parcel shrinkage constant used by supporting engines.
