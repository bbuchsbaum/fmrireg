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
  shared_estimator = c("pooled_acvf", "mean_series"),
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

  Maximum number of GLS refinement iterations. Standard design-corrected
  AR fitting uses the initial OLS residuals once and holds that estimate
  fixed for the GLS solve. Iterative refinement remains available to
  ARMA models (`q > 0`).

- pooling:

  Temporal covariance pooling scope. With built-in shared AR estimation,
  `"run"` estimates one coefficient vector per run and `"global"` pools
  across runs.

- shared_estimator:

  Spatial estimator for a shared temporal covariance. `"pooled_acvf"`
  (the default) pools residual autocovariances across voxels and
  therefore targets a typical voxel covariance. `"mean_series"` first
  averages residual values across voxels and targets the coherent
  spatial component; it is experimental, can be much more autocorrelated
  than an individual voxel. Because OLS projection is linear across
  response columns, the matching design correction remains valid after
  averaging.

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
