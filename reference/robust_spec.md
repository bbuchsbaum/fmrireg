# Robust-estimation specification

Robust-estimation specification

## Usage

``` r
robust_spec(
  type = c("none", "huber", "bisquare"),
  k_huber = 1.345,
  c_tukey = 4.685,
  max_iter = 2L,
  scale_scope = c("run", "global", "voxel"),
  reestimate_phi = FALSE
)
```

## Arguments

- type:

  Robust loss, or `"none"` for ordinary least squares.

- k_huber:

  Tuning constant for the Huber loss.

- c_tukey:

  Tuning constant for Tukey's bisquare loss.

- max_iter:

  Maximum number of robust reweighting iterations.

- scale_scope:

  Scope at which robust scale is estimated.

- reestimate_phi:

  Re-estimate temporal-noise parameters after robust reweighting.
