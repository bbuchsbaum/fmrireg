# Configuration for fmri_lm fitting

`fmri_lm_control()` creates an `fmri_lm_config` object collecting all
options for robust and autoregressive modelling. It validates inputs and
applies defaults so downstream functions receive a single structured
list.

## Usage

``` r
fmri_lm_control(
  robust_options = list(),
  ar_options = list(),
  volume_weights_options = list(),
  soft_subspace_options = list()
)
```

## Arguments

- robust_options:

  list of robust fitting options. See Details.

- ar_options:

  list of autoregressive modelling options. See Details.

- volume_weights_options:

  list of volume weighting options. See Details. For simple cases, use
  the `volume_weights` parameter in
  [`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
  instead.

- soft_subspace_options:

  list of soft subspace projection options. See Details. For simple
  cases, use the `nuisance_projection` parameter in
  [`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
  instead.

## Value

An object of class `fmri_lm_config`.

## Details

For common use cases,
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
provides convenience parameters that are easier to use than these
detailed option lists:

- `volume_weights = TRUE` enables volume weighting with defaults

- `volume_weights = "tukey"` enables with Tukey method

- `nuisance_projection = N` enables soft projection with matrix N

- `nuisance_projection = "mask.nii"` enables with mask file

Use the `*_options` lists below only when you need fine-grained control.

`robust_options` may contain:

- `type` (`FALSE`, "huber", "bisquare")

- `k_huber`

- `c_tukey`

- `max_iter`

- `scale_scope` ("run", "global")

- `reestimate_phi` (logical)

`ar_options` may contain:

- `struct` ("iid", "ar1", "ar2", "arp")

- `p` (order for "arp")

- `iter_gls` (integer number of GLS iterations)

- `global` (logical, use global phi)

- `voxelwise` (logical)

- `exact_first` (logical)

`volume_weights_options` may contain:

- `enabled` (logical, whether to compute and apply volume weights)

- `method` ("inverse_squared", "soft_threshold", "tukey")

- `threshold` (numeric, DVARS threshold for weighting)

- `weights` (optional pre-computed weight vector)

`soft_subspace_options` may contain:

- `enabled` (logical, whether to apply soft subspace projection)

- `nuisance_mask` (path to NIfTI mask or logical vector)

- `nuisance_matrix` (pre-computed nuisance timeseries matrix)

- `lambda` (numeric, "auto", or "gcv")

- `warn_redundant` (logical, warn if baseline has nuisance terms)
