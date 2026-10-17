# Fit an fMRI Linear Regression Model with a Specified Fitting Strategy

This function fits an fMRI linear regression model using the specified
`fmri_model` object and dataset. The estimation scope in `cfg` selects
the strategy (`"runwise_meta"` fits runs separately and pools them;
`"joint"` fits one chunkwise model) and `compute` the mechanics. It is
primarily an internal function used by `fmri_lm`.

## Usage

``` r
fmri_lm_fit(fmrimod, dataset, cfg, compute = compute_spec())
```

## Arguments

- fmrimod:

  An `fmri_model` object.

- dataset:

  An `fmri_frame` object containing the time-series data.

- cfg:

  An `fmri_lm_control` object containing all fitting options. See
  [`fmri_lm_control`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md).

- compute:

  A
  [`compute_spec`](https://bbuchsbaum.github.io/fmrireg/reference/compute_spec.md):
  voxel chunks, matrix or reference backend, parallelism and progress.

## Value

A fitted fMRI linear regression model with the specified fitting
strategy.

## See also

[`fmri_lm`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md),
[`fmri_model`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_model.md),
[`matrix_frame`](https://bbuchsbaum.github.io/fmrireg/reference/matrix_frame.md)
