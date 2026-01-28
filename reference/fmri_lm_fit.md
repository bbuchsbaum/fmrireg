# Fit an fMRI Linear Regression Model with a Specified Fitting Strategy

This function fits an fMRI linear regression model using the specified
`fmri_model` object, dataset, and data splitting strategy (either
`"runwise"` or `"chunkwise"`). It is primarily an internal function used
by the `fmri_lm` function.

## Usage

``` r
fmri_lm_fit(
  fmrimod,
  dataset,
  strategy = c("runwise", "chunkwise"),
  cfg,
  nchunks = 10,
  use_fast_path = FALSE,
  progress = FALSE,
  parallel_voxels = FALSE,
  ...
)
```

## Arguments

- fmrimod:

  An `fmri_model` object.

- dataset:

  An `fmri_dataset` object containing the time-series data.

- strategy:

  The data splitting strategy, either `"runwise"` or `"chunkwise"`.
  Default is `"runwise"`.

- cfg:

  An `fmri_lm_config` object containing all fitting options. See
  [`fmri_lm_control`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm_control.md).

- nchunks:

  Number of data chunks when strategy is `"chunkwise"`. Default is `10`.

- use_fast_path:

  Logical. If `TRUE`, use matrix-based computation for speed. Default is
  `FALSE`.

- progress:

  Logical. Whether to display a progress bar during model fitting.
  Default is `FALSE`.

- parallel_voxels:

  Logical. If TRUE, voxelwise AR processing within runs is parallelised
  using `future.apply`. Default is `FALSE`.

- ...:

  Additional arguments.

## Value

A fitted fMRI linear regression model with the specified fitting
strategy.

## See also

[`fmri_lm`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md),
[`fmri_model`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_model.md),
[`fmri_dataset`](https://bbuchsbaum.github.io/fmridataset/reference/fmri_dataset.html)
