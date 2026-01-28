# Perform Runwise Linear Modeling on fMRI Dataset

This function performs a runwise linear model analysis on an fMRI
dataset, running the linear model on each run separately and then
pooling results.

This function performs a runwise linear model analysis on an fMRI
dataset by running the linear model for each data run and combining the
results.

## Usage

``` r
runwise_lm(
  dset,
  model,
  contrast_objects,
  cfg,
  verbose = FALSE,
  use_fast_path = FALSE,
  progress = FALSE,
  phi_fixed = NULL,
  sigma_fixed = NULL,
  parallel_voxels = FALSE
)

runwise_lm(
  dset,
  model,
  contrast_objects,
  cfg,
  verbose = FALSE,
  use_fast_path = FALSE,
  progress = FALSE,
  phi_fixed = NULL,
  sigma_fixed = NULL,
  parallel_voxels = FALSE
)
```

## Arguments

- dset:

  An `fmri_dataset` object.

- model:

  The `fmri_model` used for the analysis.

- contrast_objects:

  The list of full contrast objects.

- cfg:

  An `fmri_lm_config` object containing all fitting options.

- verbose:

  Logical. Whether to display progress messages (default is `FALSE`).

- use_fast_path:

  Logical. Whether to use fast path computation (default is `FALSE`).

- progress:

  Logical. Display a progress bar for run processing. Default is
  `FALSE`.

- phi_fixed:

  Optional fixed AR parameters.

- sigma_fixed:

  Optional fixed robust scale estimate.

- parallel_voxels:

  Logical. If TRUE, process voxels in parallel using `future.apply`.
  Default is `FALSE`.

## Value

A list containing the combined results from runwise linear model
analysis.

A list containing the combined results from runwise linear model
analysis.
