# Perform Chunkwise Linear Modeling on fMRI Dataset

This function performs a chunkwise linear model analysis on an fMRI
dataset, splitting the dataset into chunks and running the linear model
on each chunk.

## Usage

``` r
# S3 method for class 'fmri_dataset_old'
chunkwise_lm(
  dset,
  model,
  contrast_objects,
  nchunks,
  cfg,
  verbose = FALSE,
  use_fast_path = FALSE,
  progress = FALSE,
  phi_fixed = NULL,
  sigma_fixed = NULL
)
```

## Arguments

- dset:

  An `fmri_dataset` object.

- model:

  The `fmri_model` used for the analysis.

- contrast_objects:

  The list of full contrast objects.

- nchunks:

  The number of chunks to divide the dataset into.

- cfg:

  An `fmri_lm_config` object containing all fitting options.

- verbose:

  Logical. Whether to display progress messages (default is `FALSE`).

- use_fast_path:

  Logical. If `TRUE`, use matrix-based computation for speed. Default is
  `FALSE`.

- progress:

  Logical. Display a progress bar for chunk processing. Default is
  `FALSE`.

- phi_fixed:

  Optional fixed AR parameters.

- sigma_fixed:

  Optional fixed robust scale estimate.

## Value

A list containing the unpacked chunkwise results.
