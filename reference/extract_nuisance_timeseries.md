# Extract Nuisance Timeseries from Mask

Extracts voxel timeseries from regions defined by a mask (e.g., WM/CSF).
This is the typical input for soft subspace projection.

## Usage

``` r
extract_nuisance_timeseries(dataset, mask, run = NULL)
```

## Arguments

- dataset:

  An fmri_dataset object.

- mask:

  A binary mask (logical vector or 3D array) indicating nuisance voxels,
  or a file path to a NIfTI mask.

- run:

  Optional run index to extract data from a specific run.

## Value

Matrix of nuisance timeseries (time x nuisance_voxels).
