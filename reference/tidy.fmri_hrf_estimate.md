# Tidy a smooth HRF estimate

Tidy a smooth HRF estimate

## Usage

``` r
# S3 method for class 'fmri_hrf_estimate'
tidy(x, curve = NULL, voxel = NULL, ...)
```

## Arguments

- x:

  An `fmri_hrf_estimate` object.

- curve:

  Optional curve names to retain.

- voxel:

  Optional voxel names or indices to retain.

- ...:

  Unused.

## Value

A tibble with one row per time, curve, and voxel.
