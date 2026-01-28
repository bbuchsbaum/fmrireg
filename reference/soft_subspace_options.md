# Soft Subspace Control Options

Creates a configuration object for soft subspace projection.

## Usage

``` r
soft_subspace_options(
  enabled = FALSE,
  nuisance_mask = NULL,
  nuisance_matrix = NULL,
  lambda = "auto",
  warn_redundant = TRUE
)
```

## Arguments

- enabled:

  Logical. Whether to apply soft subspace projection.

- nuisance_mask:

  Path to NIfTI mask or logical vector indicating nuisance voxels.

- nuisance_matrix:

  Pre-computed nuisance timeseries matrix (alternative to mask).

- lambda:

  Ridge penalty: numeric, "auto", or "gcv".

- warn_redundant:

  Logical. Warn if baseline model contains nuisance terms.

## Value

A list of class "soft_subspace_options".

## Examples

``` r
# Using a mask file
opts <- soft_subspace_options(
  enabled = TRUE,
  nuisance_mask = "path/to/wm_csf_mask.nii.gz",
  lambda = "auto"
)
#> Error in soft_subspace_options(enabled = TRUE, nuisance_mask = "path/to/wm_csf_mask.nii.gz",     lambda = "auto"): could not find function "soft_subspace_options"

# Using pre-computed nuisance matrix
N <- matrix(rnorm(100 * 20), 100, 20)
opts <- soft_subspace_options(
  enabled = TRUE,
  nuisance_matrix = N,
  lambda = 0.5
)
#> Error in soft_subspace_options(enabled = TRUE, nuisance_matrix = N, lambda = 0.5): could not find function "soft_subspace_options"
```
