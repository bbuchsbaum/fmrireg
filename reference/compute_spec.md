# Mechanical execution specification

Mechanical execution specification

## Usage

``` r
compute_spec(
  voxel_chunks = 10L,
  backend = c("matrix", "reference"),
  parallel = c("none", "voxels", "chunks"),
  progress = FALSE
)
```

## Arguments

- voxel_chunks:

  Number of voxel partitions for joint fitting.

- backend:

  Matrix or formula/reference implementation.

- parallel:

  Parallelization dimension.

- progress:

  Display progress.
