# Build sparse landmark weight matrix W (V x L) using k-NN heat-kernel

Build sparse landmark weight matrix W (V x L) using k-NN heat-kernel

## Usage

``` r
build_landmark_weights(coords, lcoords, k = 16L, h = NULL)
```

## Arguments

- coords:

  V x d matrix of voxel coordinates (d = 2 or 3)

- lcoords:

  L x d matrix of landmark coordinates (subset of coords or arbitrary
  within space)

- k:

  integer; number of nearest landmarks used per voxel (default 16)

- h:

  optional bandwidth; if NULL, uses median of k-NN distances

## Value

dgCMatrix W with rows summing to 1 (each voxel is a convex combination
of k landmarks)
