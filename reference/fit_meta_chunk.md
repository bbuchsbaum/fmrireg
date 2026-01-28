# Fit Meta-Analysis for a Chunk of Voxels

Fit Meta-Analysis for a Chunk of Voxels

## Usage

``` r
fit_meta_chunk(
  chunk_data,
  X,
  method,
  robust,
  weights,
  weights_custom,
  n_threads,
  combine = NULL,
  contrasts = NULL,
  return_cov = NULL
)
```

## Arguments

- chunk_data:

  List with beta and se matrices

- X:

  Design matrix

- method:

  Meta-analysis method

- robust:

  Robust estimation method

- weights:

  Weighting scheme

- weights_custom:

  Custom weights

- n_threads:

  Number of threads

## Value

List with meta-analysis results
