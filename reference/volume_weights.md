# Compute Volume Quality Weights from Data

Convenience function that computes DVARS and converts to weights in one
step. This is the main user-facing function for volume quality
weighting.

## Usage

``` r
volume_weights(
  Y,
  method = "inverse_squared",
  threshold = 1.5,
  return_dvars = FALSE
)
```

## Arguments

- Y:

  Numeric matrix of fMRI data (time x voxels).

- method:

  Weighting method passed to `dvars_to_weights`.

- threshold:

  Threshold passed to `dvars_to_weights`.

- return_dvars:

  Logical. If TRUE, return both weights and DVARS.

## Value

If return_dvars is FALSE, a numeric vector of weights. If TRUE, a list
with components "weights" and "dvars".

## Examples

``` r
set.seed(123)
Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
Y[50, ] <- Y[50, ] + 5  # Add artifact

# One-step computation of weights
result <- volume_weights(Y, return_dvars = TRUE)
cat("DVARS at artifact:", round(result$dvars[50], 2), "\n")
#> DVARS at artifact: 3.38 
cat("Weight at artifact:", round(result$weights[50], 3), "\n")
#> Weight at artifact: 0.163 

# With Tukey method for more aggressive downweighting
w_tukey <- volume_weights(Y, method = "tukey")
```
