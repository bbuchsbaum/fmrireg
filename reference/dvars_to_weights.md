# Convert DVARS to Volume Weights

Transforms DVARS quality metrics into weights for weighted least squares
fitting. Volumes with high DVARS receive lower weights, implementing
"soft scrubbing" without hard censoring thresholds.

## Usage

``` r
dvars_to_weights(
  dvars,
  method = c("inverse_squared", "soft_threshold", "tukey"),
  threshold = 1.5,
  steepness = 2
)
```

## Arguments

- dvars:

  Numeric vector of DVARS values (from `compute_dvars`).

- method:

  Character. Weighting method: "inverse_squared" (default),
  "soft_threshold", or "tukey".

- threshold:

  Numeric. For "soft_threshold", the DVARS value above which weights
  decay. Default 1.5 (1.5x median if normalized).

- steepness:

  Numeric. For "soft_threshold", controls decay rate. Default 2.

## Value

Numeric vector of weights in the interval from 0 to 1, same length as
dvars.

## Examples

``` r
set.seed(123)
Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
Y[50, ] <- Y[50, ] + 5
dvars <- compute_dvars(Y)
#> Error in compute_dvars(Y): could not find function "compute_dvars"

# Compare different weighting methods
w_inv <- dvars_to_weights(dvars, method = "inverse_squared")
#> Error in dvars_to_weights(dvars, method = "inverse_squared"): could not find function "dvars_to_weights"
w_soft <- dvars_to_weights(dvars, method = "soft_threshold")
#> Error in dvars_to_weights(dvars, method = "soft_threshold"): could not find function "dvars_to_weights"
w_tukey <- dvars_to_weights(dvars, method = "tukey")
#> Error in dvars_to_weights(dvars, method = "tukey"): could not find function "dvars_to_weights"

# Check weight at the artifact volume
cat("Weights at volume 50:\n")
#> Weights at volume 50:
cat("  inverse_squared:", round(w_inv[50], 3), "\n")
#> Error: object 'w_inv' not found
cat("  soft_threshold:", round(w_soft[50], 3), "\n")
#> Error: object 'w_soft' not found
cat("  tukey:", round(w_tukey[50], 3), "\n")
#> Error: object 'w_tukey' not found
```
