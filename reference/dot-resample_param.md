# Resample parameter vector with specified distribution

Helper for drawing per-event amplitudes/durations around base values.

## Usage

``` r
.resample_param(
  base,
  sd,
  dist = c("lognormal", "gamma", "gaussian"),
  allow_negative = FALSE
)
```

## Arguments

- base:

  Numeric vector of base values to jitter

- sd:

  Numeric standard deviation for the sampling distribution

- dist:

  Distribution name: "lognormal", "gamma", or "gaussian"

- allow_negative:

  Logical; if TRUE, allow negative draws (only used for gaussian)

## Value

A numeric vector of resampled values with the same length as base
