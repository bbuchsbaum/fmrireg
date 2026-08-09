# Volume-weight specification

Volume-weight specification

## Usage

``` r
weights_spec(
  method = c("none", "inverse_squared", "soft_threshold", "tukey"),
  threshold = 1.5,
  values = NULL
)
```

## Arguments

- method:

  Volume-weight construction rule.

- threshold:

  Positive threshold or tuning constant for the selected rule.

- values:

  Optional explicit finite volume weights.
