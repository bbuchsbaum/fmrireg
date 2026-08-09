# Variance and reference-distribution specification

Variance and reference-distribution specification

## Usage

``` r
variance_spec(
  method = c("model", "hac", "sandwich"),
  max_lag = "auto",
  taper = c("none", "tukey", "parzen"),
  debias = TRUE,
  df = c("residual", "satterthwaite")
)
```

## Arguments

- method:

  Covariance estimator: model-based, HAC, or sandwich.

- max_lag:

  Maximum HAC lag, or `"auto"` for a data-dependent lag.

- taper:

  HAC lag-window taper.

- debias:

  Apply the estimator's finite-sample correction.

- df:

  Reference degrees-of-freedom method.
