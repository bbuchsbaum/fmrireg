# Statistical fitting scope

Statistical fitting scope

## Usage

``` r
estimation_spec(
  scope = c("joint", "runwise_meta"),
  run_pooling = "fixed_effect"
)
```

## Arguments

- scope:

  Fit a single joint model or fit runs separately and combine them by
  fixed-effect meta-analysis.

- run_pooling:

  Runwise combination rule. Only fixed-effect pooling is currently
  supported.
