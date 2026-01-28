# Convert t-statistics to Effect Sizes

Helper function to convert t-statistics and df to beta and SE estimates

## Usage

``` r
t_to_beta_se(t, df, n = NULL)
```

## Arguments

- t:

  T-statistic values

- df:

  Degrees of freedom

- n:

  Sample size (optional, improves SE estimation)

## Value

List with beta and se estimates
