# Reshape Coefficient Data

This function reshapes coefficient data from wide to long format and
merges it with design information.

## Usage

``` r
reshape_coef(df, des, measure = "value")
```

## Arguments

- df:

  A data frame containing coefficient estimates.

- des:

  A data frame containing design information.

- measure:

  The name of the value column in the reshaped data. Default is
  `"value"`.

## Value

A data frame in long format with merged design information.
