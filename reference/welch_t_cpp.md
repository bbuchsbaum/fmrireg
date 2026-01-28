# Welch two-sample t-test across features

Welch two-sample t-test across features

## Usage

``` r
welch_t_cpp(Y, g_in)
```

## Arguments

- Y:

  S x P matrix

- g_in:

  Length S vector of group indicators (1/2 or 0/1)

## Value

List with muA, muB, t, df (Welch), nA, nB
