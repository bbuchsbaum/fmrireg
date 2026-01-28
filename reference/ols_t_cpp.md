# OLS t-test / ANCOVA across features

OLS t-test / ANCOVA across features

## Usage

``` r
ols_t_cpp(Y, X)
```

## Arguments

- Y:

  S x P matrix (subjects x features)

- X:

  S x K design matrix with intercept if desired

## Value

List with beta (K x P), se (K x P), t (K x P), df (scalar), ok (P)
