# Compute Effective Sample Size for Meta-Analysis

Computes the effective sample size based on the heterogeneity estimate.
This is useful for understanding the impact of between-study
heterogeneity.

## Usage

``` r
meta_effective_n(v, tau2)
```

## Arguments

- v:

  Numeric vector of within-study variances

- tau2:

  Numeric scalar; between-study variance (tau-squared)

## Value

Numeric scalar; effective sample size

## See also

[`fmri_meta`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md)

## Examples

``` r
meta_effective_n(v = rep(0.05, 3), tau2 = 0.01)
#> [1] 3
```
