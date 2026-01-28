# Convert T-statistics to Effect Sizes and Variances

Converts t-statistics and degrees of freedom to standardized mean
differences (Cohen's d) and their sampling variances for meta-analysis.

## Usage

``` r
t_to_d(t, df, n = NULL)
```

## Arguments

- t:

  Numeric vector or matrix of t-statistics

- df:

  Numeric scalar or vector; degrees of freedom (matching t)

- n:

  Numeric scalar; sample size per group (for two-sample t-tests)

## Value

List with components:

- d:

  Numeric vector or matrix; standardized mean differences

- v:

  Numeric vector or matrix; sampling variances

## See also

[`fmri_meta`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md)

## Examples

``` r
t_to_d(t = 2, df = 18)
#> $d
#> [1] 0.4588315
#> 
#> $v
#> [1] 0.05817175
#> 
t_to_d(t = 2, df = 18, n = 20)
#> $d
#> [1] 0.942809
#> 
#> $v
#> [1] 0.2246914
#> 
```
