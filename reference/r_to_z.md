# Convert Correlation to Fisher's Z

Transforms correlations to Fisher's Z scale for meta-analysis.

## Usage

``` r
r_to_z(r, n)
```

## Arguments

- r:

  Numeric vector or matrix of correlations

- n:

  Integer scalar; sample size

## Value

List with components:

- z:

  Numeric vector or matrix; Fisher's Z transformed correlations

- v:

  Numeric vector or matrix; sampling variances

## See also

[`fmri_meta`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md)

## Examples

``` r
r_to_z(r = 0.4, n = 30)
#> $z
#> [1] 0.4236489
#> 
#> $v
#> [1] 0.03703704
#> 
```
