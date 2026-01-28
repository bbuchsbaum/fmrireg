# Get Covariates

Get Covariates

## Usage

``` r
get_covariates(x)
```

## Arguments

- x:

  A group_data object

## Value

Data frame of covariates or NULL

## Examples

``` r
gd <- fmrireg:::.demo_group_data_csv()
get_covariates(gd)
#>    age
#> s1  30
#> s2  32
#> s3  34
```
