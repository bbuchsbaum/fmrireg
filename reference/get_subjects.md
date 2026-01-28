# Get Subject IDs

Get Subject IDs

## Usage

``` r
get_subjects(x)
```

## Arguments

- x:

  A group_data object

## Value

Character vector of unique subject IDs

## Examples

``` r
gd <- fmrireg:::.demo_group_data_csv()
get_subjects(gd)
#> [1] "s1" "s2" "s3"
```
