# Get Available Contrasts

Get Available Contrasts

## Usage

``` r
get_contrasts(gd)
```

## Arguments

- gd:

  A group_data_csv object

## Value

Character vector of contrast names

## Examples

``` r
gd <- fmrireg:::.demo_group_data_csv()
get_contrasts(gd)
#> [1] "A_vs_B"
```
