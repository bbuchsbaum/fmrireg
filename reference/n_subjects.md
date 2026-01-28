# Extract Number of Subjects

Extract Number of Subjects

## Usage

``` r
n_subjects(x)
```

## Arguments

- x:

  A group_data object

## Value

Integer number of unique subjects

## Examples

``` r
gd <- fmrireg:::.demo_group_data_csv()
n_subjects(gd)
#> [1] 3
```
