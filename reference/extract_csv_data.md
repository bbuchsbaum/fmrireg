# Extract Data for Meta-Analysis from CSV

Extract Data for Meta-Analysis from CSV

## Usage

``` r
extract_csv_data(gd, roi = NULL, contrast = NULL)
```

## Arguments

- gd:

  A group_data_csv object

- roi:

  Optional ROI name to extract

- contrast:

  Optional contrast name to extract

## Value

List with effect sizes and variances

## Examples

``` r
gd <- fmrireg:::.demo_group_data_csv()
extract_csv_data(gd, roi = "ROI1")
#> $beta
#> [1] 0.20 0.30 0.25
#> 
#> $se
#> [1] 0.1 0.1 0.1
#> 
#> $var
#> [1] 0.01 0.01 0.01
#> 
#> $subjects
#> [1] "s1" "s2" "s3"
#> 
```
