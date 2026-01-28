# Get Available ROIs

Get Available ROIs

## Usage

``` r
get_rois(gd)
```

## Arguments

- gd:

  A group_data_csv object

## Value

Character vector of ROI names

## Examples

``` r
gd <- fmrireg:::.demo_group_data_csv()
get_rois(gd)
#> [1] "ROI1" "ROI2"
```
