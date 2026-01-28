# Summary of Meta-Analysis Results

Summary of Meta-Analysis Results

## Usage

``` r
# S3 method for class 'fmri_meta'
summary(object, threshold = 0.05, ...)
```

## Arguments

- object:

  An fmri_meta object

- threshold:

  P-value threshold for significance (default: 0.05)

- ...:

  Additional summary arguments

## Value

A list containing summary statistics invisibly

## Examples

``` r
toy_meta <- structure(
  list(
    coefficients = matrix(c(0.3, -0.1), nrow = 1,
      dimnames = list(NULL, c("A", "B"))),
    se = matrix(c(0.05, 0.07), nrow = 1),
    method = "DL",
    robust = "none",
    formula = ~ condition,
    n_subjects = 12,
    n_rois = 1
  ),
  class = "fmri_meta"
)
summary(toy_meta, threshold = 0.1)
#> fMRI Meta-Analysis Summary
#> ==========================
#> 
#> fMRI Meta-Analysis Results
#> ==========================
#> 
#> Method: DL 
#> Robust: none 
#> Formula: ~condition 
#> Subjects: 12 
#> ROIs analyzed: 1 
#> 
#> Coefficients:
#>   A:
#>     Mean effect: 0.3 
#>     Mean SE: 0.05 
#>     Significant:1/1 (100%)
#>   B:
#>     Mean effect: -0.1 
#>     Mean SE: 0.07 
#>     Significant:0/1 (0%)
```
