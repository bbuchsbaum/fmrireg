# Print Meta-Analysis Results

Print Meta-Analysis Results

## Usage

``` r
# S3 method for class 'fmri_meta'
print(x, ...)
```

## Arguments

- x:

  An fmri_meta object

- ...:

  Additional print arguments

## Value

Invisibly returns the input object x

## Examples

``` r
toy_meta <- structure(
  list(
    coefficients = matrix(0.2, nrow = 1),
    se = matrix(0.05, nrow = 1),
    method = "DL",
    robust = "none",
    formula = ~ condition,
    n_subjects = 12,
    n_rois = 1
  ),
  class = "fmri_meta"
)
print(toy_meta)
#> fMRI Meta-Analysis Results
#> ==========================
#> 
#> Method: DL 
#> Robust: none 
#> Formula: ~condition 
#> Subjects: 12 
#> ROIs analyzed: 1 
```
