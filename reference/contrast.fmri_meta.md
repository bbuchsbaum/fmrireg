# Apply Contrast to Meta-Analysis Results

Note: Uses exact standard errors when covariance is available
(return_cov="tri") or for ROI CSV fits with non-robust estimation;
otherwise uses a diagonal variance approximation.

## Usage

``` r
# S3 method for class 'fmri_meta'
contrast(x, contrast, ...)
```

## Arguments

- x:

  An fmri_meta object

- contrast:

  Contrast specification. Can be:

  - A numeric vector of contrast weights

  - A formula (e.g., ~ groupold - groupyoung)

  - A named vector (e.g., c("groupold" = 1, "groupyoung" = -1))

- ...:

  Additional arguments

## Value

An fmri_meta_contrast object with contrast results

## Examples

``` r
toy_meta <- structure(
  list(
    coefficients = matrix(c(0.3, 0.1), nrow = 1,
      dimnames = list(NULL, c("A", "B"))),
    se = matrix(c(0.05, 0.06), nrow = 1),
    robust = "none"
  ),
  class = "fmri_meta"
)
contrast(toy_meta, c(1, -1))
#> $estimate
#> [1] 0.2
#> 
#> $se
#> [1] 0.0781025
#> 
#> $z
#> [1] 2.560738
#> 
#> $p
#>            [,1]
#> [1,] 0.01044502
#> 
#> $weights
#> [1]  1 -1
#> 
#> $name
#> [1] "custom"
#> 
#> $parent
#> fMRI Meta-Analysis Results
#> ==========================
#> 
#> Method: 
#> Robust: none 
#> Formula: NULL 
#> Subjects: 
#> 
#> attr(,"class")
#> [1] "fmri_meta_contrast"
```
