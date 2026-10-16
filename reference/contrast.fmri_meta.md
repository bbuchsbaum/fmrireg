# Apply Contrast to Meta-Analysis Results

The contrast standard error is \\\sqrt{c' Var(\hat\beta) c}\\, using the
coefficient covariance that
[`fmri_meta()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md)
stores with every fit. For objects without a stored covariance only
single-coefficient contrasts are exact; other contrasts fall back to a
diagonal approximation with a warning.

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
#> Warning: Coefficient covariance is not available for this fmri_meta object; the contrast standard error ignores covariances between coefficients and is approximate.
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
