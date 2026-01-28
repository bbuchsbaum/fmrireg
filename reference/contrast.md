# Contrast generic

Provide a contrast generic that dispatches on the first argument. Falls
back to fmridesign::contrast for non-fmri_meta classes.

## Usage

``` r
contrast(x, ...)
```

## Arguments

- x:

  object

- ...:

  passed to methods

## Value

A contrast object with computed contrast weights and statistics

## Examples

``` r
meta <- fmrireg:::.demo_fmri_meta()
#> Warning: fmri_meta.group_data_csv() is deprecated. Use fmri_meta() with a GDS-backed group_data() object.
contrast(meta, c("(Intercept)" = 1))
#> $estimate
#> [1] 0.25 0.15
#> 
#> $se
#> [1] 0.05773503 0.05773503
#> 
#> $z
#> [1] 4.330127 2.598076
#> 
#> $p
#>              [,1]
#> ROI1 1.490234e-05
#> ROI2 9.374768e-03
#> 
#> $weights
#> (Intercept) 
#>           1 
#> 
#> $name
#> [1] "(Intercept)"
#> 
#> $parent
#> fMRI Meta-Analysis Results
#> ==========================
#> 
#> Method: fe 
#> Robust: none 
#> Formula: ~1 
#> Subjects: 3 
#> ROIs analyzed: 2 
#> 
#> Heterogeneity:
#>   Mean tau^2: 0 
#>   Mean I^2: 0 %
#> 
#> attr(,"class")
#> [1] "fmri_meta_contrast"
```
