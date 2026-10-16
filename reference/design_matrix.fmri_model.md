# Design Matrix for fMRI Models

Extract the combined design matrix from an fMRI model containing both
event and baseline terms.

## Usage

``` r
# S3 method for class 'fmri_model'
design_matrix(x, blockid = NULL, ...)
```

## Arguments

- x:

  An fmri_model object

- blockid:

  Optional numeric vector specifying which blocks/runs to include

- ...:

  Additional arguments (not used)

## Value

A tibble containing the combined design matrix with event and baseline
terms

## Examples

``` r
fm <- fmrireg:::.demo_fmri_model()
head(design_matrix(fm))
#> # A tibble: 6 × 10
#>   condition_condition.A condition_condition.B base_bs1_block_1 base_bs2_block_1
#>                   <dbl>                 <dbl>            <dbl>            <dbl>
#> 1               0.00342               0                  0                0    
#> 2               0.101                 0                  0.444            0.222
#> 3               0.178                 0                  0.222            0.444
#> 4               0.228                 0                  0                0    
#> 5               0                     0.00342            0                0    
#> 6               0                     0.101              0                0    
#> # ℹ 6 more variables: base_bs3_block_1 <dbl>, base_bs1_block_2 <dbl>,
#> #   base_bs2_block_2 <dbl>, base_bs3_block_2 <dbl>, constant_1 <dbl>,
#> #   constant_2 <dbl>
```
