# Extract term matrices from fmri_model

Extract design matrices for individual terms from an fmri_model object.

## Usage

``` r
# S3 method for class 'fmri_model'
term_matrices(x, blocknum = NULL, ...)
```

## Arguments

- x:

  An fmri_model object

- blocknum:

  Optional vector of block numbers to extract matrices for

- ...:

  Additional arguments (currently unused)

## Value

A list of matrices, one for each term in the model

## Examples

``` r
fm <- fmrireg:::.demo_fmri_model()
term_matrices(fm)
#> $condition
#>      condition_condition.A condition_condition.B
#> [1,]            0.06243065            0.00000000
#> [2,]            1.13751689            0.00000000
#> [3,]            1.77810440            0.00000000
#> [4,]            2.21043230            0.00000000
#> [5,]            0.00000000            0.04828652
#> [6,]            0.00000000            1.07265083
#> [7,]            0.00000000            1.81464049
#> [8,]            0.00000000            2.37677852
#> 
#> $drift
#>      base_bs1_block_1 base_bs2_block_1 base_bs3_block_1 base_bs1_block_2
#> [1,]        0.0000000        0.0000000       0.00000000        0.0000000
#> [2,]        0.4444444        0.2222222       0.03703704        0.0000000
#> [3,]        0.2222222        0.4444444       0.29629630        0.0000000
#> [4,]        0.0000000        0.0000000       1.00000000        0.0000000
#> [5,]        0.0000000        0.0000000       0.00000000        0.0000000
#> [6,]        0.0000000        0.0000000       0.00000000        0.4444444
#> [7,]        0.0000000        0.0000000       0.00000000        0.2222222
#> [8,]        0.0000000        0.0000000       0.00000000        0.0000000
#>      base_bs2_block_2 base_bs3_block_2
#> [1,]        0.0000000       0.00000000
#> [2,]        0.0000000       0.00000000
#> [3,]        0.0000000       0.00000000
#> [4,]        0.0000000       0.00000000
#> [5,]        0.0000000       0.00000000
#> [6,]        0.2222222       0.03703704
#> [7,]        0.4444444       0.29629630
#> [8,]        0.0000000       1.00000000
#> 
#> $block
#>      constant_1 constant_2
#> [1,]          1          0
#> [2,]          1          0
#> [3,]          1          0
#> [4,]          1          0
#> [5,]          0          1
#> [6,]          0          1
#> [7,]          0          1
#> [8,]          0          1
#> 
#> attr(,"event_term_indices")
#> [1] 1 2
#> attr(,"baseline_term_indices")
#> [1]  3  4  5  6  7  8  9 10
#> attr(,"blocknum")
#> [1] 1 2
#> attr(,"varnames")
#>  [1] "condition_condition.A" "condition_condition.B" "base_bs1_block_1"     
#>  [4] "base_bs2_block_1"      "base_bs3_block_1"      "base_bs1_block_2"     
#>  [7] "base_bs2_block_2"      "base_bs3_block_2"      "constant_1"           
#> [10] "constant_2"           
```
