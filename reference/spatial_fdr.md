# Spatially-Aware Multiple Comparisons Correction

Performs spatially-aware FDR control using structure-adaptive weighted
Benjamini-Hochberg (SABHA-style) procedure. This method leverages
spatial structure in the data to increase power while controlling the
false discovery rate.

## Usage

``` r
spatial_fdr(
  z = NULL,
  p = NULL,
  group = NULL,
  alpha = 0.05,
  tau = 0.5,
  lambda = 1,
  neighbors = NULL,
  min_pi0 = 0.05,
  empirical_null = TRUE,
  verbose = FALSE
)
```

## Arguments

- z:

  Numeric vector of Z-values (one per feature). Provide either z or p,
  not both.

- p:

  Numeric vector of p-values (two-sided). Provide either z or p, not
  both.

- group:

  Integer or factor vector of group IDs for each feature (e.g., parcel
  IDs, block IDs). Must have same length as z or p.

- alpha:

  Numeric scalar; FDR level to control (default: 0.05)

- tau:

  Numeric scalar; Storey threshold in (0,1) for \\\pi_0\\ estimation
  (default: 0.5). Higher values are more conservative.

- lambda:

  Numeric scalar; smoothing strength for groupwise \\\pi_0\\ across
  neighbors (default: 1.0). Set to 0 for no smoothing, higher values for
  more smoothing.

- neighbors:

  Optional list of length G (number of groups) where each element is an
  integer vector of 1-based neighbor IDs. Used for spatial smoothing of
  \\\pi_0\\.

- min_pi0:

  Numeric scalar; lower bound for \\\pi_0\\ to stabilize weights
  (default: 0.05). Prevents infinite weights.

- empirical_null:

  Logical; if TRUE, estimate null distribution parameters (\\\mu_0\\,
  \\\sigma_0\\) from central z-values using robust estimators (default:
  TRUE).

- verbose:

  Logical; print progress messages (default: FALSE)

## Value

Object of class "spatial_fdr_result" containing:

- reject:

  Logical vector indicating rejected hypotheses (discoveries)

- q:

  Numeric vector of FDR-adjusted p-values (q-values)

- p:

  Numeric vector of two-sided p-values used for testing

- weights:

  Numeric vector of normalized weights used in weighted BH

- pi0_raw:

  Numeric vector of raw \\\pi_0\\ estimates per group

- pi0_smooth:

  Numeric vector of smoothed \\\pi_0\\ estimates per group

- threshold:

  Numeric scalar; BH threshold used for rejection

- k:

  Integer scalar; number of rejections

- mu0:

  Numeric scalar; estimated null mean (if empirical_null = TRUE)

- sigma0:

  Numeric scalar; estimated null SD (if empirical_null = TRUE)

- groups:

  Integer vector of compressed group IDs (1..G)

- group:

  Factor or integer vector of original group IDs

- G:

  Integer scalar; number of groups

- alpha:

  Numeric scalar; FDR level used

- coef_name:

  Character scalar; name of coefficient (for S3 method)

## Details

This function implements a spatially-aware multiple testing procedure
that:

1.  Estimates the proportion of null hypotheses (pi0) within spatial
    groups

2.  Optionally smooths these estimates across neighboring groups

3.  Uses the pi0 estimates to weight the Benjamini-Hochberg procedure

4.  Provides more power in regions with true signal while maintaining
    FDR control

The method is particularly effective for:

- Voxelwise analyses with spatial clustering of signal

- Parcel-based analyses with anatomical or functional grouping

- Any scenario where hypotheses can be grouped spatially or functionally

## References

- Benjamini & Hochberg (1995). Controlling the false discovery rate.

- Storey (2002). A direct approach to false discovery rates.

- Hu et al. (2010). False discovery rate control with groups (SABHA).

## See also

[`create_3d_blocks`](https://bbuchsbaum.github.io/fmrireg/reference/create_3d_blocks.md)

## Examples

``` r
# Simple synthetic example
set.seed(123)
n <- 1000
z_scores <- c(rnorm(800), rnorm(200, mean = 2))  # 200 true signals
group_ids <- rep(1:10, each = 100)  # 10 groups of 100 features each

# Basic usage without spatial smoothing
result <- spatial_fdr(z = z_scores, group = group_ids, alpha = 0.05)
summary(result)
#> Spatial FDR Results
#> ===================
#> Features: 1000 
#> Groups: 10 
#> FDR level: 0.05 
#> Discoveries: 180 ( 18 %)
#> Threshold: 0.009 
#> Empirical null: mu = 0.079 , sigma = 0.673 
#> 
#> Pi0 summary:
#>   Range: 0.08 - 0.8 
#>   Mean: 0.6 
#>   Median: 0.69 
#> 
#> Group-level discoveries:
#>   Groups with discoveries: 10 / 10 
#>   Max proportion in group: 0.69 
#>   Mean proportion in groups with signal: 0.18 

# Create simple neighbor structure (each group neighbors with adjacent groups)
neighbors <- lapply(1:10, function(i) {
  c(if(i > 1) i-1, if(i < 10) i+1)
})

# With spatial smoothing
result_smooth <- spatial_fdr(z = z_scores, group = group_ids, 
                            neighbors = neighbors, lambda = 1.0)
print(result_smooth)
#> Spatial FDR Results
#> ===================
#> Features: 1000 
#> Groups: 10 
#> FDR level: 0.05 
#> Discoveries: 185 ( 18.5 %)
#> Threshold: 0.00925 
#> Empirical null: mu = 0.079 , sigma = 0.673 
#> 
#> Pi0 summary:
#>   Range: 0.12 - 0.76 
#>   Mean: 0.597 
#>   Median: 0.688 
```
