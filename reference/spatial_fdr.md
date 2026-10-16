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

  Integer, character or factor vector of group IDs for each feature
  (e.g., parcel IDs, block IDs). Must have same length as z or p. Groups
  are numbered 1..G by sorted unique value, or by factor level (unused
  levels included) when `group` is a factor. This numbering indexes
  `neighbors`, `pi0_raw`, `pi0_smooth` and `m_g`. Features with `NA`
  group are excluded from \\\pi_0\\ estimation and receive weight 1
  before normalisation.

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

  Optional list of length G (number of groups) where element g is an
  integer vector of 1-based neighbor IDs of group g, both in the group
  numbering described under `group`. Used for spatial smoothing of
  \\\pi_0\\.

- min_pi0:

  Numeric scalar; lower bound for \\\pi_0\\ to stabilize weights
  (default: 0.05). Prevents infinite weights.

- empirical_null:

  Logical; if TRUE, estimate null distribution parameters (\\\mu_0\\,
  \\\sigma_0\\) from the z-values (default: TRUE). See Details.

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

  Integer vector of group IDs (1..G) in the numbering described under
  `group`

- G:

  Integer scalar; number of groups

- m_g:

  Integer vector; number of features per group

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

Weights are \\1/\pi_0\\ of each feature's group, normalised to sum to
the number of finite p-values, and the weighted BH procedure of
Genovese, Roeder & Wasserman (2006) is applied: reject \\p_i \le w_i
\alpha k / n\\.

**Empirical null.** When `z` is supplied with `empirical_null = TRUE`,
\\\mu_0\\ and \\\sigma_0\\ are estimated by Efron's (2004) central
matching: a quadratic is fitted by Poisson regression to the log
histogram counts of z within \\\pm 1.5\\ robust SDs of the median.
Because a minority of non-null features only changes the intercept of
the log density near the centre and not its curvature, this estimate is
consistent under \\N(0,1)\\ and stays close to it when, for example,
10\\ of features are true signals. With fewer than 500 finite z-values
(or if the fit fails) the median and normal-consistent MAD of all finite
z-values are used instead; with fewer than 10 the theoretical null is
used.

The method is particularly effective for:

- Voxelwise analyses with spatial clustering of signal

- Parcel-based analyses with anatomical or functional grouping

- Any scenario where hypotheses can be grouped spatially or functionally

## References

- Benjamini & Hochberg (1995). Controlling the false discovery rate.

- Storey (2002). A direct approach to false discovery rates.

- Hu et al. (2010). False discovery rate control with groups (SABHA).

- Genovese, Roeder & Wasserman (2006). False discovery control with
  p-value weighting.

- Efron (2004). Large-scale simultaneous hypothesis testing: the choice
  of a null hypothesis.

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
#> Discoveries: 8 ( 0.8 %)
#> Threshold: 4e-04 
#> Empirical null: mu = 0.205 , sigma = 1.153 
#> 
#> Pi0 summary:
#>   Range: 0.32 - 1 
#>   Mean: 0.868 
#>   Median: 1 
#> 
#> Group-level discoveries:
#>   Groups with discoveries: 2 / 10 
#>   Max proportion in group: 0.04 
#>   Mean proportion in groups with signal: 0.04 

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
#> Discoveries: 9 ( 0.9 %)
#> Threshold: 0.00045 
#> Empirical null: mu = 0.205 , sigma = 1.153 
#> 
#> Pi0 summary:
#>   Range: 0.34 - 1 
#>   Mean: 0.867 
#>   Median: 1 
```
