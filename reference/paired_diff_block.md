# Helper Functions for fmri_ttest

Support functions for paired differences, sign flipping, and wrapper
functions for OLS and meta-analysis with voxelwise covariates. Compute
Paired Within-Subject Differences

## Usage

``` r
paired_diff_block(blkA, blkB, rho = 0)
```

## Arguments

- blkA:

  First group_data block

- blkB:

  Second group_data block

- rho:

  Optional within-subject correlation between A and B. Can be:

  - Scalar: Same correlation for all subjects and features

  - Vector of length S: Per-subject correlations

  - Vector of length P: Per-feature correlations

  - Matrix (S x P): Subject-feature specific correlations

  Default is 0 (independence).

## Value

A new block with Y = Y_A - Y_B and propagated variance if available

## Details

Creates within-subject differences (A - B) for paired t-tests from two
blocks with identical subjects and features.
