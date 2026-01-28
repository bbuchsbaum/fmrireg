# Soft Subspace Projection for Nuisance Removal

Ridge-regularized projection to remove nuisance variance without
explicitly choosing components or adding confound regressors. This is
"CompCor without choosing components."

## Conceptual Overview

Traditional CompCor extracts K principal components from white
matter/CSF and regresses them out. This requires choosing K, which is
arbitrary.

Soft subspace projection instead treats the entire WM/CSF timeseries as
a nuisance basis and removes it with shrinkage. Each nuisance direction
is partially removed proportional to its variance, avoiding the hard
keep/delete decision of PCA truncation.

## The Projection Operator

Given nuisance matrix N (time x nuisance_voxels), the soft projection
is: \$\$P\_\lambda = I - N (N^T N + \lambda I)^{-1} N^T\$\$

Applied to data Y and design X:

- `Y_clean = P_lambda %*% Y`

- `X_clean = P_lambda %*% X` (important to avoid bias)

## Lambda Selection

The shrinkage parameter lambda controls aggressiveness:

- Small lambda: aggressive removal (risk removing signal)

- Large lambda: gentle removal (risk leaving nuisance)

Three selection methods are available:

- "auto":

  Singular value heuristic: `lambda = median(d^2)` where d are singular
  values of N. Fast, stable, no tuning. Components with variance below
  median are heavily shrunk; above median lightly shrunk.

- "gcv":

  Generalized cross-validation minimizes `RSS/(1-df/n)^2` for ridge
  regression of Y on N. Finds lambda giving best leave-one-out
  prediction without actually doing LOO. Requires Y. Computational cost
  is O(k) per evaluation where k = min(n, p).

- numeric:

  User-specified lambda value.

## Typical Usage

The soft subspace projection workflow has three steps:

1.  Extract nuisance timeseries from WM/CSF mask (or provide
    pre-computed)

2.  Create the soft projection operator

3.  Apply to both data Y and design matrix X before GLM fitting

For use within
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md),
see the convenience parameter `nuisance_projection` or the more detailed
`soft_subspace_options`.

## When to Use

Soft subspace projection is most beneficial when:

- You have physiological noise from WM/CSF that isn't captured by motion
  parameters

- Traditional CompCor requires arbitrary component selection

- You want to avoid adding many confound regressors to the design

Consider alternatives when:

- Motion parameters alone sufficiently control artifacts

- You prefer explicit confound regressors for interpretability

- Data has very few timepoints relative to nuisance dimensions
