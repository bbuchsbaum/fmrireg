# Create Soft Projection Operator

Computes the soft (ridge-regularized) projection matrix that removes
variance explainable by the nuisance subspace while avoiding
overfitting.

## Usage

``` r
soft_projection(N, lambda = "auto", Y = NULL)
```

## Arguments

- N:

  Nuisance matrix (time x nuisance_voxels). Typically extracted from
  white matter and CSF voxels. Can have many columns (thousands);
  computational cost depends on min(nrow, ncol) due to SVD.

- lambda:

  Ridge penalty controlling shrinkage strength:

  "auto"

  :   (Default) Uses median of squared singular values. Fast, stable,
      requires no tuning. Recommended for most use cases.

  "gcv"

  :   Generalized cross-validation. Optimizes prediction of Y from N
      with ridge penalty. Requires Y parameter. More principled but
      slower.

  numeric

  :   User-specified value. Larger = less aggressive removal.

- Y:

  Optional data matrix for GCV-based lambda selection. Required if
  `lambda = "gcv"`, ignored otherwise.

## Value

A list with class `"soft_projection"` containing:

- P_lambda:

  Function that applies the projection to a matrix

- lambda:

  The selected/specified lambda value

- method:

  Method used: "singular_value_heuristic", "gcv", or "user_specified"

- effective_df:

  Effective degrees of freedom removed (sum of shrinkage factors)

- n_nuisance:

  Number of nuisance columns in N

- n_timepoints:

  Number of timepoints

## Examples

``` r
# Create nuisance matrix (e.g., from WM/CSF voxels)
set.seed(123)
N <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
proj <- soft_projection(N, lambda = "auto")
#> Error in soft_projection(N, lambda = "auto"): could not find function "soft_projection"

# Apply to data and design
Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
Y_clean <- proj$P_lambda(Y)
#> Error in proj$P_lambda: object of type 'closure' is not subsettable

# Full workflow: project both data and design
X <- cbind(1, rnorm(100), rnorm(100))  # intercept + 2 predictors
cleaned <- apply_soft_projection(proj, Y, X)
#> Error in apply_soft_projection(proj, Y, X): could not find function "apply_soft_projection"
# Now fit GLM with cleaned$Y and cleaned$X

# Using GCV for lambda selection (data-driven)
proj_gcv <- soft_projection(N, lambda = "gcv", Y = Y)
#> Error in soft_projection(N, lambda = "gcv", Y = Y): could not find function "soft_projection"
print(proj_gcv)  # Shows selected lambda and effective df
#> Error: object 'proj_gcv' not found
```
