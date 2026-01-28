# Apply Soft Projection to Data and Design

Applies the soft projection to both the data matrix Y and design matrix
X. Both must be projected to avoid bias in coefficient estimates.

## Usage

``` r
apply_soft_projection(proj, Y, X)
```

## Arguments

- proj:

  A soft_projection object from
  [`soft_projection()`](https://bbuchsbaum.github.io/fmrireg/reference/soft_projection.md).

- Y:

  Data matrix (time x voxels).

- X:

  Design matrix (time x predictors).

## Value

List with projected Y and X.

## Examples

``` r
set.seed(123)
N <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
X <- cbind(1, rnorm(100))

proj <- soft_projection(N, lambda = "auto")
#> Error in soft_projection(N, lambda = "auto"): could not find function "soft_projection"
cleaned <- apply_soft_projection(proj, Y, X)
#> Error in apply_soft_projection(proj, Y, X): could not find function "apply_soft_projection"
# Use cleaned$Y and cleaned$X for GLM fitting
```
