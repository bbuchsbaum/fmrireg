# Mixed Model Solver using Rcpp and roptim

This function solves a mixed model using Rcpp and roptim for
optimization. It estimates variance components in a mixed model,
potentially speeding up computations compared to the pure R
implementation.

## Usage

``` r
mixed_solve_cpp(
  y,
  Z = NULL,
  K = NULL,
  X = NULL,
  method = "REML",
  bounds = c(1e-09, 1e+09),
  SE = FALSE,
  return_Hinv = FALSE
)
```

## Arguments

- y:

  Response vector.

- Z:

  Design matrix for random effects (default: identity matrix of size n).

- K:

  Kinship matrix (default: NULL).

- X:

  Design matrix for fixed effects (default: vector of ones).

- method:

  Optimization method, either "REML" or "ML" (default: "REML").

- bounds:

  Bounds for the optimizer (default: c(1e-9, 1e9)).

- SE:

  Logical, whether to return standard errors (default: FALSE).

- return_Hinv:

  Logical, whether to return the inverse of H (default: FALSE).

## Value

A list containing:

- Vu:

  Estimated variance component for random effects.

- Ve:

  Estimated variance component for residuals.

- beta:

  Estimated fixed effects coefficients.

- u:

  Estimated random effects coefficients.

- LL:

  Log-likelihood of the model.

- beta.SE:

  Standard errors of fixed effects coefficients (if SE = TRUE).

- u.SE:

  Standard errors of random effects coefficients (if SE = TRUE).

- Hinv:

  Inverse of H (if return_Hinv = TRUE).

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage with random data
set.seed(123)
n <- 100
y <- rnorm(n)
Z <- matrix(rnorm(n * 5), n, 5)
K <- diag(5)
X <- matrix(1, n, 1)
result <- mixed_solve_cpp(y, Z, K, X)
} # }
```
