#' Fractional Ridge Regression
#'
#' Performs fractional ridge regression, approximating regularization parameters to match desired
#' fractions of the Ordinary Least Squares (OLS) parameter vector length.
#'
#' @param X A numeric matrix of shape \code{(n, p)} representing the design matrix,
#'   where \code{n} is the number of observations and \code{p} is the number of predictors.
#' @param y A numeric vector or matrix of shape \code{(n, b)} representing the response variable(s),
#'   where \code{b} is the number of targets. If \code{y} is a vector, it is treated as a single target.
#' @param fracs A numeric vector specifying the desired fractions of the OLS parameter vector length.
#'   Must be sorted in increasing order and be between 0 and 1. Default is \code{seq(0.1, 1.0, by = 0.1)}.
#' @param tol A numeric value specifying the tolerance below which singular values are considered zero.
#'   Must be positive. Default is \code{1e-10}.
#'
#' @return An object of class \code{fracridge} containing:
#' \describe{
#'   \item{\code{coef}}{A numeric array of shape \code{(p, length(fracs), b)} containing the estimated coefficients
#'     for each fraction and target. If \code{y} is a vector, the third dimension is dropped.}
#'   \item{\code{alpha}}{A numeric matrix of shape \code{(length(fracs), b)} containing the regularization parameters
#'     associated with each solution.}
#'   \item{\code{fracs}}{The fractions used in the regression.}
#'   \item{\code{rank}}{The effective rank of the design matrix \code{X}.}
#'   \item{\code{singular_values}}{The singular values of \code{X}.}
#' }
#'
#' @details
#' Fractional ridge regression finds ridge regression solutions where the parameter vector lengths
#' equal specified fractions of the OLS solution length. This is achieved by adjusting the 
#' regularization parameter \code{alpha} for each target fraction.
#'
#' The function uses Singular Value Decomposition (SVD) for efficient computation. Singular values
#' below \code{tol} are treated as zero, and a warning is issued if rank deficiency is detected.
#'
#' The regularization path is computed using a grid of alpha values on a log scale, with
#' interpolation used to find the exact alpha values that achieve the target fractions.
#'
#' @examples
#' \dontrun{
#' # Generate random data
#' set.seed(0)
#' X <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' y <- rnorm(100)
#'
#' # Perform fractional ridge regression
#' result <- fracridge(X, y, fracs = seq(0.3, 0.9, by = 0.1))
#'
#' # Inspect the coefficients and regularization parameters
#' print(result$coef)
#' print(result$alpha)
#'
#' # Compare with OLS solution length
#' ols_coef <- solve(crossprod(X), crossprod(X, y))
#' ols_norm <- sqrt(sum(ols_coef^2))
#' coef_norm <- sqrt(colSums(result$coef^2))
#' print(coef_norm / ols_norm)  # Should approximately equal fracs
#' }
#'
#' @export
fracridge <- function(X, y, fracs = seq(0.1, 1.0, by = 0.1), tol = 1e-10) {
  # Input validation
  if (!is.matrix(X)) {
    stop("X must be a numeric matrix")
  }
  if (!(is.numeric(y) || is.matrix(y))) {
    stop("y must be a numeric vector or matrix")
  }
  if (!is.numeric(fracs) || any(fracs <= 0) || any(fracs > 1) || any(diff(fracs) < 0)) {
    stop("fracs must be a sorted numeric vector with values between 0 and 1")
  }
  if (!is.numeric(tol) || tol <= 0) {
    stop("tol must be a positive number")
  }
  
  # Ensure y is a matrix
  if (is.vector(y)) {
    y <- as.matrix(y)
    single_target <- TRUE
  } else {
    y <- as.matrix(y)
    single_target <- FALSE
  }
  
  n <- nrow(X)
  p <- ncol(X)
  b <- ncol(y)
  f <- length(fracs)
  
  # Perform SVD
  svd_result <- svd(X)
  U <- svd_result$u
  S <- svd_result$d
  V <- svd_result$v
  
  # Check rank
  rank <- sum(S > tol)
  if (rank < p) {
    warning(sprintf("X is rank deficient. Effective rank %d < %d columns", rank, p))
  }
  
  # Compute OLS coefficients in rotated space
  S_inv <- ifelse(S > tol, 1/S, 0)
  UTy <- crossprod(U, y)
  beta_rotated <- S_inv * UTy
  
  # Compute OLS coefficient norms
  ols_norm <- sqrt(colSums(beta_rotated^2))
  
  # Create alpha grid for interpolation (log-spaced)
  alpha_grid <- 10^seq(-3, 4, length.out = 100)
  alpha_grid <- c(0, alpha_grid)  # Include unregularized solution
  
  # Pre-allocate arrays
  coef_array <- array(0, dim = c(p, f, b))
  alpha_matrix <- matrix(0, nrow = f, ncol = b)
  
  # Compute solutions for each target
  for (i in seq_len(b)) {
    # Compute coefficient norms for each alpha
    beta_norms <- sapply(alpha_grid, function(alpha) {
      sqrt(sum((S^2 / (S^2 + alpha) * beta_rotated[, i])^2))
    })
    beta_norms <- beta_norms / beta_norms[1]  # Normalize
    
    # Interpolate to find target alphas
    target_alphas <- approx(x = rev(beta_norms), y = rev(log1p(alpha_grid)),
                           xout = fracs, rule = 2)$y
    target_alphas <- exp(target_alphas) - 1
    
    # Store alphas
    alpha_matrix[, i] <- target_alphas
    
    # Compute and store coefficients
    for (j in seq_len(f)) {
      scale <- S^2 / (S^2 + target_alphas[j])
      coef_array[, j, i] <- V %*% (scale * beta_rotated[, i])
    }
  }
  
  # Format output
  if (single_target) {
    coef_array <- drop(coef_array)
    alpha_matrix <- drop(alpha_matrix)
  }
  
  dimnames(coef_array) <- list(
    colnames(X),
    sprintf("frac_%.2f", fracs),
    if (single_target) NULL else sprintf("target_%d", seq_len(b))
  )
  
  structure(
    list(
      coef = coef_array,
      alpha = alpha_matrix,
      fracs = fracs,
      rank = rank,
      singular_values = S
    ),
    class = "fracridge"
  )
}

#' Predict Method for fracridge Objects
#'
#' Generates predictions from a \code{fracridge} model.
#'
#' @param object An object of class \code{fracridge}.
#' @param newdata A numeric matrix of new data with the same number of columns as the training data.
#' @param fracs A numeric vector specifying which fractions to use for prediction.
#'   Must be a subset of the fractions used in the training. Default is \code{object$fracs}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric matrix of predictions with rows corresponding to observations in \code{newdata}
#'   and columns corresponding to the specified fractions and targets.
#'
#' @examples
#' \dontrun{
#' # Generate random data
#' set.seed(0)
#' X <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' y <- rnorm(100)
#'
#' # Perform fractional ridge regression
#' result <- fracridge(X, y, fracs = 0.3)
#'
#' # Generate new data
#' X_new <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)
#'
#' # Make predictions
#' preds <- predict(result, X_new)
#' }
#'
#' @export
predict.fracridge <- function(object, newdata, fracs = object$fracs, ...) {
  if (!inherits(object, "fracridge")) {
    stop("object must be of class 'fracridge'.")
  }
  
  if (!is.matrix(newdata)) {
    newdata <- as.matrix(newdata)
  }
  
  p <- ncol(newdata)
  if (p != dim(object$coef)[1]) {
    stop("Number of columns in newdata must match the number of predictors in the model.")
  }
  
  # Find indices of requested fractions
  frac_indices <- match(fracs, object$fracs)
  if (any(is.na(frac_indices))) {
    stop("All requested fracs must be present in the model.")
  }
  
  # Extract relevant coefficients
  coef_subset <- object$coef[, frac_indices, , drop = FALSE]
  
  # Compute predictions
  preds <- newdata %*% coef_subset
  
  # If single target, drop the third dimension
  if (length(dim(preds)) == 2 && dim(coef_subset)[3] == 1) {
    preds <- matrix(preds, ncol = length(fracs))
    colnames(preds) <- paste0("frac_", fracs)
  } else {
    # Flatten the last two dimensions
    preds <- matrix(preds, ncol = length(fracs) * dim(coef_subset)[3])
    colnames(preds) <- as.vector(outer(paste0("frac_", fracs), 
                                       if (dim(coef_subset)[3] > 1) paste0("_target_", 1:dim(coef_subset)[3]) else "", 
                                       paste0))
  }
  
  return(preds)
}

#' Print Method for fracridge Objects
#'
#' @param x An object of class \code{fracridge}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \dontrun{
#' # Generate random data
#' set.seed(0)
#' X <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' y <- rnorm(100)
#'
#' # Perform fractional ridge regression
#' result <- fracridge(X, y, fracs = 0.3)
#'
#' # Print the result
#' print(result)
#' }
#'
#' @export
print.fracridge <- function(x, ...) {
  cat("Fractional Ridge Regression Model\n")
  cat("Fractions used:", paste(x$fracs, collapse = ", "), "\n")
  cat("Regularization parameters:\n")
  print(x$alpha)
  cat("Effective rank:", x$rank, "\n")
  cat("Singular values:\n")
  print(x$singular_values)
  invisible(x)
}

#' Fractional Ridge Regression with Cross-Validation
#'
#' Performs fractional ridge regression with cross-validation to select the optimal fraction.
#'
#' @param X A numeric matrix of shape \code{(n, p)} representing the design matrix,
#'   where \code{n} is the number of observations and \code{p} is the number of predictors.
#' @param y A numeric vector or matrix of shape \code{(n, b)} representing the response variable(s),
#'   where \code{b} is the number of targets. If \code{y} is a vector, it is treated as a single target.
#' @param frac_grid A numeric vector specifying the grid of fractions to consider for cross-validation.
#'   Must be sorted in increasing order and be between 0 and 1. Default is \code{seq(0.1, 1.0, by = 0.1)}.
#' @param tol A numeric value specifying the tolerance below which singular values are considered zero.
#'   Must be positive. Default is \code{1e-10}.
#' @param cv An integer specifying the number of cross-validation folds, or an object to be used
#'   as a cross-validation generator (e.g., \code{cv = 5}). Default is \code{5}.
#' @param scoring A character string or function to evaluate the predictions on the test set.
#'   Default is \code{"r2"}.
#'
#' @return An object of class \code{fracridge_cv} containing:
#' \describe{
#'   \item{\code{best_frac}}{The fraction that achieved the best cross-validated score.}
#'   \item{\code{coef}}{A numeric array of shape \code{(p, b)} containing the estimated coefficients
#'     for the best fraction and each target.}
#'   \item{\code{alpha}}{A numeric vector containing the regularization parameters associated with the best fraction
#'     for each target.}
#'   \item{\code{cv_results}}{A data frame containing cross-validation scores for each fraction.}
#' }
#'
#' @details
#' This function performs fractional ridge regression across a grid of fractions and selects the
#' fraction that maximizes the cross-validated performance metric (e.g., R-squared).
#'
#' @examples
#' \dontrun{
#' # Generate random data
#' set.seed(1)
#' X <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' y <- rnorm(100)
#'
#' # Perform cross-validated fractional ridge regression
#' cv_result <- fracridge_cv(X, y, frac_grid = seq(0.1, 1.0, by = 0.1))
#'
#' # Print the best fraction
#' print(cv_result$best_frac)
#'
#' # Inspect coefficients
#' print(cv_result$coef)
#' }
#'
#' @export
fracridge_cv <- function(X, y, frac_grid = seq(0.1, 1.0, by = 0.1), tol = 1e-10,
                         cv = 5, scoring = "r2") {
  # Check input types
  if (!is.matrix(X)) {
    stop("X must be a numeric matrix.")
  }
  if (!(is.numeric(y) || is.matrix(y))) {
    stop("y must be a numeric vector or matrix.")
  }
  
  # Ensure y is a matrix
  if (is.vector(y)) {
    y <- matrix(y, ncol = 1)
    single_target <- TRUE
  } else {
    single_target <- FALSE
  }
  
  # Validate frac_grid
  if (!is.numeric(frac_grid) || any(frac_grid <= 0) || any(frac_grid > 1) || any(diff(frac_grid) < 0)) {
    stop("frac_grid must be a sorted numeric vector with values between 0 and 1.")
  }
  frac_grid <- as.numeric(frac_grid)
  
  n <- nrow(X)
  p <- ncol(X)
  b <- ncol(y)
  f <- length(frac_grid)
  
  # Define cross-validation folds
  set.seed(123)  # For reproducibility
  folds <- sample(rep(1:cv, length.out = n))
  
  # Initialize storage for scores
  scores <- matrix(0, nrow = f, ncol = b)
  rownames(scores) <- frac_grid
  colnames(scores) <- if (single_target) "target" else paste0("target_", 1:b)
  
  for (fold in 1:cv) {
    test_idx <- which(folds == fold)
    train_idx <- setdiff(1:n, test_idx)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx, , drop = FALSE]
    X_test <- X[test_idx, , drop = FALSE]
    y_test <- y[test_idx, , drop = FALSE]
    
    # Fit fractional ridge on training data
    model <- fracridge(X_train, y_train, fracs = frac_grid, tol = tol)
    
    for (ii in 1:f) {
      frac <- frac_grid[ii]
      preds <- predict(model, X_test, fracs = frac)
      
      if (single_target) {
        r2 <- cor(y_test, preds)^2
      } else {
        r2 <- colSums((preds - y_test)^2)
        # You can use other metrics here
        # For simplicity, using R-squared
        ss_res <- colSums((y_test - preds)^2)
        ss_tot <- colSums((y_test - colMeans(y_test))^2)
        r2 <- 1 - ss_res / ss_tot
      }
      scores[ii, ] <- scores[ii, ] + r2
    }
  }
  
  # Average scores over folds
  avg_scores <- scores / cv
  
  # Select the fraction with the highest average score
  best_frac_idx <- which.max(avg_scores)
  best_frac <- frac_grid[best_frac_idx]
  
  # Fit the model on the entire dataset with the best fraction
  final_model <- fracridge(X, y, fracs = best_frac, tol = tol)
  
  # Extract coefficients and alphas
  coef_best <- final_model$coef
  alpha_best <- final_model$alpha
  
  if (single_target) {
    coef_best <- coef_best[, , 1, drop = TRUE]
    alpha_best <- alpha_best[, 1]
  }
  
  # Compile cross-validation results
  cv_results <- as.data.frame(avg_scores)
  cv_results$frac <- frac_grid
  cv_results <- cv_results[, c("frac", if (single_target) "target" else paste0("target_", 1:b))]
  
  # Return as a list with class 'fracridge_cv'
  result <- list(
    best_frac = best_frac,
    coef = coef_best,
    alpha = alpha_best,
    cv_results = cv_results
  )
  
  class(result) <- "fracridge_cv"
  return(result)
}

#' Print Method for fracridge_cv Objects
#'
#' @param x An object of class \code{fracridge_cv}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \dontrun{
#' # Generate random data
#' set.seed(1)
#' X <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' y <- rnorm(100)
#'
#' # Perform cross-validated fractional ridge regression
#' cv_result <- fracridge_cv(X, y, frac_grid = seq(0.1, 1.0, by = 0.1))
#'
#' # Print the result
#' print(cv_result)
#' }
#'
#' @export
print.fracridge_cv <- function(x, ...) {
  cat("Cross-Validated Fractional Ridge Regression\n")
  cat("Best fraction:", x$best_frac, "\n")
  cat("Regularization parameters for the best fraction:\n")
  print(x$alpha)
  cat("Cross-Validation Results:\n")
  print(x$cv_results)
  invisible(x)
}

#' Calculate Vector Length
#'
#' Computes the Euclidean (L2) norm of a vector or the norm across a specified axis of a matrix.
#'
#' @param vec A numeric vector or matrix.
#' @param axis An integer specifying the axis to compute the norm across.
#'   \code{1} for columns and \code{2} for rows. Default is \code{1}.
#'
#' @return A numeric vector containing the Euclidean norms.
#'
#' @examples
#' vec <- c(3, 4)
#' vec_len(vec)  # Should return 5
#'
#' mat <- matrix(1:6, nrow = 2, ncol = 3)
#' vec_len(mat, axis = 1)  # Norms of rows
#' vec_len(mat, axis = 2)  # Norms of columns
#'
#' @export
vec_len <- function(vec, axis = 1) {
  if (!is.matrix(vec)) {
    return(sqrt(sum(vec^2)))
  }
  
  if (axis == 1) {
    return(sqrt(colSums(vec^2)))
  } else if (axis == 2) {
    return(sqrt(rowSums(vec^2)))
  } else {
    stop("axis must be 1 (columns) or 2 (rows).")
  }
}