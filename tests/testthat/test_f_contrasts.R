# Test F-contrasts and statistical inference
# Testing compute_f_statistic() and related functionality

library(fmrireg)
library(testthat)

# Helper function to create factorial design data
create_factorial_data <- function(n_per_cell = 20, sigma = 1) {
  # 2x2 factorial design
  factor_a <- rep(c("A1", "A2"), each = n_per_cell * 2)
  factor_b <- rep(rep(c("B1", "B2"), each = n_per_cell), 2)
  
  # Design matrix with interaction
  X <- model.matrix(~ factor_a * factor_b)
  
  # True parameters: intercept, main A, main B, interaction
  beta_true <- c(10, 2, -1, 0.5)
  
  # Generate data
  n <- length(factor_a)
  Y <- X %*% beta_true + rnorm(n, sd = sigma)
  
  list(
    X = X,
    Y = Y,
    factor_a = factor_a,
    factor_b = factor_b,
    beta_true = beta_true,
    n = n
  )
}

# Helper to fit model using integrated solver
fit_model <- function(X, Y, config = NULL) {
  if (is.null(config)) {
    config <- fmrireg:::fmri_lm_config(robust = FALSE)
  }
  
  # Ensure Y is matrix
  if (!is.matrix(Y)) {
    Y <- matrix(Y, ncol = 1)
  }
  
  fmrireg:::solve_integrated_glm(
    X = X,
    Y = Y,
    config = config
  )
}

test_that("F-statistic computation for factorial designs", {
  set.seed(123)
  dat <- create_factorial_data(n_per_cell = 30)
  
  # Fit model
  fit_result <- fit_model(dat$X, dat$Y)
  
  # Test main effect of factor A (contrasts rows 2)
  contrast_a <- matrix(c(0, 1, 0, 0), nrow = 1)
  f_result_a <- fmrireg:::compute_f_statistic(fit_result, contrast_a)
  
  expect_true(!is.null(f_result_a$fstat))
  expect_true(!is.null(f_result_a$pvalue))
  expect_equal(f_result_a$df1, 1)
  expect_equal(f_result_a$df2, dat$n - ncol(dat$X))
  
  # F-statistic should be positive
  expect_true(all(f_result_a$fstat > 0))
  
  # Test main effect of factor B
  contrast_b <- matrix(c(0, 0, 1, 0), nrow = 1)
  f_result_b <- fmrireg:::compute_f_statistic(fit_result, contrast_b)
  
  expect_true(all(f_result_b$fstat > 0))
  
  # Test interaction
  contrast_int <- matrix(c(0, 0, 0, 1), nrow = 1)
  f_result_int <- fmrireg:::compute_f_statistic(fit_result, contrast_int)
  
  expect_true(all(f_result_int$fstat > 0))
})

test_that("Multi-row contrast matrices", {
  set.seed(234)
  dat <- create_factorial_data(n_per_cell = 25)
  fit_result <- fit_model(dat$X, dat$Y)
  
  # Test both main effects simultaneously (2 df test)
  contrast_main <- rbind(
    c(0, 1, 0, 0),  # Main effect A
    c(0, 0, 1, 0)   # Main effect B
  )
  
  f_result <- fmrireg:::compute_f_statistic(fit_result, contrast_main)
  
  expect_equal(f_result$df1, 2)  # 2 restrictions
  expect_equal(f_result$df2, dat$n - ncol(dat$X))
  expect_true(all(f_result$fstat > 0))
  
  # Test all effects (3 df test)
  contrast_all <- rbind(
    c(0, 1, 0, 0),  # Main effect A
    c(0, 0, 1, 0),  # Main effect B  
    c(0, 0, 0, 1)   # Interaction
  )
  
  f_result_all <- fmrireg:::compute_f_statistic(fit_result, contrast_all)
  
  expect_equal(f_result_all$df1, 3)
  expect_true(all(f_result_all$fstat > 0))
})

test_that("Single-row F-contrasts match t-statistics squared", {
  set.seed(345)
  n <- 100
  X <- cbind(1, rnorm(n), rnorm(n))
  beta <- c(2, 1, -0.5)
  Y <- X %*% beta + rnorm(n, sd = 1.5)
  
  fit_result <- fit_model(X, Y)
  
  # Test single parameter
  contrast <- matrix(c(0, 1, 0), nrow = 1)
  
  # Compute t-statistic via compute_contrast
  t_result <- fmrireg:::compute_contrast(fit_result, contrast)
  
  # Compute F-statistic
  f_result <- fmrireg:::compute_f_statistic(fit_result, contrast)
  
  # F = t^2 for single df in numerator
  expect_equal(f_result$fstat[1], t_result$tstat[1]^2, tolerance = 1e-10)
  
  # p-values should match
  p_from_t <- 2 * pt(-abs(t_result$tstat[1]), df = fit_result$dfres)
  p_from_f <- pf(f_result$fstat[1], 1, fit_result$dfres, lower.tail = FALSE)
  expect_equal(p_from_t, p_from_f, tolerance = 1e-10)
})

test_that("Singular contrast matrix handling", {
  set.seed(456)
  dat <- create_factorial_data()
  fit_result <- fit_model(dat$X, dat$Y)
  
  # Create singular contrast matrix (second row is multiple of first)
  contrast_singular <- rbind(
    c(0, 1, 0, 0),
    c(0, 2, 0, 0)  # 2 * first row
  )
  
  # Should handle with warning
  expect_warning(
    f_result <- fmrireg:::compute_f_statistic(fit_result, contrast_singular),
    "singular"
  )
  
  # Result should still be computed using generalized inverse
  expect_true(!is.null(f_result$fstat))
  expect_true(all(f_result$fstat >= 0))
  
  # Create truly independent but nearly singular matrix
  contrast_near_singular <- rbind(
    c(0, 1, 0, 0),
    c(0, 1.0001, 0, 0)  # Almost the same
  )
  
  # Should work but might warn
  f_result2 <- suppressWarnings(
    fmrireg:::compute_f_statistic(fit_result, contrast_near_singular)
  )
  
  expect_true(!is.null(f_result2$fstat))
})

test_that("Validation against stats::anova() for known cases", {
  set.seed(567)
  dat <- create_factorial_data(n_per_cell = 40)
  
  # Fit using lm for comparison
  df <- data.frame(
    Y = dat$Y,
    factor_a = dat$factor_a,
    factor_b = dat$factor_b
  )
  
  # For intercept + single factor model, results should match anova()
  lm_simple <- lm(Y ~ factor_a, data = df)
  anova_simple <- anova(lm_simple)
  
  # Our F-test for just factor_a in simple model
  X_simple <- model.matrix(~ factor_a, data = df)
  fit_simple <- fit_model(X_simple, dat$Y)
  contrast_simple <- matrix(c(0, 1), nrow = 1)
  f_simple <- fmrireg:::compute_f_statistic(fit_simple, contrast_simple)
  
  expect_equal(f_simple$fstat[1], anova_simple["factor_a", "F value"], 
               tolerance = 0.001)
})

test_that("Integration with existing contrast infrastructure", {
  set.seed(678)
  dat <- create_factorial_data()
  fit_result <- fit_model(dat$X, dat$Y)
  
  # Test compute_voxelwise_contrasts with F-tests
  contrast_list <- list(
    # Single row (t-test)
    t_test = matrix(c(0, 1, 0, 0), nrow = 1),
    # Multi-row (F-test) 
    f_test = rbind(
      c(0, 1, 0, 0),
      c(0, 0, 1, 0)
    )
  )
  
  results <- fmrireg:::compute_voxelwise_contrasts(fit_result, contrast_list)
  
  expect_length(results, 2)
  
  # First should be t-test result
  expect_true("tstat" %in% names(results[[1]]))
  
  # Second should be F-test result
  expect_true("fstat" %in% names(results[[2]]))
  expect_equal(results[[2]]$df1, 2)
})

test_that("Edge cases for F-contrasts", {
  set.seed(789)
  n <- 50
  X <- cbind(1, rnorm(n))
  Y <- X %*% c(1, 2) + rnorm(n)
  
  fit_result <- fit_model(X, Y)
  
  # Empty contrast matrix - will error in matrix multiplication
  expect_error(
    fmrireg:::compute_f_statistic(fit_result, matrix(nrow = 0, ncol = 2))
  )
  
  # Wrong number of columns - will error in matrix multiplication
  expect_error(
    fmrireg:::compute_f_statistic(fit_result, matrix(c(1, 2, 3), nrow = 1))
  )
  
  # All zero contrast
  contrast_zero <- matrix(c(0, 0), nrow = 1)
  f_zero <- suppressWarnings(
    fmrireg:::compute_f_statistic(fit_result, contrast_zero)
  )
  
  # F-stat should be 0 for null contrast
  expect_equal(f_zero$fstat[1], 0, tolerance = 1e-10)
  expect_equal(f_zero$pvalue[1], 1, tolerance = 1e-10)
})

test_that("F-tests with rank-deficient contrasts", {
  set.seed(890)
  # Create rank-deficient design
  n <- 60
  group <- rep(1:3, each = 20)
  X <- model.matrix(~ 0 + factor(group))  # No intercept, 3 groups
  
  # Add intercept to make it rank-deficient
  X <- cbind(1, X)
  
  # True effects
  Y <- 10 + 2*(group == 2) - 1*(group == 3) + rnorm(n)
  
  # Fit will handle rank deficiency
  fit_result <- fit_model(X, Y)
  
  # Test contrasts that respect the constraint
  # Compare group 2 vs group 1
  contrast_2v1 <- matrix(c(0, -1, 1, 0), nrow = 1)
  
  f_result <- suppressWarnings(
    fmrireg:::compute_f_statistic(fit_result, contrast_2v1)
  )
  
  expect_true(!is.null(f_result$fstat))
  expect_true(f_result$fstat[1] > 0)  # Should detect the difference
})

test_that("F-contrasts with multiple voxels", {
  set.seed(901)
  dat <- create_factorial_data(n_per_cell = 20)
  
  # Create multi-voxel data
  n_voxels <- 10
  Y_multi <- matrix(nrow = nrow(dat$X), ncol = n_voxels)
  
  # Different effect sizes per voxel
  for (v in 1:n_voxels) {
    effect_scale <- v / n_voxels
    beta_v <- dat$beta_true * effect_scale
    Y_multi[, v] <- dat$X %*% beta_v + rnorm(nrow(dat$X))
  }
  
  # Fit all voxels
  fit_result <- fit_model(dat$X, Y_multi)
  
  # Test main effects jointly
  contrast_main <- rbind(
    c(0, 1, 0, 0),
    c(0, 0, 1, 0)
  )
  
  f_result <- fmrireg:::compute_f_statistic(fit_result, contrast_main)
  
  expect_equal(ncol(f_result$fstat), n_voxels)
  expect_equal(ncol(f_result$pvalue), n_voxels)
  
  # F-stats should generally increase with effect size (but allow for some variation)
  # Check that correlation between effect size and F-stat is positive
  effect_sizes <- 1:n_voxels / n_voxels
  cor_f_effect <- cor(f_result$fstat[1,], effect_sizes)
  expect_true(cor_f_effect > 0.5)  # Strong positive correlation
  
  # All should be non-negative
  expect_true(all(f_result$fstat >= 0))
})

test_that("F-contrasts with robust fitting", {
  set.seed(102)
  dat <- create_factorial_data(n_per_cell = 25)
  
  # Add outliers
  outlier_idx <- c(5, 15, 25, 35)
  dat$Y[outlier_idx] <- dat$Y[outlier_idx] + 10 * c(1, -1, 1, -1)
  
  # Fit with robust estimation
  config_robust <- fmrireg:::fmri_lm_config(
    robust = "huber",
    ar_options = list(cor_struct = "none")
  )
  
  fit_robust <- fit_model(dat$X, dat$Y, config_robust)
  
  # Test main effects
  contrast_main <- rbind(
    c(0, 1, 0, 0),
    c(0, 0, 1, 0)
  )
  
  f_result <- fmrireg:::compute_f_statistic(fit_robust, contrast_main)
  
  # Should work with robust fits
  expect_true(!is.null(f_result$fstat))
  expect_true(all(f_result$fstat > 0))
  
  # Check that effective_df exists (robust fitting may adjust it)
  expect_true(!is.null(fit_robust$effective_df) || !is.null(fit_robust$dfres))
})

test_that("F-contrasts match compute_contrast for single rows", {
  set.seed(103)
  n <- 80
  X <- cbind(1, rnorm(n), rnorm(n), rnorm(n))
  beta <- c(1, 2, -1, 0.5)
  Y <- X %*% beta + rnorm(n, sd = 2)
  
  fit_result <- fit_model(X, Y)
  
  # Test each parameter
  for (j in 2:ncol(X)) {
    contrast <- matrix(0, nrow = 1, ncol = ncol(X))
    contrast[1, j] <- 1
    
    # Get t-test result
    t_result <- fmrireg:::compute_contrast(fit_result, contrast)
    
    # Get F-test result  
    f_result <- fmrireg:::compute_f_statistic(fit_result, contrast)
    
    # Verify F = t^2
    expect_equal(f_result$fstat[1], t_result$tstat[1]^2, 
                 tolerance = 1e-10,
                 info = paste("Failed for parameter", j))
    
    # Verify df
    expect_equal(f_result$df1, 1)
    expect_equal(f_result$df2, t_result$df)
  }
})

test_that("Orthogonal polynomial contrasts", {
  set.seed(104)
  # Ordered factor with 4 levels
  n_per_level <- 25
  levels <- rep(1:4, each = n_per_level)
  X <- model.matrix(~ factor(levels))
  
  # Linear trend in means: 0, 1, 2, 3
  Y <- levels - 1 + rnorm(length(levels), sd = 0.5)
  
  fit_result <- fit_model(X, Y)
  
  # Check dimensions match
  expect_equal(ncol(X), 4)  # Intercept + 3 dummy variables
  
  # Orthogonal polynomial contrasts for 4 levels (excluding intercept)
  # Linear: -3, -1, 1, 3
  # Quadratic: 1, -1, -1, 1
  # Cubic: -1, 3, -3, 1
  
  contrast_poly <- rbind(
    c(0, -3, -1, 1),    # Linear (adjusted for 4 columns)
    c(0, 1, -1, -1),    # Quadratic (adjusted)
    c(0, -1, 3, -3)     # Cubic (adjusted)
  ) / sqrt(20)  # Normalize
  
  f_poly <- fmrireg:::compute_f_statistic(fit_result, contrast_poly)
  
  expect_equal(f_poly$df1, 3)
  expect_true(f_poly$fstat[1] > 0)
  
  # Test just linear trend
  contrast_linear <- matrix(c(0, -3, -1, 1) / sqrt(20), nrow = 1)
  f_linear <- fmrireg:::compute_f_statistic(fit_result, contrast_linear)
  
  # Linear trend should be highly significant given the data
  expect_true(f_linear$pvalue[1] < 0.001)
})

test_that("F-contrasts handle numerical edge cases", {
  set.seed(105)
  n <- 50
  X <- cbind(1, rnorm(n))
  
  # Test with very small variance
  Y_small_var <- X %*% c(1, 2) + rnorm(n, sd = 0.001)
  fit_small <- fit_model(X, Y_small_var)
  
  contrast <- matrix(c(0, 1), nrow = 1)
  f_result <- fmrireg:::compute_f_statistic(fit_small, contrast)
  
  # Should get very large F-statistic
  expect_true(f_result$fstat[1] > 100)
  expect_true(f_result$pvalue[1] < 0.001)
  
  # Test with no variance (perfect fit)
  Y_no_var <- X %*% c(1, 2)
  fit_perfect <- fit_model(X, Y_no_var)
  
  # Sigma2 should be very small
  expect_true(fit_perfect$sigma2 < 1e-10)
  
  # F-stat computation should still work
  f_perfect <- fmrireg:::compute_f_statistic(fit_perfect, contrast)
  expect_true(is.finite(f_perfect$fstat[1]))
})

test_that("F-contrasts with AR corrections", {
  set.seed(106)
  n <- 200
  X <- cbind(1, rnorm(n))
  
  # Generate AR(1) errors
  ar_coef <- 0.7
  e <- arima.sim(list(ar = ar_coef), n = n)
  Y <- X %*% c(1, 3) + as.vector(e)
  
  # Fit with AR correction
  config_ar <- fmrireg:::fmri_lm_config(
    robust = FALSE,
    ar_options = list(cor_struct = "ar1", iter = 2)
  )
  
  fit_ar <- fit_model(X, Y, config_ar)
  
  # Test slope parameter
  contrast <- matrix(c(0, 1), nrow = 1)
  f_result <- fmrireg:::compute_f_statistic(fit_ar, contrast)
  
  # Should detect the effect even with AR structure
  expect_true(f_result$pvalue[1] < 0.001)
  
  # Compare with OLS (no AR correction)
  fit_ols <- fit_model(X, Y)
  f_ols <- fmrireg:::compute_f_statistic(fit_ols, contrast)
  
  # AR-corrected model should have different (usually larger) standard errors
  expect_true(abs(f_result$fstat[1] - f_ols$fstat[1]) > 0.1)
})

test_that("Compute_voxelwise_contrasts correctly routes to F-tests", {
  set.seed(107)
  dat <- create_factorial_data(n_per_cell = 15)
  fit_result <- fit_model(dat$X, dat$Y)
  
  # Mixed list of t and F contrasts
  contrast_list <- list(
    # Single row - should use compute_contrast
    t_test_a = matrix(c(0, 1, 0, 0), nrow = 1),
    
    # Multi-row - should use compute_f_statistic
    f_test_main = rbind(
      c(0, 1, 0, 0),
      c(0, 0, 1, 0)
    ),
    
    # Another single row
    t_test_b = matrix(c(0, 0, 1, 0), nrow = 1),
    
    # Three-row F-test
    f_test_all = rbind(
      c(0, 1, 0, 0),
      c(0, 0, 1, 0),
      c(0, 0, 0, 1)
    )
  )
  
  results <- fmrireg:::compute_voxelwise_contrasts(fit_result, contrast_list)
  
  # Check correct routing
  expect_true("tstat" %in% names(results[[1]]))  # t-test
  expect_true("fstat" %in% names(results[[2]]))  # F-test
  expect_true("tstat" %in% names(results[[3]]))  # t-test
  expect_true("fstat" %in% names(results[[4]]))  # F-test
  
  # Check dimensions
  expect_equal(results[[2]]$df1, 2)
  expect_equal(results[[4]]$df1, 3)
})