# Tests for Group Analysis Framework

test_that("group_data constructor works with different formats", {
  skip_if_not_installed("fmristore")
  
  # Test CSV format detection
  csv_data <- data.frame(
    subject = rep(c("s1", "s2", "s3"), each = 2),
    roi = rep(c("roi1", "roi2"), 3),
    beta = rnorm(6),
    se = runif(6, 0.1, 0.5)
  )
  
  gd <- group_data(csv_data, format = "csv",
                   effect_cols = c(beta = "beta", se = "se"),
                   subject_col = "subject",
                   roi_col = "roi")
  
  expect_s3_class(gd, "group_data_csv")
  expect_s3_class(gd, "group_data")
  expect_equal(n_subjects(gd), 3)
  expect_equal(get_rois(gd), c("roi1", "roi2"))
})

test_that("group_data_csv handles ROI data correctly", {
  # Create test data
  test_data <- data.frame(
    subject = rep(paste0("sub-", 1:10), each = 3),
    roi = rep(c("ACC", "PFC", "V1"), 10),
    beta = rnorm(30, mean = 0.5),
    se = runif(30, 0.1, 0.3),
    group = rep(c("control", "patient"), each = 15),
    age = rep(rnorm(10, 40, 10), each = 3)
  )
  
  gd <- group_data_from_csv(
    test_data,
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    covariate_cols = c("group", "age")
  )
  
  expect_equal(n_subjects(gd), 10)
  expect_equal(length(get_rois(gd)), 3)
  expect_equal(names(get_covariates(gd)), c("group", "age"))
  
  # Test data extraction
  acc_data <- extract_csv_data(gd, roi = "ACC")
  expect_equal(length(acc_data$beta), 10)
  expect_equal(length(acc_data$se), 10)
})

test_that("fmri_meta fits basic meta-analysis models", {
  # Create simple test data
  n_subjects <- 20
  test_data <- data.frame(
    subject = paste0("sub-", 1:n_subjects),
    roi = "whole_brain",
    beta = rnorm(n_subjects, mean = 0.3, sd = 0.5),
    se = runif(n_subjects, 0.1, 0.3),
    group = factor(rep(c("A", "B"), each = 10))
  )
  
  gd <- group_data_from_csv(
    test_data,
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    covariate_cols = "group"
  )
  
  # Test fixed-effects model
  fit_fe <- fmri_meta(gd, method = "fe", verbose = FALSE)
  expect_s3_class(fit_fe, "fmri_meta")
  expect_equal(fit_fe$method, "fe")
  expect_true(all(fit_fe$tau2 == 0 | is.na(fit_fe$tau2)))
  
  # Test random-effects model
  fit_re <- fmri_meta(gd, method = "pm", verbose = FALSE)
  expect_s3_class(fit_re, "fmri_meta")
  expect_equal(fit_re$method, "pm")
  
  # Test with covariates
  fit_cov <- fmri_meta(gd, formula = ~ 1 + group, method = "fe", verbose = FALSE)
  expect_equal(ncol(fit_cov$coefficients), 2)
  expect_true("(Intercept)" %in% colnames(fit_cov$coefficients))
  expect_true("groupB" %in% colnames(fit_cov$coefficients))
})

test_that("meta-analysis methods work correctly", {
  # Create test data with known effect
  n_subjects <- 30
  true_effect <- 0.5
  test_data <- data.frame(
    subject = paste0("sub-", 1:n_subjects),
    roi = "test_roi",
    beta = rnorm(n_subjects, mean = true_effect, sd = 0.2),
    se = rep(0.2, n_subjects)
  )
  
  gd <- group_data_from_csv(
    test_data,
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi"
  )
  
  # Fit model
  fit <- fmri_meta(gd, method = "fe", verbose = FALSE)
  
  # Test coefficient extraction
  coefs <- coef(fit)
  expect_equal(dim(coefs), c(1, 1))
  expect_true(abs(coefs[1, 1] - true_effect) < 0.1)  # Should be close to true effect
  
  # Test standard errors
  ses <- se(fit)
  expect_equal(dim(ses), c(1, 1))
  expect_true(ses[1, 1] > 0)
  
  # Test z-scores
  z <- zscores(fit)
  expect_equal(dim(z), c(1, 1))
  expect_equal(z[1, 1], coefs[1, 1] / ses[1, 1])
  
  # Test p-values
  p <- pvalues(fit)
  expect_equal(dim(p), c(1, 1))
  expect_true(p[1, 1] >= 0 && p[1, 1] <= 1)
})

test_that("contrast application works", {
  # Create test data with two groups
  n_per_group <- 15
  test_data <- data.frame(
    subject = paste0("sub-", 1:(2 * n_per_group)),
    roi = "test_roi",
    beta = c(rnorm(n_per_group, 0.3, 0.2), rnorm(n_per_group, 0.6, 0.2)),
    se = rep(0.15, 2 * n_per_group),
    group = factor(rep(c("control", "treatment"), each = n_per_group))
  )
  
  gd <- group_data_from_csv(
    test_data,
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    covariate_cols = "group"
  )
  
  fit <- fmri_meta(gd, formula = ~ 1 + group, method = "fe", verbose = FALSE)
  
  # Test contrast with named vector
  contrast_result <- contrast(fit, c("grouptreatment" = 1))
  expect_s3_class(contrast_result, "fmri_meta_contrast")
  expect_true(length(contrast_result$estimate) == 1)
  expect_true(contrast_result$se > 0)
})

test_that("heterogeneity statistics are computed correctly", {
  # Create data with heterogeneity
  n_subjects <- 25
  test_data <- data.frame(
    subject = paste0("sub-", 1:n_subjects),
    roi = "test_roi",
    beta = rnorm(n_subjects, mean = 0.4, sd = 0.5),  # High between-study variance
    se = runif(n_subjects, 0.05, 0.15)  # Low within-study variance
  )
  
  gd <- group_data_from_csv(
    test_data,
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi"
  )
  
  # Random-effects should detect heterogeneity
  fit_re <- fmri_meta(gd, method = "pm", verbose = FALSE)
  expect_true(fit_re$tau2[1] > 0)  # Should detect between-study variance
  expect_true(fit_re$I2[1] > 0)     # Should have some heterogeneity
  
  # Fixed-effects should still compute Q and I2
  fit_fe <- fmri_meta(gd, method = "fe", verbose = FALSE)
  expect_equal(fit_fe$tau2[1], 0)   # No tau2 for fixed-effects
  expect_true(!is.na(fit_fe$Q[1]))  # Q should be computed
})

test_that("robust estimation down-weights outliers", {
  skip("Robust estimation test - requires full implementation")
  
  # Create data with outliers
  n_subjects <- 20
  test_data <- data.frame(
    subject = paste0("sub-", 1:n_subjects),
    roi = "test_roi",
    beta = c(rnorm(18, 0.3, 0.1), 2.5, -2.0),  # Two outliers
    se = rep(0.1, n_subjects)
  )
  
  gd <- group_data_from_csv(
    test_data,
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi"
  )
  
  # Compare robust vs non-robust
  fit_standard <- fmri_meta(gd, method = "fe", robust = "none", verbose = FALSE)
  fit_robust <- fmri_meta(gd, method = "fe", robust = "huber", verbose = FALSE)
  
  # Robust estimate should be closer to the non-outlier mean (0.3)
  expect_true(abs(coef(fit_robust)[1, 1] - 0.3) < 
              abs(coef(fit_standard)[1, 1] - 0.3))
})

test_that("print and summary methods work", {
  # Create minimal test data
  test_data <- data.frame(
    subject = paste0("sub-", 1:10),
    roi = "test",
    beta = rnorm(10),
    se = runif(10, 0.1, 0.3)
  )
  
  gd <- group_data_from_csv(
    test_data,
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi"
  )
  
  fit <- fmri_meta(gd, method = "fe", verbose = FALSE)
  
  # Test print
  expect_output(print(fit), "fMRI Meta-Analysis Results")
  expect_output(print(fit), "Method: fe")
  
  # Test summary
  expect_output(summary(fit), "fMRI Meta-Analysis Summary")
  expect_output(summary(fit), "Coefficients:")
})

test_that("tidy method returns proper tibble", {
  skip_if_not_installed("tibble")
  
  test_data <- data.frame(
    subject = rep(paste0("sub-", 1:10), each = 2),
    roi = rep(c("ACC", "PFC"), 10),
    beta = rnorm(20),
    se = runif(20, 0.1, 0.3),
    group = rep(c("A", "B"), each = 10)
  )
  
  gd <- group_data_from_csv(
    test_data,
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    covariate_cols = "group"
  )
  
  fit <- fmri_meta(gd, formula = ~ 1 + group, method = "fe", verbose = FALSE)
  
  tidy_result <- tidy(fit)
  expect_s3_class(tidy_result, "tbl_df")
  expect_true("term" %in% names(tidy_result))
  expect_true("estimate" %in% names(tidy_result))
  expect_true("std.error" %in% names(tidy_result))
  expect_true("p.value" %in% names(tidy_result))
})