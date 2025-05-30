# Test bootstrap functionality

library(fmrireg)
library(testthat)

test_that("multiresponse_bootstrap_lm works with basic inputs", {
  # Create simple test data
  n <- 100
  p <- 3
  X <- cbind(1, matrix(rnorm(n * (p-1)), n, p-1))
  Y <- matrix(rnorm(n * 5), n, 5)  # 5 response variables
  
  # Create data environment
  data_env <- list(.y = Y)
  
  # Create proper contrast objects
  con1 <- list(
    term = NULL,
    name = "effect1",
    weights = c(0, 1, 0),
    condnames = NULL,
    contrast_spec = NULL
  )
  class(con1) <- c("contrast", "list")
  attr(con1, "colind") <- 1:3  # Applies to all 3 columns of design matrix
  
  con2 <- list(
    term = NULL,
    name = "effect2", 
    weights = c(0, 0, 1),
    condnames = NULL,
    contrast_spec = NULL
  )
  class(con2) <- c("contrast", "list")
  attr(con2, "colind") <- 1:3  # Applies to all 3 columns of design matrix
  
  conlist <- list(con1 = con1, con2 = con2)
  
  # Run bootstrap
  result <- fmrireg:::multiresponse_bootstrap_lm(
    modmat = X,
    data_env = data_env,
    conlist = conlist,
    vnames = c("intercept", "x1", "x2"),
    fcon = NULL,
    nboot = 10,  # Small number for testing
    block_size = 20,
    event_indices = 2:3
  )
  
  expect_true(result$bootstrap)
  expect_equal(result$nboot, 10)
  expect_true(!is.null(result$con_cov))
  expect_true(!is.null(result$beta_cov))
  expect_equal(length(result$beta_cov), 2)  # Length should match event_indices length
})

test_that("multiresponse_bootstrap_lm validates inputs correctly", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- matrix(rnorm(30), 10, 3)
  data_env <- list(.y = Y)
  
  # Test missing inputs
  expect_error(
    fmrireg:::multiresponse_bootstrap_lm(data_env = data_env, conlist = list()),
    "Either 'form' or 'modmat' must be provided"
  )
  
  # Test invalid block size
  expect_error(
    fmrireg:::multiresponse_bootstrap_lm(modmat = X, data_env = data_env, block_size = -1),
    "block_size must be a positive integer"
  )
  
  # Test invalid nboot
  expect_error(
    fmrireg:::multiresponse_bootstrap_lm(modmat = X, data_env = data_env, nboot = 0),
    "nboot must be a positive integer"
  )
  
  # Test dimension mismatch
  Y_wrong <- matrix(rnorm(20), 5, 4)  # Wrong number of rows
  expect_error(
    fmrireg:::multiresponse_bootstrap_lm(modmat = X, data_env = list(.y = Y_wrong)),
    "Number of rows in modmat"
  )
})

test_that("multiresponse_bootstrap_lm handles edge cases", {
  # Small dataset
  n <- 20
  p <- 2
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * 2), n, 2)
  
  # Block size larger than data
  result <- fmrireg:::multiresponse_bootstrap_lm(
    modmat = X,
    data_env = list(.y = Y),
    conlist = list(),
    vnames = c("intercept", "x1"),
    fcon = NULL,
    nboot = 5,
    block_size = 30,  # Larger than n
    event_indices = 2
  )
  
  expect_true(result$bootstrap)
  
  # Single response
  Y_single <- matrix(rnorm(n), n, 1)
  result_single <- fmrireg:::multiresponse_bootstrap_lm(
    modmat = X,
    data_env = list(.y = Y_single),
    conlist = list(),
    vnames = c("intercept", "x1"),
    fcon = NULL,
    nboot = 5,
    block_size = 10,
    event_indices = 2
  )
  
  expect_equal(length(result_single$beta_cov), 1)  # Should have 1 element for event_indices = 2
})

test_that("bootstrap results are stable with set seed", {
  # Create data once
  n <- 50
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * 3), n, 3)
  
  # Create proper contrast object
  con <- list(
    term = NULL,
    name = "slope",
    weights = c(0, 1),
    condnames = NULL,
    contrast_spec = NULL
  )
  class(con) <- c("contrast", "list")
  attr(con, "colind") <- 1:2
  conlist <- list(con = con)
  
  set.seed(123)
  result1 <- fmrireg:::multiresponse_bootstrap_lm(
    modmat = X,
    data_env = list(.y = Y),
    conlist = conlist,
    vnames = c("intercept", "x"),
    fcon = NULL,
    nboot = 20,
    block_size = 10,
    event_indices = 2
  )
  
  set.seed(123)
  result2 <- fmrireg:::multiresponse_bootstrap_lm(
    modmat = X,
    data_env = list(.y = Y),
    conlist = conlist,
    vnames = c("intercept", "x"),
    fcon = NULL,
    nboot = 20,
    block_size = 10,
    event_indices = 2
  )
  
  # Results should be identical with same seed
  expect_equal(result1$con_cov, result2$con_cov)
  expect_equal(result1$beta_cov, result2$beta_cov)
})