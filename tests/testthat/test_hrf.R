library(testthat)

test_that("HRF_GAMMA has correct structure and properties", {
  # Test basic structure
  expect_true(inherits(HRF_GAMMA, "HRF"))
  expect_equal(attr(HRF_GAMMA, "name"), "gamma")
  expect_equal(attr(HRF_GAMMA, "param_names"), c("shape", "rate"))
  
  # Test function evaluation
  t <- seq(0, 20, by=0.5)
  result <- HRF_GAMMA(t)
  expect_true(is.numeric(result))
  expect_equal(length(result), length(t))
  expect_true(all(result >= 0))  # Gamma HRF should be non-negative
})

test_that("HRF_SPMG1 has correct structure and properties", {
  # Test basic structure
  expect_true(inherits(HRF_SPMG1, "HRF"))
  expect_equal(attr(HRF_SPMG1, "name"), "SPMG1")
  expect_equal(attr(HRF_SPMG1, "param_names"), c("A1", "A2"))
  
  # Test function evaluation
  t <- seq(0, 30, by=0.5)
  result <- HRF_SPMG1(t)
  expect_true(is.numeric(result))
  expect_equal(length(result), length(t))
  expect_equal(result[t < 0], rep(0, sum(t < 0)))  # Should be 0 for negative time
  
  # Test peak timing (should peak around 5-6 seconds)
  peak_time <- t[which.max(result)]
  expect_true(peak_time >= 4 && peak_time <= 7)
})

test_that("HRF_SPMG2 has correct structure and properties", {
  # Test basic structure
  expect_true(inherits(HRF_SPMG2, "HRF"))
  expect_equal(attr(HRF_SPMG2, "name"), "SPMG2")
  expect_equal(nbasis(HRF_SPMG2), 2)  # Should have 2 basis functions
  
  # Test function evaluation
  t <- seq(0, 30, by=0.5)
  result <- HRF_SPMG2(t)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), length(t))
  expect_equal(ncol(result), 2)  # Should return 2 columns for canonical and temporal derivative
})

test_that("HRF_GAUSSIAN has correct structure and properties", {
  # Test basic structure
  expect_true(inherits(HRF_GAUSSIAN, "HRF"))
  expect_equal(attr(HRF_GAUSSIAN, "name"), "gaussian")
  expect_equal(attr(HRF_GAUSSIAN, "param_names"), c("mean", "sd"))
  
  # Test function evaluation
  t <- seq(0, 20, by=0.5)
  result <- HRF_GAUSSIAN(t)
  expect_true(is.numeric(result))
  expect_equal(length(result), length(t))
  expect_true(all(result >= 0))  # Gaussian HRF should be non-negative
})

test_that("HRF_BSPLINE has correct structure and properties", {
  # Test basic structure
  expect_true(inherits(HRF_BSPLINE, "HRF"))
  expect_equal(attr(HRF_BSPLINE, "name"), "bspline")
  expect_equal(nbasis(HRF_BSPLINE), 5)  # Default number of basis functions
  
  # Test function evaluation
  t <- seq(0, 20, by=0.5)
  result <- HRF_BSPLINE(t)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), length(t))
  expect_equal(ncol(result), 5)  # Should return 5 columns for basis functions
})

test_that("evaluate.HRF handles different duration scenarios", {
  t <- seq(0, 20, by=0.2)
  
  # Test zero duration
  result1 <- evaluate(HRF_SPMG1, t, duration=0)
  expect_true(is.numeric(result1))
  expect_equal(length(result1), length(t))
})

test_that("gen_hrf handles lag and width correctly", {
  # Test lag
  hrf_lag <- gen_hrf(HRF_SPMG1, lag = 2)
  t <- seq(0, 20, by = 0.5)
  result_lag <- hrf_lag(t)
  result_no_lag <- HRF_SPMG1(t)
  
  # Peak should be shifted by lag
  peak_lag <- t[which.max(result_lag)]
  peak_no_lag <- t[which.max(result_no_lag)]
  expect_equal(peak_lag - peak_no_lag, 2)
  
  # Test width (block duration)
  hrf_block <- gen_hrf(HRF_SPMG1, width = 3)
  result_block <- hrf_block(t)
  
  # Block HRF should have wider response
  width_block <- sum(result_block > 0)
  width_no_block <- sum(result_no_lag > 0)
  expect_true(width_block > width_no_block)
  
  # Test combined lag and width
  hrf_both <- gen_hrf(HRF_SPMG1, lag = 2, width = 3)
  result_both <- hrf_both(t)
  peak_both <- t[which.max(result_both)]
  expect_true(peak_both > peak_no_lag)
})

test_that("gen_hrf_set combines HRFs correctly", {
  # Create basis set
  hrf1 <- gen_hrf(HRF_SPMG1, lag = 0)
  hrf2 <- gen_hrf(HRF_SPMG1, lag = 2)
  hrf3 <- gen_hrf(HRF_SPMG1, lag = 4)
  hrf_set <- gen_hrf_set(hrf1, hrf2, hrf3, name = "test_set")
  
  # Test structure
  expect_true(inherits(hrf_set, "HRF"))
  expect_equal(nbasis(hrf_set), 3)
  expect_equal(attr(hrf_set, "name"), "test_set")
  
  # Test evaluation
  t <- seq(0, 20, by = 0.5)
  result <- hrf_set(t)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(length(t), 3))
  
  # Test peaks are correctly shifted
  peaks <- apply(result, 2, function(x) t[which.max(x)])
  expect_equal(diff(peaks), c(2, 2))
})

test_that("evaluate.HRF handles different durations and summation correctly", {
  t <- seq(0, 20, by = 0.2)
  
  # Test non-zero duration
  result_dur <- evaluate(HRF_SPMG1, t, duration = 2)
  result_no_dur <- evaluate(HRF_SPMG1, t, duration = 0)
  
  # Response should be larger with duration
  expect_true(max(result_dur) > max(result_no_dur))
  
  # Test summation
  result_sum <- evaluate(HRF_SPMG1, t, duration = 2, summate = TRUE)
  result_no_sum <- evaluate(HRF_SPMG1, t, duration = 2, summate = FALSE)
  expect_false(identical(result_sum, result_no_sum))
  
  # Test precision effects
  result_fine <- evaluate(HRF_SPMG1, t, duration = 2, precision = 0.1)
  result_coarse <- evaluate(HRF_SPMG1, t, duration = 2, precision = 0.5)
  expect_false(identical(result_fine, result_coarse))
})

test_that("gen_empirical_hrf creates valid HRF", {
  # Create simple empirical HRF
  t <- seq(0, 20, by = 0.5)
  y <- dnorm(t, mean = 6, sd = 2)
  hrf <- gen_empirical_hrf(t, y, name = "test_empirical")
  
  # Test structure
  expect_true(inherits(hrf, "HRF"))
  expect_equal(attr(hrf, "name"), "test_empirical")
  expect_equal(nbasis(hrf), 1)
  
  # Test interpolation
  new_t <- seq(0, 20, by = 0.3)
  result <- hrf(new_t)
  expect_equal(length(result), length(new_t))
  expect_true(all(result >= 0))
  
  # Test extrapolation
  extended_t <- c(-2, t, 22)
  result_ext <- hrf(extended_t)
  expect_equal(result_ext[1], 0)  # Left extrapolation
  expect_equal(result_ext[length(result_ext)], 0)  # Right extrapolation
})

test_that("HRF objects maintain correct attributes", {
  # Test basic HRF attributes
  t <- seq(0, 20, by = 0.5)
  
  hrfs <- list(
    HRF_SPMG1 = HRF_SPMG1,
    HRF_SPMG2 = HRF_SPMG2,
    HRF_GAMMA = HRF_GAMMA,
    HRF_GAUSSIAN = HRF_GAUSSIAN
  )
  
  for (name in names(hrfs)) {
    hrf <- hrfs[[name]]
    expect_true(inherits(hrf, "HRF"))
    expect_true(is.function(hrf))
    expect_true(!is.null(attr(hrf, "span")))
    expect_true(!is.null(attr(hrf, "nbasis")))
    expect_true(!is.null(attr(hrf, "name")))
    
    # Test evaluation produces correct dimensions
    result <- hrf(t)
    if (attr(hrf, "nbasis") == 1) {
      expect_true(is.numeric(result))
      expect_equal(length(result), length(t))
    } else {
      expect_true(is.matrix(result))
      expect_equal(nrow(result), length(t))
      expect_equal(ncol(result), attr(hrf, "nbasis"))
    }
  }
})
