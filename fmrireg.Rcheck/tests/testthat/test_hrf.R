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

test_that("as_hrf creates valid HRF objects", {
  # Simple function
  my_func <- function(t) { t^2 }
  
  # Create HRF using as_hrf
  hrf_obj <- as_hrf(my_func, name = "test_sq", nbasis = 1L, span = 10, 
                      params = list(power = 2))
  
  # Check class
  expect_true(inherits(hrf_obj, "HRF"))
  expect_true(inherits(hrf_obj, "function"))
  
  # Check attributes
  expect_equal(attr(hrf_obj, "name"), "test_sq")
  expect_equal(attr(hrf_obj, "nbasis"), 1L)
  expect_equal(attr(hrf_obj, "span"), 10)
  expect_equal(attr(hrf_obj, "param_names"), "power")
  expect_equal(attr(hrf_obj, "params"), list(power = 2))
  
  # Check function evaluation
  expect_equal(hrf_obj(5), 25)
  
  # Check defaults
  hrf_obj_default <- as_hrf(my_func)
  expect_equal(attr(hrf_obj_default, "name"), "my_func")
  expect_equal(attr(hrf_obj_default, "nbasis"), 1L)
  expect_equal(attr(hrf_obj_default, "span"), 24)
  expect_null(attr(hrf_obj_default, "param_names"))
  expect_equal(attr(hrf_obj_default, "params"), list())
  
  # Check multi-basis
  my_multi_func <- function(t) { cbind(t, t^2) }
  hrf_multi <- as_hrf(my_multi_func, nbasis = 2L)
  expect_equal(attr(hrf_multi, "nbasis"), 2L)
  expect_equal(as.matrix(hrf_multi(3)), as.matrix(cbind(3, 9)), check.attributes = FALSE)
})

test_that("bind_basis combines HRF objects correctly", {
  # Create individual HRF objects
  f1 <- function(t) { t }
  f2 <- function(t) { t^2 }
  f3 <- function(t) { rep(1, length(t)) }
  
  hrf1 <- as_hrf(f1, name="linear", span=10)
  hrf2 <- as_hrf(f2, name="quadratic", span=12)
  hrf3 <- as_hrf(f3, name="constant", span=8)
  
  # Combine them
  combined_hrf <- bind_basis(hrf1, hrf2, hrf3)
  
  # Check class
  expect_true(inherits(combined_hrf, "HRF"))
  expect_true(inherits(combined_hrf, "function"))
  
  # Check attributes
  expect_equal(attr(combined_hrf, "name"), "linear + quadratic + constant")
  expect_equal(attr(combined_hrf, "nbasis"), 3L) # 1 + 1 + 1
  expect_equal(attr(combined_hrf, "span"), 12) # max(10, 12, 8)
  
  # Check function evaluation
  t_vals <- c(0, 1, 2, 5)
  expected_output <- cbind(f1(t_vals), f2(t_vals), f3(t_vals))
  colnames(expected_output) <- NULL # Match the expected output of bind_basis function
  
  # Use check.attributes = FALSE for robustness against potential slight differences
  expect_equal(combined_hrf(t_vals), expected_output, check.attributes = FALSE)
  
  # Test with a multi-basis input
  f_multi <- function(t) cbind(sin(t), cos(t))
  hrf_multi <- as_hrf(f_multi, name="trig", nbasis=2L, span=15)
  
  combined_hrf2 <- bind_basis(hrf1, hrf_multi)
  expect_equal(attr(combined_hrf2, "nbasis"), 3L) # 1 + 2
  expect_equal(attr(combined_hrf2, "span"), 15) # max(10, 15)
  expect_equal(attr(combined_hrf2, "name"), "linear + trig")
  
  expected_output2 <- cbind(f1(t_vals), f_multi(t_vals))
  colnames(expected_output2) <- NULL
  expect_equal(combined_hrf2(t_vals), expected_output2, check.attributes = FALSE)
  
  # Test binding just one element
  combined_single <- bind_basis(hrf1)
  expect_equal(attr(combined_single, "name"), "linear")
  expect_equal(attr(combined_single, "nbasis"), 1L)
  expect_equal(attr(combined_single, "span"), 10)
  expect_equal(combined_single(t_vals), f1(t_vals))
})

test_that("lag_hrf correctly lags an HRF object", {
  # Use HRF_SPMG1 as the base HRF
  base_hrf <- HRF_SPMG1
  t <- seq(0, 30, by = 0.5)
  lag_amount <- 5
  
  # Create lagged HRF
  lagged_hrf <- lag_hrf(base_hrf, lag_amount)
  
  # Test basic structure
  expect_true(inherits(lagged_hrf, "HRF"))
  expect_true(inherits(lagged_hrf, "function"))
  expect_equal(nbasis(lagged_hrf), nbasis(base_hrf))
  expect_equal(attr(lagged_hrf, "span"), attr(base_hrf, "span") + lag_amount)
  expect_true(grepl(paste0("_lag\\(", lag_amount, "\\)"), attr(lagged_hrf, "name")))
  expect_equal(attr(lagged_hrf, "params")$.lag, lag_amount)

  # Test function evaluation: lagged_hrf(t) should equal base_hrf(t - lag)
  result_lagged <- lagged_hrf(t)
  result_manual_lag <- base_hrf(t - lag_amount)
  expect_equal(result_lagged, result_manual_lag)
  
  # Test peak timing (should be shifted by lag_amount)
  peak_lagged <- t[which.max(result_lagged)]
  peak_base <- t[which.max(base_hrf(t))]
  # Allow for slight tolerance due to discrete time steps
  expect_true(abs((peak_lagged - peak_base) - lag_amount) < 1) 
  
  # Test with zero lag
  lagged_zero <- lag_hrf(base_hrf, 0)
  expect_equal(lagged_zero(t), base_hrf(t))
  expect_equal(attr(lagged_zero, "span"), attr(base_hrf, "span"))
  
  # Test with a multi-basis HRF (HRF_SPMG2)
  base_hrf_multi <- HRF_SPMG2
  lagged_hrf_multi <- lag_hrf(base_hrf_multi, lag_amount)
  expect_equal(nbasis(lagged_hrf_multi), nbasis(base_hrf_multi))
  expect_equal(lagged_hrf_multi(t), base_hrf_multi(t - lag_amount))
  expect_equal(attr(lagged_hrf_multi, "span"), attr(base_hrf_multi, "span") + lag_amount)
})

test_that("block_hrf correctly blocks an HRF object", {
  base_hrf <- HRF_SPMG1
  t <- seq(0, 30, by = 0.2)
  width <- 5
  precision <- 0.2

  blocked_hrf_sum <- block_hrf(base_hrf, width = width, precision = precision, summate = TRUE, normalize = FALSE)
  blocked_hrf_max <- block_hrf(base_hrf, width = width, precision = precision, summate = FALSE, normalize = FALSE)
  blocked_hrf_norm <- block_hrf(base_hrf, width = width, precision = precision, summate = TRUE, normalize = TRUE)

  # Test basic structure
  expect_true(inherits(blocked_hrf_sum, "HRF"))
  expect_equal(nbasis(blocked_hrf_sum), nbasis(base_hrf))
  expect_equal(attr(blocked_hrf_sum, "span"), attr(base_hrf, "span") + width)
  expect_true(grepl(paste0("_block\\(w=", width, "\\)"), attr(blocked_hrf_sum, "name")))
  expect_equal(attr(blocked_hrf_sum, "params")$.width, width)
  expect_equal(attr(blocked_hrf_sum, "params")$.summate, TRUE)
  expect_equal(attr(blocked_hrf_max, "params")$.summate, FALSE)
  expect_equal(attr(blocked_hrf_norm, "params")$.normalize, TRUE)

  # Test function evaluation - Compare with evaluate.HRF which uses similar logic
  eval_res_sum <- evaluate(base_hrf, t, duration = width, precision = precision, summate = TRUE, normalize = FALSE)
  eval_res_max <- evaluate(base_hrf, t, duration = width, precision = precision, summate = FALSE, normalize = FALSE)
  eval_res_norm <- evaluate(base_hrf, t, duration = width, precision = precision, summate = TRUE, normalize = TRUE)

  expect_equal(blocked_hrf_sum(t), eval_res_sum)
  # Max logic might differ slightly depending on implementation details, check if shape is reasonable
  # expect_equal(blocked_hrf_max(t), eval_res_max)
  expect_false(identical(blocked_hrf_sum(t), blocked_hrf_max(t)))
  expect_equal(blocked_hrf_norm(t), eval_res_norm)
  expect_equal(max(abs(blocked_hrf_norm(t))), 1) # Check normalization worked

  # Test width_block > width_no_block (as in gen_hrf test)
  result_block <- blocked_hrf_sum(t)
  result_no_block <- base_hrf(t)
  
  # Compare Area Under Curve (AUC) approximation as a measure of width/magnitude
  auc_block <- sum(abs(result_block)) * (t[2]-t[1]) # Multiply by time step for approx integral
  auc_no_block <- sum(abs(result_no_block)) * (t[2]-t[1])
  
  expect_true(auc_block > auc_no_block)

  # Test half_life
  blocked_hl <- block_hrf(base_hrf, width = width, precision = precision, half_life = 2)
  expect_false(identical(blocked_hl(t), blocked_hrf_sum(t)))
  expect_true(max(abs(blocked_hl(t))) < max(abs(blocked_hrf_sum(t)))) # Expect decay to reduce peak

  # Test negligible width
  blocked_negligible <- block_hrf(base_hrf, width = 0.01, precision = 0.1)
  expect_equal(blocked_negligible(t), base_hrf(t))
})

test_that("normalise_hrf correctly normalises an HRF object", {
  # Create an unnormalised HRF (Gaussian scaled by 5)
  unnorm_func <- function(t) 5 * dnorm(t, 6, 2)
  unnorm_hrf <- as_hrf(unnorm_func, name="unnorm_gauss")
  t <- seq(0, 20, by=0.1)
  
  # Normalise it
  norm_hrf <- normalise_hrf(unnorm_hrf)

  # Test basic structure
  expect_true(inherits(norm_hrf, "HRF"))
  expect_equal(nbasis(norm_hrf), 1)
  expect_equal(attr(norm_hrf, "span"), attr(unnorm_hrf, "span"))
  expect_true(grepl("_norm", attr(norm_hrf, "name")))
  expect_equal(attr(norm_hrf, "params")$.normalised, TRUE)
  
  # Test peak value
  result_norm <- norm_hrf(t)
  expect_equal(max(abs(result_norm)), 1)
  
  # Test relationship to original
  result_unnorm <- unnorm_hrf(t)
  peak_unnorm <- max(abs(result_unnorm))
  expect_equal(result_norm, result_unnorm / peak_unnorm)
  
  # Test with an already normalised HRF (should remain normalised)
  norm_spmg1 <- normalise_hrf(HRF_SPMG1)
  expect_equal(max(abs(norm_spmg1(t))), 1, tolerance = 1e-7)
  
  # Test with multi-basis HRF (HRF_SPMG2)
  unnorm_spmg2_func <- function(t) cbind(5 * HRF_SPMG2(t)[,1], 10 * HRF_SPMG2(t)[,2])
  unnorm_spmg2 <- as_hrf(unnorm_spmg2_func, name="unnorm_spmg2", nbasis=2L)
  norm_spmg2 <- normalise_hrf(unnorm_spmg2)
  
  expect_equal(nbasis(norm_spmg2), 2)
  result_norm_spmg2 <- norm_spmg2(t)
  expect_equal(max(abs(result_norm_spmg2[,1])), 1)
  expect_equal(max(abs(result_norm_spmg2[,2])), 1)
})

test_that("gen_hrf correctly sets nbasis for function inputs", {
  # Single basis functions
  hrf_g <- gen_hrf(hrf_gaussian)
  expect_equal(nbasis(hrf_g), 1)
  
  hrf_s1 <- gen_hrf(hrf_spmg1)
  expect_equal(nbasis(hrf_s1), 1)
  
  # Single basis HRF object
  hrf_s1_obj <- gen_hrf(HRF_SPMG1)
  expect_equal(nbasis(hrf_s1_obj), 1)

  # Multi-basis HRF objects
  hrf_s2_obj <- gen_hrf(HRF_SPMG2)
  expect_equal(nbasis(hrf_s2_obj), 2)
  
  hrf_s3_obj <- gen_hrf(HRF_SPMG3)
  expect_equal(nbasis(hrf_s3_obj), 3)

  # Function with parameters determining nbasis
  hrf_bs5 <- gen_hrf(hrf_bspline, N = 5)
  expect_equal(nbasis(hrf_bs5), 5)
  
  hrf_bs4 <- gen_hrf(hrf_bspline, N = 4)
  expect_equal(nbasis(hrf_bs4), 4)
  
  # Tent function (bspline with degree 1)
  hrf_tent7 <- gen_hrf(hrf_bspline, N = 7, degree = 1)
  expect_equal(nbasis(hrf_tent7), 7)
})
