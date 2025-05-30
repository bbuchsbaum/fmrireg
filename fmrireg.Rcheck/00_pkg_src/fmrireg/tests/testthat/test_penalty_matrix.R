library(testthat)
library(fmrireg)

# Test helper functions ----

test_that("first_difference creates correct difference matrix", {
  # Test with d = 1 (edge case)
  D1 <- fmrireg:::first_difference(1)
  expect_equal(nrow(D1), 0)
  expect_equal(ncol(D1), 1)
  
  # Test with d = 2
  D2 <- fmrireg:::first_difference(2)
  expected2 <- matrix(c(-1, 1), nrow = 1)
  expect_equal(D2, expected2)
  
  # Test with d = 3
  D3 <- fmrireg:::first_difference(3)
  expected3 <- matrix(c(-1, 1, 0, 0, -1, 1), nrow = 2, byrow = TRUE)
  expect_equal(D3, expected3)
  
  # Test with d = 4
  D4 <- fmrireg:::first_difference(4)
  expected4 <- matrix(c(-1, 1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1), nrow = 3, byrow = TRUE)
  expect_equal(D4, expected4)
})

test_that("second_difference creates correct difference matrix", {
  # Test with d = 1 (edge case)
  D1 <- fmrireg:::second_difference(1)
  expect_equal(nrow(D1), 0)
  expect_equal(ncol(D1), 1)
  
  # Test with d = 2 (edge case)
  D2 <- fmrireg:::second_difference(2)
  expect_equal(nrow(D2), 0)
  expect_equal(ncol(D2), 2)
  
  # Test with d = 3
  D3 <- fmrireg:::second_difference(3)
  expected3 <- matrix(c(1, -2, 1), nrow = 1)
  expect_equal(D3, expected3)
  
  # Test with d = 4
  D4 <- fmrireg:::second_difference(4)
  expected4 <- matrix(c(1, -2, 1, 0, 0, 1, -2, 1), nrow = 2, byrow = TRUE)
  expect_equal(D4, expected4)
  
  # Test with d = 5
  D5 <- fmrireg:::second_difference(5)
  expected5 <- matrix(c(1, -2, 1, 0, 0, 0, 1, -2, 1, 0, 0, 0, 1, -2, 1), nrow = 3, byrow = TRUE)
  expect_equal(D5, expected5)
})

# Test penalty matrix for different HRF types ----

test_that("penalty_matrix works for HRF_SPMG1 (canonical)", {
  R <- penalty_matrix(HRF_SPMG1)
  
  # Should be 1x1 identity matrix
  expect_equal(dim(R), c(1, 1))
  expect_equal(R[1,1], 1)
  expect_true(isSymmetric(R))
})

test_that("penalty_matrix works for HRF_SPMG2 (canonical + temporal derivative)", {
  R <- penalty_matrix(HRF_SPMG2)
  
  # Should be 2x2 diagonal matrix with differential shrinkage
  expect_equal(dim(R), c(2, 2))
  expect_true(isSymmetric(R))
  expect_true(all(diag(R) > 0))
  expect_equal(R[1,1], 1)  # Canonical not penalized
  expect_true(R[2,2] > R[1,1])  # Derivative penalized more
  expect_equal(R[1,2], 0)  # Off-diagonal should be zero
  expect_equal(R[2,1], 0)
})

test_that("penalty_matrix works for HRF_SPMG3 (canonical + both derivatives)", {
  R <- penalty_matrix(HRF_SPMG3)
  
  # Should be 3x3 diagonal matrix
  expect_equal(dim(R), c(3, 3))
  expect_true(isSymmetric(R))
  expect_true(all(diag(R) > 0))
  expect_equal(R[1,1], 1)  # Canonical not penalized
  expect_true(R[2,2] > R[1,1])  # Temporal derivative penalized
  expect_true(R[3,3] > R[1,1])  # Dispersion derivative penalized
  
  # Off-diagonals should be zero
  expect_equal(R[1,2], 0)
  expect_equal(R[1,3], 0)
  expect_equal(R[2,3], 0)
})

test_that("penalty_matrix works for HRF_FIR", {
  R <- penalty_matrix(HRF_FIR)
  
  # Should be symmetric and positive semi-definite
  expect_true(isSymmetric(R))
  expect_equal(dim(R), c(nbasis(HRF_FIR), nbasis(HRF_FIR)))
  
  # Check that it's a valid penalty matrix (positive semi-definite)
  eigenvals <- eigen(R, only.values = TRUE)$values
  expect_true(all(eigenvals >= -1e-10))  # Allow for numerical precision
  
  # Should have banded structure (non-zero near diagonal)
  # For second-order differences, should be pentadiagonal
  expect_true(all(R[abs(row(R) - col(R)) > 2] == 0))
})

test_that("penalty_matrix works for HRF_BSPLINE", {
  R <- penalty_matrix(HRF_BSPLINE)
  
  # Should be symmetric and positive semi-definite
  expect_true(isSymmetric(R))
  expect_equal(dim(R), c(nbasis(HRF_BSPLINE), nbasis(HRF_BSPLINE)))
  
  # Check that it's a valid penalty matrix
  eigenvals <- eigen(R, only.values = TRUE)$values
  expect_true(all(eigenvals >= -1e-10))
  
  # Should have banded structure like FIR
  expect_true(all(R[abs(row(R) - col(R)) > 2] == 0))
})

# Test penalty matrix parameters ----

test_that("penalty_matrix respects order parameter for FIR", {
  # Test first-order differences
  R1 <- penalty_matrix(HRF_FIR, order = 1)
  expect_true(isSymmetric(R1))
  
  # Test second-order differences  
  R2 <- penalty_matrix(HRF_FIR, order = 2)
  expect_true(isSymmetric(R2))
  
  # They should be different
  expect_false(identical(R1, R2))
  
  # First-order should have tridiagonal structure
  expect_true(all(R1[abs(row(R1) - col(R1)) > 1] == 0))
  
  # Second-order should have pentadiagonal structure
  expect_true(all(R2[abs(row(R2) - col(R2)) > 2] == 0))
})

test_that("penalty_matrix respects scale parameter", {
  # Test with scaling
  R_scaled <- penalty_matrix(HRF_FIR, scale = TRUE)
  
  # Test without scaling
  R_unscaled <- penalty_matrix(HRF_FIR, scale = FALSE)
  
  # They should be different
  expect_false(identical(R_scaled, R_unscaled))
  
  # Scaled version should have diagonal elements closer to 1
  if (any(diag(R_scaled) > 0)) {
    expect_true(abs(mean(diag(R_scaled)[diag(R_scaled) > 0]) - 1) < 0.1)
  }
})

test_that("penalty_matrix respects shrink_deriv parameter for SPMG", {
  # Test with default shrinkage
  R_default <- penalty_matrix(HRF_SPMG2)
  
  # Test with custom shrinkage
  R_custom <- penalty_matrix(HRF_SPMG2, shrink_deriv = 10)
  
  # They should be different
  expect_false(identical(R_default, R_custom))
  
  # Custom should have larger penalty on derivative
  expect_true(R_custom[2,2] > R_default[2,2])
  expect_equal(R_custom[1,1], R_default[1,1])  # Canonical unchanged
})

test_that("penalty_matrix respects shrink_disp parameter for SPMG3", {
  # Test with default (shrink_disp = shrink_deriv)
  R_default <- penalty_matrix(HRF_SPMG3, shrink_deriv = 4)
  
  # Test with custom dispersion shrinkage
  R_custom <- penalty_matrix(HRF_SPMG3, shrink_deriv = 4, shrink_disp = 8)
  
  # They should be different
  expect_false(identical(R_default, R_custom))
  
  # Custom should have different penalty on dispersion derivative
  expect_true(R_custom[3,3] != R_default[3,3])
  expect_equal(R_custom[1,1], R_default[1,1])  # Canonical unchanged
  expect_equal(R_custom[2,2], R_default[2,2])  # Temporal derivative unchanged
})

# Test edge cases ----

test_that("penalty_matrix handles small basis sets gracefully", {
  # Save original nbasis function
  original_nbasis <- getFromNamespace("nbasis.HRF", "fmrireg")
  
  # Test with 1 basis function
  mock_hrf_1 <- HRF_SPMG1
  attr(mock_hrf_1, "name") <- "fir"
  
  # Mock nbasis to return 1
  mock_nbasis_1 <- function(x, ...) 1L
  assignInNamespace("nbasis.HRF", mock_nbasis_1, "fmrireg")
  
  R <- penalty_matrix(mock_hrf_1)
  expect_equal(dim(R), c(1, 1))
  expect_equal(R[1,1], 1)
  
  # Test with 2 basis functions for second-order differences
  mock_hrf_2 <- HRF_SPMG1
  attr(mock_hrf_2, "name") <- "fir"
  
  # Mock nbasis to return 2
  mock_nbasis_2 <- function(x, ...) 2L
  assignInNamespace("nbasis.HRF", mock_nbasis_2, "fmrireg")
  
  # Should fall back to identity for second-order with only 2 basis functions
  R <- penalty_matrix(mock_hrf_2, order = 2)
  expect_equal(R, diag(2))
  
  # Restore original function
  assignInNamespace("nbasis.HRF", original_nbasis, "fmrireg")
})

test_that("penalty_matrix handles unknown HRF types", {
  # Create a mock HRF with unknown name
  mock_hrf <- HRF_SPMG1
  attr(mock_hrf, "name") <- "unknown_type"
  
  R <- penalty_matrix(mock_hrf)
  
  # Should fall back to identity matrix
  expect_equal(R, diag(nbasis(mock_hrf)))
})

test_that("penalty_matrix handles HRF without name attribute", {
  # Create a mock HRF without name
  mock_hrf <- HRF_SPMG1
  attr(mock_hrf, "name") <- NULL
  
  R <- penalty_matrix(mock_hrf)
  
  # Should fall back to identity matrix
  expect_equal(R, diag(nbasis(mock_hrf)))
})

test_that("penalty_matrix validates order parameter", {
  # Test that invalid order values throw errors
  expect_error(penalty_matrix(HRF_FIR, order = 3))
  expect_error(penalty_matrix(HRF_FIR, order = 0))
  expect_error(penalty_matrix(HRF_FIR, order = -1))
  
  # Test that valid orders work
  expect_no_error(penalty_matrix(HRF_FIR, order = 1))
  expect_no_error(penalty_matrix(HRF_FIR, order = 2))
})

# Test mathematical properties ----

test_that("penalty matrices are positive semi-definite", {
  hrfs <- list(
    HRF_SPMG1, HRF_SPMG2, HRF_SPMG3,
    HRF_FIR, HRF_BSPLINE
  )
  
  for (hrf in hrfs) {
    R <- penalty_matrix(hrf)
    eigenvals <- eigen(R, only.values = TRUE)$values
    expect_true(all(eigenvals >= -1e-10), 
                info = paste("HRF:", attr(hrf, "name")))
  }
})

test_that("penalty matrices are symmetric", {
  hrfs <- list(
    HRF_SPMG1, HRF_SPMG2, HRF_SPMG3,
    HRF_FIR, HRF_BSPLINE
  )
  
  for (hrf in hrfs) {
    R <- penalty_matrix(hrf)
    expect_true(isSymmetric(R), 
                info = paste("HRF:", attr(hrf, "name")))
  }
})

test_that("penalty matrices have correct dimensions", {
  hrfs <- list(
    HRF_SPMG1, HRF_SPMG2, HRF_SPMG3,
    HRF_FIR, HRF_BSPLINE
  )
  
  for (hrf in hrfs) {
    R <- penalty_matrix(hrf)
    d <- nbasis(hrf)
    expect_equal(dim(R), c(d, d), 
                info = paste("HRF:", attr(hrf, "name")))
  }
})

# Test integration with other functions ----

test_that("penalty_matrix works with hrfspec objects", {
  # Test that the method exists in the methods list
  methods_list <- methods("penalty_matrix")
  expect_true("penalty_matrix.hrfspec" %in% methods_list)
})

# Test specific penalty structures ----

test_that("FIR penalty matrix has expected structure", {
  # Create a simple FIR-like HRF for testing
  mock_fir <- HRF_FIR
  
  R <- penalty_matrix(mock_fir, order = 2, scale = FALSE)
  d <- nbasis(mock_fir)
  
  # Should be banded with specific pattern
  # Main diagonal should be positive
  expect_true(all(diag(R) > 0))
  
  # Should have the characteristic pattern of D^T D where D is second difference
  # For second differences: main diagonal has 1, 5, 6, 6, ..., 6, 5, 1
  # Off-diagonals have specific patterns
  if (d >= 3) {
    # Check that it has the right sparsity pattern
    expect_true(all(R[abs(row(R) - col(R)) > 2] == 0))
  }
})

test_that("SPMG penalty matrices have diagonal structure", {
  # Test SPMG2
  R2 <- penalty_matrix(HRF_SPMG2, shrink_deriv = 4)
  expect_true(all(R2[row(R2) != col(R2)] == 0))  # Off-diagonal zeros
  expect_equal(R2[1,1], 1)
  expect_equal(R2[2,2], 4)
  
  # Test SPMG3 - verify it works with custom parameters
  R3 <- penalty_matrix(HRF_SPMG3, shrink_deriv = 4, shrink_disp = 6)
  expect_equal(dim(R3), c(3, 3))  # Should be 3x3
  expect_true(all(R3[row(R3) != col(R3)] == 0))  # Off-diagonal zeros
  expect_equal(R3[1,1], 1)
  expect_equal(R3[2,2], 4)
  expect_equal(R3[3,3], 6)
  
  # Test SPMG3 with default parameters (shrink_disp = shrink_deriv)
  R3_default <- penalty_matrix(HRF_SPMG3, shrink_deriv = 4)
  expect_equal(dim(R3_default), c(3, 3))
  expect_equal(R3_default[1,1], 1)
  expect_equal(R3_default[2,2], 4)
  expect_equal(R3_default[3,3], 4)  # Should equal shrink_deriv when shrink_disp is NULL
})

# Test name-based dispatch ----

test_that("name-based dispatch works correctly", {
  # Test different name variations
  test_names <- list(
    c("fir", "FIR"),
    c("bspline", "BSPLINE", "bs", "BS"),
    c("spmg1", "SPMG1", "spmg2", "SPMG2"),
    c("fourier", "FOURIER", "cosine", "COSINE")
  )
  
  for (names_group in test_names) {
    # Create mock HRFs with different name cases
    for (name in names_group) {
      mock_hrf <- HRF_SPMG1
      attr(mock_hrf, "name") <- name
      
      # Should not error
      R <- penalty_matrix(mock_hrf)
      expect_true(is.matrix(R))
      expect_true(isSymmetric(R))
    }
  }
}) 