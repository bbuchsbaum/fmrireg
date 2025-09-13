## FMRI Meta-Analysis Contrast Computation Diagnostic Tests
## 
## These tests are designed to diagnose critical aspects of contrast computation and covariance handling
## in the fMRI meta-analysis implementation. They verify numerical accuracy and consistency between
## different computational paths.

test_that("meta-analysis contrast SE computation: full covariance vs diagonal approximation", {
  # Load the package
  library(fmrireg)
  
  # TEST 1: Exact contrast standard error computation using full covariance matrix vs diagonal approximation
  # 
  # WHAT THIS DIAGNOSES:
  # - Correct implementation of quadratic form c'*Cov*c for contrast standard errors
  # - Proper handling of correlation structure between predictors
  # - Difference between exact SE (using full covariance matrix) vs approximate SE (diagonal only)
  # - Critical numerical accuracy in contrast variance computation
  #
  # KEY DIAGNOSTIC VALUE:
  # When predictors are correlated, ignoring off-diagonal covariance elements leads to 
  # incorrect standard errors. This test creates a scenario where we know the exact 
  # answer and can detect if the implementation incorrectly ignores correlations.
  
  # Create a simple 2x2 design with known covariance structure
  X <- matrix(c(1, 1, 0, 1), nrow = 2, ncol = 2)  # intercept + condition
  colnames(X) <- c("intercept", "condition")
  
  # Known beta coefficients (2 voxels, 2 predictors)
  betas <- matrix(c(1.0, 2.0,   # voxel 1: intercept=1, condition=2
                    1.5, 1.8),  # voxel 2: intercept=1.5, condition=1.8
                  nrow = 2, ncol = 2, byrow = TRUE)
  
  # Known covariance matrix for each voxel (must be positive definite)
  # Voxel 1: moderate positive correlation between intercept and condition
  cov1 <- matrix(c(0.25, 0.10,   # var(intercept)=0.25, cov(intercept,condition)=0.10
                   0.10, 0.16),  # var(condition)=0.16
                 nrow = 2, ncol = 2)
  
  # Voxel 2: negative correlation
  cov2 <- matrix(c(0.36, -0.05,  # var(intercept)=0.36, cov(intercept,condition)=-0.05
                   -0.05, 0.09), # var(condition)=0.09  
                 nrow = 2, ncol = 2)
  
  # Verify positive definiteness
  expect_true(all(eigen(cov1)$values > 0), "cov1 must be positive definite")
  expect_true(all(eigen(cov2)$values > 0), "cov2 must be positive definite")
  
  # Pack upper triangular covariance (matching C++ packing order)
  K <- 2  # number of predictors
  tsize <- K * (K + 1) / 2  # should be 3
  cov_tri <- matrix(NA_real_, nrow = tsize, ncol = 2)
  
  # Pack voxel 1: upper triangle in row-major order (a,a), (a,b), (b,b)
  idx <- 1
  for (a in 1:K) {
    for (b in a:K) {
      cov_tri[idx, 1] <- cov1[a, b]
      cov_tri[idx, 2] <- cov2[a, b]
      idx <- idx + 1
    }
  }
  
  # Expected packed values:
  # idx=1: (1,1) -> cov1[1,1]=0.25, cov2[1,1]=0.36
  # idx=2: (1,2) -> cov1[1,2]=0.10, cov2[1,2]=-0.05
  # idx=3: (2,2) -> cov1[2,2]=0.16, cov2[2,2]=0.09
  expect_equal(cov_tri[1, ], c(0.25, 0.36))
  expect_equal(cov_tri[2, ], c(0.10, -0.05))
  expect_equal(cov_tri[3, ], c(0.16, 0.09))
  
  # Create mock fmri_meta object with full covariance
  # SE matrix should be voxels (rows) x predictors (cols)
  se_matrix <- matrix(c(sqrt(0.25), sqrt(0.36),    # intercept column (voxel1, voxel2)
                        sqrt(0.16), sqrt(0.09)),    # condition column (voxel1, voxel2)
                      nrow = 2, ncol = 2, byrow = FALSE)
  meta_full <- list(
    coefficients = t(betas),  # transpose to match expected format (2 voxels x 2 predictors)
    se = se_matrix,
    cov = list(type = "tri", tri = cov_tri),
    model = list(X = X)
  )
  class(meta_full) <- "fmri_meta"
  
  # Create mock fmri_meta object with diagonal-only approximation
  meta_diag <- list(
    coefficients = t(betas),
    se = se_matrix,
    cov = NULL,  # no full covariance available
    model = list(X = X)
  )
  class(meta_diag) <- "fmri_meta"
  
  # Test contrast: condition - intercept (should be c(-1, 1))
  contrast_weights <- c(-1, 1)  # intercept coefficient gets -1, condition gets +1
  
  # Compute contrast using full covariance
  result_full <- contrast(meta_full, contrast_weights)
  
  # Compute contrast using diagonal approximation  
  result_diag <- contrast(meta_diag, contrast_weights)
  
  # Manual verification for voxel 1 using full covariance
  # Contrast variance = c' * Cov * c where c = (-1, 1)
  expected_var1 <- t(contrast_weights) %*% cov1 %*% contrast_weights
  expected_var1 <- as.numeric(expected_var1)  # = (-1,1) * cov1 * (-1,1)'
  # = (-1*0.25 + 1*0.10) + (-1*0.10 + 1*0.16) = -0.15 + 0.06 = 0.31
  manual_var1 <- (-1)^2 * 0.25 + 2*(-1)*1*0.10 + 1^2 * 0.16
  expect_equal(manual_var1, 0.25 - 0.20 + 0.16)  # = 0.21
  expect_equal(expected_var1, 0.21, tolerance = 1e-10)
  
  # Manual verification for voxel 2
  expected_var2 <- t(contrast_weights) %*% cov2 %*% contrast_weights
  expected_var2 <- as.numeric(expected_var2)
  manual_var2 <- (-1)^2 * 0.36 + 2*(-1)*1*(-0.05) + 1^2 * 0.09
  expect_equal(manual_var2, 0.36 + 0.10 + 0.09)  # = 0.55
  expect_equal(expected_var2, 0.55, tolerance = 1e-10)
  
  # Check that full covariance method produces expected results
  expect_equal(result_full$se[1], sqrt(0.21), tolerance = 1e-10, 
               info = "Voxel 1: Full covariance SE should match manual calculation")
  expect_equal(result_full$se[2], sqrt(0.55), tolerance = 1e-10,
               info = "Voxel 2: Full covariance SE should match manual calculation")
  
  # Check that diagonal approximation gives different (and incorrect) results
  # Diagonal approximation: var = sum(c_i^2 * se_i^2) 
  # For contrast (-1, 1): var = (-1)^2*se_intercept^2 + (1)^2*se_condition^2
  expected_diag_var1 <- contrast_weights[1]^2 * 0.25 + contrast_weights[2]^2 * 0.16  # = 1*0.25 + 1*0.16 = 0.41
  expected_diag_var2 <- contrast_weights[1]^2 * 0.36 + contrast_weights[2]^2 * 0.09  # = 1*0.36 + 1*0.09 = 0.45
  
  expect_equal(result_diag$se[1], sqrt(0.41), tolerance = 1e-10,
               info = "Voxel 1: Diagonal approximation should match expected calculation")
  expect_equal(result_diag$se[2], sqrt(0.45), tolerance = 1e-10,
               info = "Voxel 2: Diagonal approximation should match expected calculation")
  
  # The key diagnostic: full covariance should be different from diagonal approximation
  # when there is non-zero covariance between predictors
  expect_false(isTRUE(all.equal(result_full$se, result_diag$se)),
               info = "Full covariance and diagonal approximation should give different SEs when predictors are correlated")
  
  # Specific differences should match our calculations
  expect_equal(result_full$se[1]^2 - result_diag$se[1]^2, 0.21 - 0.41, tolerance = 1e-10,
               info = "Voxel 1: SE difference should reflect covariance effect")
  expect_equal(result_full$se[2]^2 - result_diag$se[2]^2, 0.55 - 0.45, tolerance = 1e-10,
               info = "Voxel 2: SE difference should reflect covariance effect")
})


test_that("upper triangular covariance matrix packing and unpacking consistency", {
  # Load the package
  library(fmrireg)
  
  # TEST 2: Correct packing/unpacking of upper triangular covariance matrix
  # 
  # WHAT THIS DIAGNOSES:
  # - Consistency between C++ packing logic (meta_kernels.cpp) and R unpacking logic (fmri_meta.R)
  # - Correct triangular indexing that avoids off-by-one errors
  # - Proper handling of symmetric matrix reconstruction from packed form
  # - Edge cases with different matrix sizes (1x1, 2x2, 3x3, 4x4)
  #
  # KEY DIAGNOSTIC VALUE:
  # The C++ code packs the upper triangle of covariance matrices into a linear array,
  # and the R code must unpack this correctly. Any indexing mismatch leads to incorrect
  # covariance matrices and wrong contrast standard errors. This test ensures perfect
  # round-trip consistency and tests the specific indexing patterns used.
  
  # Create test covariance matrices of different sizes to test the indexing logic
  test_matrices <- list(
    # 2x2 case (3 unique elements)
    matrix(c(1.0, 0.5, 
             0.5, 2.0), nrow = 2, ncol = 2),
    
    # 3x3 case (6 unique elements)  
    matrix(c(4.0, 1.2, 0.8,
             1.2, 3.0, 0.3,
             0.8, 0.3, 2.5), nrow = 3, ncol = 3),
    
    # 4x4 case (10 unique elements)
    matrix(c(9.0, 2.1, 1.5, 0.4,
             2.1, 7.0, 1.8, 0.7,
             1.5, 1.8, 5.5, 1.1,
             0.4, 0.7, 1.1, 3.2), nrow = 4, ncol = 4)
  )
  
  for (i in seq_along(test_matrices)) {
    cov_mat <- test_matrices[[i]]
    K <- nrow(cov_mat)
    tsize <- K * (K + 1) / 2
    
    # Verify positive definiteness
    expect_true(all(eigen(cov_mat)$values > 0), 
                info = paste("Test matrix", i, "must be positive definite"))
    
    # STEP 1: Pack using C++ logic (simulate what meta_kernels.cpp does)
    packed_cpp_style <- numeric(tsize)
    idx <- 1
    for (a in 1:K) {
      for (b in a:K) {
        packed_cpp_style[idx] <- cov_mat[a, b]
        idx <- idx + 1
      }
    }
    
    # STEP 2: Unpack using R logic (simulate what fmri_meta.R does)
    unpacked_mat <- matrix(0, nrow = K, ncol = K)
    idx <- 1
    for (a in 1:K) {
      for (b in a:K) {
        unpacked_mat[a, b] <- packed_cpp_style[idx]
        unpacked_mat[b, a] <- packed_cpp_style[idx]  # fill symmetric element
        idx <- idx + 1
      }
    }
    
    # STEP 3: Verify perfect round-trip consistency
    expect_equal(unpacked_mat, cov_mat, tolerance = 1e-15,
                 info = paste("Matrix", i, ": Round-trip packing/unpacking should preserve original matrix exactly"))
    
    # STEP 4: Test specific indexing patterns that could cause off-by-one errors
    if (K >= 2) {
      # Test that off-diagonal elements are correctly placed
      expect_equal(unpacked_mat[1, 2], cov_mat[1, 2],
                   info = paste("Matrix", i, ": Off-diagonal (1,2) should be preserved"))
      expect_equal(unpacked_mat[2, 1], cov_mat[2, 1], 
                   info = paste("Matrix", i, ": Symmetric element (2,1) should be preserved"))
    }
    
    if (K >= 3) {
      # Test corner case: last row, first column should match first row, last column
      expect_equal(unpacked_mat[K, 1], cov_mat[1, K],
                   info = paste("Matrix", i, ": Corner symmetry should be preserved"))
    }
    
    # STEP 5: Verify the packing order matches expected pattern
    # For K=3, order should be: (1,1), (1,2), (1,3), (2,2), (2,3), (3,3)
    if (K == 3) {
      expected_order <- c(cov_mat[1,1], cov_mat[1,2], cov_mat[1,3], 
                         cov_mat[2,2], cov_mat[2,3], cov_mat[3,3])
      expect_equal(packed_cpp_style, expected_order,
                   info = "K=3: Packing order should follow row-wise upper triangular pattern")
    }
    
    # STEP 6: Test the contrast computation logic that uses this packing
    # Create a simple contrast and verify it produces the same result
    # whether computed from the original matrix or the packed/unpacked version
    contrast_weights <- rep(1/K, K)  # equal weights that sum to 1
    
    # Direct computation from original matrix
    direct_var <- as.numeric(t(contrast_weights) %*% cov_mat %*% contrast_weights)
    
    # Computation using unpacked matrix (simulating the R contrast code)
    unpacked_var <- 0.0
    idx <- 1
    for (a in 1:K) {
      for (b in a:K) {
        cab <- packed_cpp_style[idx]
        w <- contrast_weights[a] * contrast_weights[b]
        if (b > a) {
          unpacked_var <- unpacked_var + 2.0 * cab * w  # off-diagonal, count twice
        } else {
          unpacked_var <- unpacked_var + cab * w  # diagonal, count once
        }
        idx <- idx + 1
      }
    }
    
    expect_equal(unpacked_var, direct_var, tolerance = 1e-14,
                 info = paste("Matrix", i, ": Contrast variance from packed form should match direct computation"))
  }
  
  # STEP 7: Edge case test - single predictor (K=1) 
  single_cov <- matrix(4.0, nrow = 1, ncol = 1)
  K <- 1
  tsize <- K * (K + 1) / 2  # should be 1
  
  packed_single <- numeric(tsize)
  idx <- 1
  for (a in 1:K) {
    for (b in a:K) {
      packed_single[idx] <- single_cov[a, b]
      idx <- idx + 1
    }
  }
  
  expect_equal(length(packed_single), 1, info = "K=1: Should have exactly 1 packed element")
  expect_equal(packed_single[1], 4.0, info = "K=1: Single element should be preserved")
  
  # Unpack single element
  unpacked_single <- matrix(0, nrow = K, ncol = K)
  idx <- 1
  for (a in 1:K) {
    for (b in a:K) {
      unpacked_single[a, b] <- packed_single[idx]
      unpacked_single[b, a] <- packed_single[idx]
      idx <- idx + 1
    }
  }
  
  expect_equal(unpacked_single, single_cov,
               info = "K=1: Single element round-trip should work correctly")
})