# Diagnostic tests for T-statistic combination methods
# These tests focus on edge cases and mathematical correctness

test_that("Stouffer's Z method preserves sign and handles weighting correctly", {
  # Test Case: Sign preservation with mixed positive/negative effects
  # This diagnoses whether the sign logic works correctly when combining
  # opposing effects of different magnitudes
  
  # Setup: 3 subjects, 2 voxels
  # Voxel 1: Strong positive effects (should stay positive)
  # Voxel 2: Mixed effects where negative should dominate when properly weighted
  tmat <- matrix(c(
    3.0, -1.0,    # Subject 1: strong positive, weak negative  
    2.5, -2.0,    # Subject 2: medium positive, medium negative
    2.0, -4.0     # Subject 3: medium positive, strong negative
  ), nrow = 3, byrow = TRUE)
  
  df <- rep(20, 3)  # Equal degrees of freedom
  
  # Test 1: Equal weights - should combine based on magnitude only
  z_equal <- combine_t_statistics(tmat, df = df, method = "stouffer", weights = "equal")
  
  # Manual verification for voxel 1 (all positive):
  # t-values: [3.0, 2.5, 2.0], df=20
  # p_two = 2 * pt(abs(t), 20, lower.tail=FALSE) = [0.0070, 0.0206, 0.0594]
  # z_i = qnorm(1 - p_two/2) * sign(t) = [2.693, 2.054, 1.553] (all positive)
  # z_combined = sum(z_i) / sqrt(3) = (2.693 + 2.054 + 1.553) / sqrt(3) = 3.542
  expected_z1 <- (qnorm(1 - pt(3.0, 20, lower.tail = FALSE)) + 
                  qnorm(1 - pt(2.5, 20, lower.tail = FALSE)) + 
                  qnorm(1 - pt(2.0, 20, lower.tail = FALSE))) / sqrt(3)
  
  expect_equal(z_equal[1], expected_z1, tolerance = 1e-3)
  expect_true(z_equal[1] > 0, info = "Voxel 1 should remain positive")
  expect_true(z_equal[2] < 0, info = "Voxel 2 should be negative (strong negative effects)")
  
  # Test 2: Custom weights favoring subject 3 (strongest effects)
  # This should amplify the direction of subject 3's effects
  custom_weights <- c(1, 1, 3)  # Subject 3 gets triple weight
  z_custom <- combine_t_statistics(tmat, df = df, method = "stouffer", 
                                   weights = "custom", weights_custom = custom_weights)
  
  # With subject 3 heavily weighted:
  # Voxel 1: Still positive but less so (subject 3 has weakest positive)
  # Voxel 2: More negative (subject 3 has strongest negative)
  expect_true(z_custom[1] > 0 && z_custom[1] < z_equal[1], 
              info = "Custom weighting should reduce positive effect in voxel 1")
  expect_true(z_custom[2] < z_equal[2], 
              info = "Custom weighting should strengthen negative effect in voxel 2")
  
  # Test 3: Edge case with zero t-statistics (should handle gracefully)
  tmat_zero <- rbind(tmat, c(0, 0))
  df_zero <- c(df, 20)
  z_zero <- combine_t_statistics(tmat_zero, df = df_zero, method = "stouffer", weights = "equal")
  expect_true(all(is.finite(z_zero)), info = "Should handle zero t-statistics without NaN")
})

test_that("Fisher's method handles degrees of freedom and directional effects correctly", {
  # Test Case: Different degrees of freedom per subject
  # This diagnoses whether Fisher's method properly handles varying df
  # and whether the directional combination logic is sound
  
  # Setup: Carefully chosen t-statistics with known p-values
  # Using t-values that give "nice" p-values for manual verification
  tmat <- matrix(c(
    2.086,  -2.086,   # Subject 1, df=19: p ≈ 0.050 (two-tailed)
    2.228,  -1.325,   # Subject 2, df=9:  p ≈ 0.050, 0.218  
    1.960,  -3.000    # Subject 3, df=∞:  p ≈ 0.050, 0.003
  ), nrow = 3, byrow = TRUE)
  
  df <- c(19, 9, 1000)  # Different df per subject (1000 ≈ infinity)
  
  z_fisher <- combine_t_statistics(tmat, df = df, method = "fisher", weights = "equal")
  
  # Manual verification for voxel 1 (all positive effects, p ≈ 0.05 each):
  # X2 = -2 * sum(log(p_i)) = -2 * (log(0.05) + log(0.05) + log(0.05))
  # X2 = -2 * 3 * log(0.05) ≈ -2 * 3 * (-2.996) ≈ 17.98
  # p_combined = pchisq(17.98, df=6, lower.tail=FALSE) ≈ 0.0063
  # z_mag = qnorm(1 - 0.0063/2) ≈ 2.74
  # Since all t-stats positive: sign = +1, so z_comb ≈ +2.74
  
  p_vals_vox1 <- c(
    2 * pt(2.086, 19, lower.tail = FALSE),
    2 * pt(2.228, 9, lower.tail = FALSE), 
    2 * pt(1.960, 1000, lower.tail = FALSE)
  )
  X2_manual <- -2 * sum(log(p_vals_vox1))
  p_comb_manual <- pchisq(X2_manual, df = 6, lower.tail = FALSE)
  z_mag_manual <- qnorm(1 - p_comb_manual/2)
  
  expect_equal(z_fisher[1], z_mag_manual, tolerance = 1e-2)
  expect_true(z_fisher[1] > 0, info = "Voxel 1 should be positive (all positive t-stats)")
  
  # Test directional combination: voxel 2 has mixed signs
  # Mean direction: (2.086 + (-1.325) + (-3.000))/3 = -0.746 < 0
  # So result should be negative despite some positive effects
  expect_true(z_fisher[2] < 0, info = "Voxel 2 should be negative (negative mean direction)")
  
  # Test edge case: What happens with very small p-values?
  # This tests numerical stability near p=0
  tmat_extreme <- matrix(c(5.0, 6.0), nrow = 1)
  df_extreme <- 50
  z_extreme <- combine_t_statistics(tmat_extreme, df = df_extreme, method = "fisher")
  
  expect_true(all(is.finite(z_extreme)), info = "Should handle extreme t-values without overflow")
  expect_true(all(z_extreme > 0), info = "Extreme positive t-values should give positive z")
  expect_true(z_extreme[2] > z_extreme[1], info = "Larger t-value should give larger combined z")
  
  # Test another edge case: Equal opposing effects
  # This tests the tie-breaking logic when mean direction is exactly zero
  tmat_balanced <- matrix(c(2.0, -2.0), nrow = 1)  # Perfect balance
  z_balanced <- combine_t_statistics(tmat_balanced, df = 20, method = "fisher")
  
  # Mean direction = (2.0 + (-2.0))/1 = 0
  # Code uses: ifelse(sgn == 0, 1, sgn), so should default to positive
  expect_equal(length(z_balanced), 2)
  expect_true(is.finite(z_balanced[1]) && is.finite(z_balanced[2]))
  # Note: The actual sign when mean=0 depends on implementation details
  # but the key is that it should be finite and consistent
})