test_that("fmri_meta_fit FE recovers simple means and shapes", {
  skip_on_cran()

  # 4 subjects, 2 features, intercept-only design
  Y <- matrix(c(1,1,1,1, 0,2,0,2), nrow = 4, byrow = FALSE)  # subjects x features
  V <- matrix(1, nrow = 4, ncol = 2)                         # unit variances
  X <- matrix(1, nrow = 4, ncol = 1)                         # intercept

  fit <- fmrireg:::fmri_meta_fit(Y = Y, V = V, X = X, method = "fe", robust = "none")

  # beta, se, z have predictors x features
  expect_equal(dim(fit$beta), c(1, 2))
  expect_equal(dim(fit$se), c(1, 2))
  expect_true(all(is.finite(fit$beta)))
  expect_true(all(is.finite(fit$se)))

  # Weighted mean with unit variances equals simple mean
  # Feature 1: mean = 1; Feature 2: mean = 1
  expect_equal(as.numeric(fit$beta), c(1, 1), tolerance = 1e-8)

  # SE for intercept under FE with unit variances: sqrt(1/sum(w)) = 1/sqrt(4) = 0.5
  expect_equal(as.numeric(fit$se), c(0.5, 0.5), tolerance = 1e-8)

  # FE uses tau2 = 0 and reports Q/I2/df
  expect_true(all(fit$tau2 == 0))
  expect_true(all(fit$df == (nrow(Y) - ncol(X))))
  expect_type(fit$ok, "logical")
  expect_identical(fit$method, "fe")
  expect_identical(fit$robust, "none")
})

test_that("fmri_meta_fit DL/PM estimate tau2 on heterogeneous feature", {
  skip_on_cran()

  # Same setup as above
  Y <- matrix(c(1,1,1,1, 0,2,0,2), nrow = 4, byrow = FALSE)
  V <- matrix(1, nrow = 4, ncol = 2)
  X <- matrix(1, nrow = 4, ncol = 1)

  # For feature 2, FE Q0 = 4, df = 3 -> expected tau2 ≈ 1/3 for DL/PM
  target_tau2 <- 1/3

  fit_dl <- fmrireg:::fmri_meta_fit(Y = Y, V = V, X = X, method = "dl", robust = "none")
  fit_pm <- fmrireg:::fmri_meta_fit(Y = Y, V = V, X = X, method = "pm", robust = "none")

  # Feature 1 perfectly homogeneous => tau2 ~ 0
  expect_equal(fit_dl$tau2[1], 0, tolerance = 1e-10)
  expect_equal(fit_pm$tau2[1], 0, tolerance = 1e-10)

  # Feature 2 heterogeneous => tau2 close to 1/3
  expect_equal(fit_dl$tau2[2], target_tau2, tolerance = 1e-8)
  expect_equal(fit_pm$tau2[2], target_tau2, tolerance = 1e-6)

  # RE increases SE relative to FE for heterogeneous feature
  fit_fe <- fmrireg:::fmri_meta_fit(Y = Y, V = V, X = X, method = "fe", robust = "none")
  expect_gt(fit_dl$se[1,2], fit_fe$se[1,2])
  expect_gt(fit_pm$se[1,2], fit_fe$se[1,2])
})

test_that("fmri_meta_fit robust huber downweights outlier", {
  skip_on_cran()

  # Strong outlier should pull FE estimate; robust should reduce impact
  Y <- matrix(c(0,0,0,10), ncol = 1)  # subjects x 1 feature
  V <- matrix(1, nrow = 4, ncol = 1)
  X <- matrix(1, nrow = 4, ncol = 1)

  fe <- fmrireg:::fmri_meta_fit(Y = Y, V = V, X = X, method = "fe", robust = "none")
  rb <- fmrireg:::fmri_meta_fit(Y = Y, V = V, X = X, method = "fe", robust = "huber")

  # FE mean = 2.5; robust should be strictly smaller due to downweighting
  expect_equal(as.numeric(fe$beta), 2.5, tolerance = 1e-12)
  expect_lt(as.numeric(rb$beta), 2.5)
})

test_that("fmri_meta_fit validates shapes and handles invalid variances", {
  skip_on_cran()

  Y <- matrix(rnorm(6), nrow = 3, ncol = 2)
  V <- matrix(1, nrow = 3, ncol = 2)
  X <- cbind(1, rnorm(3))

  # Shape mismatches
  expect_error(fmrireg:::fmri_meta_fit(Y = Y, V = V[1:2,], X = X))
  expect_error(fmrireg:::fmri_meta_fit(Y = Y, V = V, X = X[1:2,]))

  # Non-positive variances trigger a warning and are treated as missing
  V_bad <- V; V_bad[1,1] <- 0
  expect_warning(fmrireg:::fmri_meta_fit(Y = Y, V = V_bad, X = X), "Non-positive variances")
})

test_that("tau2 estimation algorithms DL vs PM mathematical correctness", {
  skip_on_cran()
  
  # DIAGNOSTIC TEST 1: Verify DL and PM tau2 estimation algorithms
  # 
  # This test uses a carefully constructed 3-study intercept-only meta-analysis
  # where both DL and PM should converge to the same analytical solution.
  #
  # Setup: 3 studies with effects [1, 3, 5] and variances [4, 1, 4]
  # 
  # Hand calculations:
  # - FE weights: w = 1/v = [0.25, 1, 0.25]
  # - FE pooled estimate: sum(w*y)/sum(w) = (0.25*1 + 1*3 + 0.25*5)/(0.25+1+0.25) = 4.5/1.5 = 3
  # - FE residuals: r = y - 3 = [-2, 0, 2] 
  # - Q0 = sum(w * r²) = 0.25*4 + 1*0 + 0.25*4 = 2
  # - df = 3 - 1 = 2
  # - For DL: C = sum(w) - sum(w²*t_i) where t_i = x_i'(X'WX)^(-1)x_i
  #   With intercept-only: (X'WX)^(-1) = 1/sum(w) = 1/1.5 = 2/3
  #   So t_i = 1 * (2/3) * 1 = 2/3 for all i
  #   C = 1.5 - (0.25² + 1² + 0.25²) * (2/3) = 1.5 - 1.125*(2/3) = 1.5 - 0.75 = 0.75
  # - DL tau2 = (Q0 - df)/C = (2 - 2)/0.75 = 0
  # - PM should also give tau2 = 0 since Q(0) = Q0 = 2 = df
  
  Y <- matrix(c(1, 3, 5), ncol = 1)  # 3 studies, 1 feature
  V <- matrix(c(4, 1, 4), ncol = 1)  # Different variances
  X <- matrix(1, nrow = 3, ncol = 1) # Intercept only
  
  fit_dl <- fmrireg:::fmri_meta_fit(Y = Y, V = V, X = X, method = "dl", robust = "none")
  fit_pm <- fmrireg:::fmri_meta_fit(Y = Y, V = V, X = X, method = "pm", robust = "none")
  
  # Both methods should give tau2 = 0 (exactly, within numerical precision)
  expect_equal(fit_dl$tau2[1], 0, tolerance = 1e-12)
  expect_equal(fit_pm$tau2[1], 0, tolerance = 1e-12)
  
  # Both should give same beta estimate = 3
  expect_equal(as.numeric(fit_dl$beta), 3, tolerance = 1e-12)
  expect_equal(as.numeric(fit_pm$beta), 3, tolerance = 1e-12)
  
  # Q_fe should equal 2 as calculated
  expect_equal(fit_dl$Q_fe[1], 2, tolerance = 1e-12)
  expect_equal(fit_pm$Q_fe[1], 2, tolerance = 1e-12)
  
  # Now test case where tau2 > 0: increase Q0 by making effects more heterogeneous
  Y2 <- matrix(c(0, 3, 6), ncol = 1)  # More spread: [-3, 0, 3] residuals from mean=3
  fit_dl2 <- fmrireg:::fmri_meta_fit(Y = Y2, V = V, X = X, method = "dl", robust = "none")
  fit_pm2 <- fmrireg:::fmri_meta_fit(Y = Y2, V = V, X = X, method = "pm", robust = "none")
  
  # New Q0 = 0.25*9 + 1*0 + 0.25*9 = 4.5, so DL tau2 = (4.5 - 2)/0.75 = 3.33...
  expected_tau2_dl <- (4.5 - 2) / 0.75
  expect_equal(fit_dl2$tau2[1], expected_tau2_dl, tolerance = 1e-10)
  
  # PM solves Q(tau2) = df = 2. With tau2=5: Q(5) = 2 exactly 
  # (verified manually: weights 1/9, 1/6, 1/9 give weighted mean 3, residuals [-3,0,3], Q=2)
  expected_tau2_pm <- 5.0
  expect_equal(fit_pm2$tau2[1], expected_tau2_pm, tolerance = 1e-5)
  
  # This demonstrates that DL and PM can give different tau2 estimates - this is expected!
  expect_true(abs(fit_dl2$tau2[1] - fit_pm2$tau2[1]) > 1)
})

test_that("missing and infinite value filtering logic robustness", {
  skip_on_cran()
  
  # DIAGNOSTIC TEST 2: Test robustness of C++ filtering logic in fit_one()
  #
  # The C++ code uses: intersect(intersect(finite_y, finite_v), positive_v)
  # This test verifies edge cases in this filtering logic
  
  # Case 1: Mixed missing/infinite values should be properly filtered
  Y <- matrix(c(1, NA, 3, Inf, -Inf, 2), ncol = 1)
  V <- matrix(c(1, 1, NA, 1, 1, 0), ncol = 1)  # Note: V[6] = 0 (non-positive)
  X <- matrix(1, nrow = 6, ncol = 1)
  
  # Only rows 1 and 4 should be valid: finite Y, finite positive V
  # Row 1: Y=1, V=1 ✓
  # Row 2: Y=NA ✗  
  # Row 3: Y=3, V=NA ✗
  # Row 4: Y=Inf ✗
  # Row 5: Y=-Inf ✗  
  # Row 6: Y=2, V=0 ✗ (non-positive)
  # Wait, let me recalculate: only row 1 is fully valid
  Y <- matrix(c(1, NA, 3, 4, 5, 2), ncol = 1)  # Make row 4 finite
  V <- matrix(c(1, 1, NA, 1, 1, 0), ncol = 1)
  
  # Valid rows: 1 (Y=1,V=1), 4 (Y=4,V=1), 5 (Y=5,V=1)
  # This gives 3 valid observations with 1 parameter -> df = 2
  
  fit <- fmrireg:::fmri_meta_fit(Y = Y, V = V, X = X, method = "fe", robust = "none")
  expect_true(fit$ok[1])  # Should succeed
  expect_equal(fit$df[1], 2)  # 3 valid obs - 1 parameter
  
  # Expected beta = weighted mean of [1,4,5] with weights [1,1,1] = 10/3
  expect_equal(as.numeric(fit$beta), 10/3, tolerance = 1e-12)
  
  # Case 2: Boundary condition - exactly p valid observations (should fail)
  Y_boundary <- matrix(c(1, NA, NA), ncol = 1)
  V_boundary <- matrix(c(1, 1, 1), ncol = 1)
  X_boundary <- matrix(1, nrow = 3, ncol = 1)
  
  # Only 1 valid observation, but need > p observations (p=1 for intercept)
  fit_boundary <- fmrireg:::fmri_meta_fit(Y = Y_boundary, V = V_boundary, X = X_boundary, 
                                         method = "fe", robust = "none")
  expect_false(fit_boundary$ok[1])  # Should fail
  expect_true(is.na(fit_boundary$beta[1]))
  
  # Case 3: All infinite variances (edge case)
  Y_inf <- matrix(c(1, 2, 3), ncol = 1)
  V_inf <- matrix(c(Inf, Inf, Inf), ncol = 1)
  
  fit_inf <- fmrireg:::fmri_meta_fit(Y = Y_inf, V = V_inf, X = X_boundary, 
                                    method = "fe", robust = "none")
  expect_false(fit_inf$ok[1])  # Should fail - no positive finite variances
  
  # Case 4: Recovery with sufficient valid data after extensive filtering
  Y_recover <- matrix(c(NA, Inf, -Inf, 1, 2, 3, 4), ncol = 1)
  V_recover <- matrix(c(1, NA, 0, 1, 1, 1, 1), ncol = 1)
  X_recover <- matrix(1, nrow = 7, ncol = 1)
  
  # Valid rows: 4,5,6,7 (4 observations > 1 parameter)
  fit_recover <- fmrireg:::fmri_meta_fit(Y = Y_recover, V = V_recover, X = X_recover, 
                                        method = "dl", robust = "none")
  expect_true(fit_recover$ok[1])
  expect_equal(fit_recover$df[1], 3)  # 4 valid - 1 parameter  
  expect_equal(as.numeric(fit_recover$beta), mean(c(1,2,3,4)), tolerance = 1e-12)
})
