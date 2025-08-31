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
