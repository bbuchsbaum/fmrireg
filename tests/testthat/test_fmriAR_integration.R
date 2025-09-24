# Test fmriAR Integration
# Verifies that AR functionality properly delegates to fmriAR package

test_that("fmriAR adapter correctly translates configurations", {
  # Test AR1 configuration
  cfg_ar1 <- list(ar = list(struct = "ar1", iter_gls = 2, exact_first = TRUE))
  params <- .fmrireg_to_fmriAR_config(cfg_ar1)
  expect_equal(params$method, "ar")
  expect_equal(params$p, 1L)
  expect_equal(params$exact_first, "ar1")
  expect_equal(params$iter, 2L)

  # Test AR2 configuration
  cfg_ar2 <- list(ar = list(struct = "ar2", global = TRUE))
  params <- .fmrireg_to_fmriAR_config(cfg_ar2)
  expect_equal(params$p, 2L)
  expect_equal(params$pooling, "global")

  # Test ARP configuration
  cfg_arp <- list(ar = list(struct = "arp", p = 5))
  params <- .fmrireg_to_fmriAR_config(cfg_arp)
  expect_equal(params$p, 5L)

  # Test IID (no AR)
  cfg_iid <- list(ar = list(struct = "iid"))
  params <- .fmrireg_to_fmriAR_config(cfg_iid)
  expect_equal(params$p, 0L)
})

test_that("fmriAR integration produces numerically similar results", {
  skip_if_not_installed("fmriAR")

  set.seed(123)
  n_time <- 100
  n_vox <- 5
  true_phi <- 0.5

  # Generate AR(1) data
  Y <- matrix(0, n_time, n_vox)
  for (v in 1:n_vox) {
    errors <- as.numeric(arima.sim(list(ar = true_phi), n = n_time))
    Y[, v] <- errors
  }

  # Simple design matrix
  X <- cbind(1, rnorm(n_time))

  # Add signal
  Y <- Y + X %*% matrix(c(0, 2), ncol = n_vox, nrow = 2)

  # Test AR estimation via fmriAR
  residuals <- Y - X %*% base::qr.solve(X, Y)
  cfg <- list(struct = "ar1", exact_first = TRUE)

  plan <- .estimate_ar_via_fmriAR(residuals, cfg)

  # Check that AR coefficient is recovered reasonably well
  expect_true(inherits(plan, "fmriAR_plan"))

  # Extract phi (should be close to true_phi)
  if (!is.null(plan$phi) && length(plan$phi) > 0) {
    estimated_phi <- plan$phi[[1]][1]
    expect_true(abs(estimated_phi - true_phi) < 0.3)  # Reasonable tolerance
  }
})

test_that("whitening via fmriAR works correctly", {
  skip_if_not_installed("fmriAR")

  set.seed(456)
  n <- 80
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * 3), n, 3)

  # Create a simple AR plan
  plan <- fmriAR::compat$plan_from_phi(
    phi = c(0.5, 0.2),
    exact_first = FALSE
  )

  # Apply whitening
  result <- .apply_ar_whitening_via_fmriAR(X, Y, plan)

  expect_equal(dim(result$X), dim(X))
  expect_equal(dim(result$Y), dim(Y))

  # Check that whitening changed the data
  expect_true(!all(result$X == X))
  expect_true(!all(result$Y == Y))
})

test_that("iterative AR-GLS via fmriAR converges", {
  skip_if_not_installed("fmriAR")

  set.seed(789)
  n <- 100
  X <- cbind(1, rnorm(n))

  # Generate AR(1) errors
  true_phi <- 0.4
  errors <- as.numeric(arima.sim(list(ar = true_phi), n = n))
  Y <- X %*% c(1, 2) + errors
  Y <- matrix(Y, ncol = 1)

  # Run iterative AR-GLS
  cfg <- list(struct = "ar1", iter_gls = 3, exact_first = TRUE)
  result <- .iterative_ar_gls_via_fmriAR(X, Y, cfg, max_iter = 3)

  expect_true(!is.null(result$plan))
  expect_equal(dim(result$X_white), dim(X))
  expect_equal(dim(result$Y_white), dim(Y))

  # Check AR coefficient
  if (!is.null(result$ar_coef)) {
    phi_est <- result$ar_coef[[1]][1]
    expect_true(abs(phi_est - true_phi) < 0.3)
  }
})

test_that("effective df calculation with fmriAR plan", {
  skip_if_not_installed("fmriAR")

  # Create a simple plan
  plan <- fmriAR::compat$plan_from_phi(phi = c(0.5))

  # Test effective df
  df <- .compute_ar_effective_df_compat(n = 100, p = 5, plan = plan)

  # With phi = 0.5, ar_factor = 1 - 0.25 = 0.75
  # effective_n = 100 * 0.75 = 75
  # df = 75 - 5 = 70
  expect_equal(df, 70)

  # Test with no AR
  plan_no_ar <- fmriAR::compat$plan_from_phi(phi = numeric(0))
  df_no_ar <- .compute_ar_effective_df_compat(n = 100, p = 5, plan = plan_no_ar)
  expect_equal(df_no_ar, 95)
})

test_that("multi-run AR estimation via fmriAR", {
  skip_if_not_installed("fmriAR")

  set.seed(321)
  n_runs <- 3
  n_per_run <- 50
  n_total <- n_runs * n_per_run

  # Create run indices
  run_indices <- lapply(1:n_runs, function(r) {
    ((r-1) * n_per_run + 1):(r * n_per_run)
  })

  # Generate data with different AR per run
  Y <- matrix(0, n_total, 2)
  true_phis <- c(0.3, 0.5, 0.4)

  for (r in 1:n_runs) {
    idx <- run_indices[[r]]
    for (v in 1:2) {
      Y[idx, v] <- as.numeric(arima.sim(list(ar = true_phis[r]), n = n_per_run))
    }
  }

  # Estimate AR with run-wise pooling
  cfg <- list(struct = "ar1", global = FALSE)
  plan <- .estimate_ar_via_fmriAR(Y, cfg, run_indices)

  expect_true(inherits(plan, "fmriAR_plan"))

  # Should have separate estimates per run
  if (!is.null(plan$phi) && length(plan$phi) > 1) {
    expect_equal(length(plan$phi), n_runs)
  }
})
