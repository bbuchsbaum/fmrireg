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

test_that("fixed-order design correction uses a separate covariance-tail budget", {
  skip_if_not_installed("fmriAR", minimum_version = "0.3.3")

  set.seed(124)
  n <- 80L
  X <- cbind(1, rnorm(n), as.numeric(scale(seq_len(n))))
  residuals <- qr.resid(qr(X), matrix(rnorm(n * 5L), nrow = n))

  for (order in 1:2) {
    cfg <- list(struct = if (order == 1L) "ar1" else "ar2")
    actual <- .estimate_ar_via_fmriAR(residuals, cfg, design = X)
    budget <- .ar_correction_lag_budget(X, target_order = order)
    expect_gt(budget, order)
    oracle <- fmriAR::fit_noise(
      resid = residuals,
      method = "ar",
      p = order,
      p_max = order,
      pooling = "global",
      exact_first = "ar1",
      design = X,
      correction_max_lag = budget
    )
    expect_equal(actual$phi, oracle$phi, tolerance = 1e-12)
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
    theta = numeric(0),
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

test_that("iterative and integrated solvers honor the shared AR estimand", {
  skip_if_not_installed("fmriAR", minimum_version = "0.3.3")

  set.seed(790)
  n <- 300L
  X <- cbind(1, as.numeric(scale(seq_len(n))))
  common <- as.numeric(arima.sim(list(ar = 0.75), n = n))
  opposing <- 3 * as.numeric(arima.sim(list(ar = 0.05), n = n))
  Y <- cbind(common + opposing, common - opposing)

  pooled_cfg <- list(
    struct = "ar1", iter_gls = 1L, exact_first = FALSE,
    shared_estimator = "pooled_acvf"
  )
  mean_cfg <- pooled_cfg
  mean_cfg$shared_estimator <- "mean_series"

  pooled <- .iterative_ar_gls_via_fmriAR(X, Y, pooled_cfg, max_iter = 1L)
  mean_fit <- .iterative_ar_gls_via_fmriAR(X, Y, mean_cfg, max_iter = 1L)

  residuals <- Y - X %*% base::qr.solve(X, Y)
  expect_identical(.shared_ar_correction_design(residuals, pooled_cfg, X), X)
  expect_identical(.shared_ar_correction_design(residuals, mean_cfg, X), X)
  mean_oracle <- .estimate_ar_via_fmriAR(
    matrix(rowMeans(residuals), ncol = 1L), mean_cfg, design = X
  )
  expect_equal(mean_fit$plan$phi, mean_oracle$phi, tolerance = 1e-12)
  expect_gt(
    abs(as.numeric(mean_fit$plan$phi[[1]]) -
          as.numeric(pooled$plan$phi[[1]])),
    0.2
  )

  control <- fmri_lm_control(
    noise = noise_spec("ar1", shared_estimator = "mean_series")
  )
  integrated <- solve_integrated_glm(X, Y, control)
  expect_equal(integrated$ar_coef, mean_fit$ar_coef, tolerance = 1e-12)
})

test_that("global shared AR preserves run boundaries and design correction", {
  skip_if_not_installed("fmriAR", minimum_version = "0.3.3")

  set.seed(791)
  run_length <- 80L
  run_indices <- list(seq_len(run_length), run_length + seq_len(run_length))
  X <- cbind(
    intercept = 1,
    trend = rep(as.numeric(scale(seq_len(run_length))), 2L)
  )
  residuals <- qr.resid(qr(X), matrix(rnorm(160L * 4L), nrow = 160L))
  cfg <- list(
    struct = "ar1", global = TRUE, exact_first = FALSE,
    shared_estimator = "pooled_acvf"
  )

  expected <- .estimate_ar_via_fmriAR(
    residuals, cfg, run_indices = run_indices, design = X
  )
  actual <- .estimate_shared_ar_parameters(
    residuals, 1L, cfg, run_indices = run_indices, design = X
  )

  expect_equal(actual, as.numeric(expected$phi[[1]]), tolerance = 1e-12)
})

test_that("effective df calculation with fmriAR plan", {
  skip_if_not_installed("fmriAR")

  # Create a simple plan
  plan <- fmriAR::compat$plan_from_phi(phi = c(0.5), theta = numeric(0))

  # Test effective df
  df <- .compute_ar_effective_df_compat(n = 100, p = 5, plan = plan)

  # Match adapter's ACF-based effective-n calculation.
  rho <- 0.5^(seq_len(99))
  k <- seq_along(rho)
  denom <- 1 + 2 * sum((1 - k / 100) * rho)
  expected_df <- max(min(100 / denom, 100) - 5, 1)
  expect_equal(df, expected_df)

  # Test with no AR (construct minimal plan to avoid compat warnings).
  plan_no_ar <- structure(list(phi = list(numeric(0))), class = "fmriAR_plan")
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
