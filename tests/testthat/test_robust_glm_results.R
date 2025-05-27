context("Robust GLM integration")

simulate_spike_dataset <- function(n_time = 40, n_vox = 3, spike = FALSE, seed = 1) {
  set.seed(seed)
  onsets <- c(5, 15, 25, 35)
  event_tab <- data.frame(onset = onsets, cond = factor("A"), run = 1)
  base <- matrix_dataset(matrix(rnorm(n_time * n_vox, sd = 0.1), n_time, n_vox),
                         TR = 1, run_length = n_time, event_table = event_tab)
  model <- create_fmri_model(onset ~ hrf(cond), block = ~ run,
                             dataset = base, durations = 0)
  X <- design_matrix(model)
  ev_cols <- unlist(attr(model$event_model$design_matrix, "col_indices"))
  beta <- rep(0, ncol(X))
  beta[ev_cols] <- 1
  Y <- X %*% beta + matrix(rnorm(n_time * n_vox, sd = 0.1), n_time, n_vox)
  if (spike) {
    Y[10, ] <- Y[10, ] + 10
  }
  dset <- matrix_dataset(Y, TR = 1, run_length = n_time, event_table = event_tab)
  list(dset = dset, beta_true = 1, ev_cols = ev_cols)
}

# No outlier: robust vs OLS should match

test_that("robust=TRUE matches OLS on clean data", {
  sim <- simulate_spike_dataset(spike = FALSE)
  mod_ols <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = sim$dset,
                     durations = 0, use_fast_path = TRUE)
  mod_rb <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = sim$dset,
                    durations = 0, use_fast_path = TRUE,
                    robust = TRUE)
  b_ols <- as.numeric(coef(mod_ols)[, 1])
  b_rb  <- as.numeric(coef(mod_rb)[, 1])
  expect_equal(b_rb, b_ols, tolerance = 0.05)
  se_ols <- as.numeric(standard_error(mod_ols)[, 1])
  se_rb  <- as.numeric(standard_error(mod_rb)[, 1])
  expect_equal(se_rb, se_ols, tolerance = 0.05)
})

# Outlier: robust should be closer to clean fit

test_that("robust fitting downweights spikes", {
  sim_clean <- simulate_spike_dataset(spike = FALSE, seed = 2)
  sim_spike <- simulate_spike_dataset(spike = TRUE, seed = 2)
  mod_true <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = sim_clean$dset,
                      durations = 0, use_fast_path = TRUE)
  mod_ols <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = sim_spike$dset,
                     durations = 0, use_fast_path = TRUE)
  mod_rb <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = sim_spike$dset,
                    durations = 0, use_fast_path = TRUE,
                    robust = TRUE)
  b_true <- as.numeric(coef(mod_true)[, 1])
  err_ols <- abs(as.numeric(coef(mod_ols)[, 1]) - b_true)
  err_rb  <- abs(as.numeric(coef(mod_rb)[, 1]) - b_true)
  expect_lt(err_rb, err_ols)
  se_ols <- as.numeric(standard_error(mod_ols)[, 1])
  se_rb  <- as.numeric(standard_error(mod_rb)[, 1])
  expect_lt(se_rb, se_ols)
})

# Compare against MASS::rlm for single voxel

test_that("fast_rlm_run approximates MASS::rlm", {
  set.seed(10)
  n <- 30
  X <- cbind(1, rnorm(n))
  beta <- c(0.5, 2)
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.1))
  y[5] <- y[5] + 5
  proj <- fmrireg:::.fast_preproject(X)
  fit_fast <- fmrireg:::fast_rlm_run(X, matrix(y, ncol = 1), proj,
                                     psi = "huber", max_it = 2)
  fit_mass <- MASS::rlm(X, y, psi = MASS::psi.huber, k = 1.345)
  expect_equal(as.numeric(fit_fast$betas), as.numeric(fit_mass$coef), tolerance = 0.05)
})
