simulate_glm_ar_dataset <- function(ar_coeff = numeric(), n_runs = 2, n_time = 120,
                                    n_vox = 3, amplitudes = c(1, 0.8, 0.6), seed = 123) {
  set.seed(seed)
  TR <- 1
  run_length <- rep(n_time, n_runs)
  sframe <- fmrihrf::sampling_frame(run_length, TR)

  onset_cond1 <- c(6, 18, 30, 42)
  onset_cond2 <- c(12, 24, 36, 48)
  event_table <- do.call(rbind, lapply(seq_len(n_runs), function(r) {
    df <- data.frame(
      run = r,
      onset = c(onset_cond1, onset_cond2),
      condition = factor(c(rep("condition1", length(onset_cond1)),
                           rep("condition2", length(onset_cond2))))
    )
    df[order(df$onset), ]
  }))

  time_points <- fmrihrf::samples(sframe, global = TRUE)
  global_onsets <- fmrihrf::global_onsets(sframe, event_table$onset, event_table$run)
  reg1 <- fmrihrf::regressor(global_onsets[event_table$condition == "condition1"], fmrihrf::HRF_SPMG1, amplitude = 1)
  reg2 <- fmrihrf::regressor(global_onsets[event_table$condition == "condition2"], fmrihrf::HRF_SPMG1, amplitude = 1.5)
  base_signal <- fmrihrf::evaluate(reg1, time_points) + fmrihrf::evaluate(reg2, time_points)

  block_ids <- fmrihrf::blockids(sframe)
  datamat <- matrix(0, nrow = length(time_points), ncol = n_vox)
  for (v in seq_len(n_vox)) {
    scaled_signal <- base_signal * amplitudes[pmin(v, length(amplitudes))]
    noise_vec <- numeric(length(time_points))
    for (r in seq_len(n_runs)) {
      idx <- which(block_ids == r)
      noise_run <- if (length(ar_coeff) == 0) {
        rnorm(length(idx), sd = 0.3)
      } else {
        as.numeric(arima.sim(model = list(ar = ar_coeff), n = length(idx), sd = 0.3))
      }
      noise_vec[idx] <- noise_run
    }
    datamat[, v] <- scaled_signal + noise_vec
  }

  fmridataset::matrix_dataset(datamat, TR = TR, run_length = run_length, event_table = event_table)
}

test_that("fmri_lm chunkwise recovers AR coefficients across structures", {
  skip_if_not_installed("fmriAR")

  cases <- list(
    list(struct = "ar1", phi = c(0.5), args = list(), tol = 0.15),
    list(struct = "ar2", phi = c(0.6, -0.2), args = list(), tol = 0.25),
    list(struct = "arp", phi = c(0.5, -0.15, 0.08), args = list(p = 3), tol = 0.3)
  )

  for (cfg in cases) {
    dset <- simulate_glm_ar_dataset(ar_coeff = cfg$phi, n_runs = 3)
    ar_opts <- c(list(struct = cfg$struct, iter_gls = 2), cfg$args)
    fit <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset,
                   strategy = "chunkwise", nchunks = 1, use_fast_path = FALSE,
                   ar_options = ar_opts)
    phi_hat <- ar_parameters(fit)
    expect_equal(length(phi_hat), length(cfg$phi))
    expect_equal(as.numeric(phi_hat), cfg$phi, tolerance = cfg$tol)
  }
})

test_that("fmri_lm runwise recovers AR1 coefficients", {
  skip_if_not_installed("fmriAR")
  dset <- simulate_glm_ar_dataset(ar_coeff = 0.4, n_runs = 2)
  fit <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset,
                 strategy = "runwise", use_fast_path = FALSE,
                 ar_options = list(struct = "ar1", iter_gls = 2))
  phi_hat <- ar_parameters(fit)
  expect_equal(length(phi_hat), 1)
  expect_equal(as.numeric(phi_hat), 0.4, tolerance = 0.2)
})

test_that("latent sketch engine integrates fmriAR whitening", {
  skip_if_not_installed("fmriAR")
  dset <- simulate_glm_ar_dataset(ar_coeff = 0.3, n_runs = 2, n_time = 80)
  fit <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset,
                 engine = "latent_sketch",
                 lowrank = list(time_sketch = list(method = "gaussian", m = 64L)),
                 ar_options = list(struct = "ar1", iter_gls = 2))

  tidy_stats <- tidy(fit, type = "estimates")
  expect_true(all(is.finite(tidy_stats$estimate)))
  expect_true(all(is.finite(tidy_stats$std_error)))

  phi_avg <- ar_parameters(fit)
  expect_true(is.null(phi_avg) || is.double(phi_avg))
  phi_raw <- ar_parameters(fit, scope = "raw")
  expect_true(is.list(phi_raw))
})
