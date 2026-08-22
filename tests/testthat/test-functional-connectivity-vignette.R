.simulate_ar1_errors <- function(n_time, n_voxels, rho, innovation_sd = 1) {
  innovations <- matrix(
    stats::rnorm(n_time * n_voxels, sd = innovation_sd),
    nrow = n_time, ncol = n_voxels
  )
  errors <- matrix(0, nrow = n_time, ncol = n_voxels)
  errors[1, ] <- innovations[1, ] / sqrt(1 - rho^2)
  for (time_index in 2:n_time) {
    errors[time_index, ] <-
      rho * errors[time_index - 1L, ] + innovations[time_index, ]
  }
  errors
}

.run_connectivity_fixture <- function(simulation_seed, n_time = 96L,
                                      n_voxels = 64L, n_active = 12L,
                                      effect = 0.8, rho = 0.3) {
  set.seed(simulation_seed)
  seed_ts <- as.numeric(base::scale(stats::arima.sim(
    model = list(ar = 0.5), n = n_time
  )))
  null_probe <- as.numeric(base::scale(stats::arima.sim(
    model = list(ar = 0.4), n = n_time
  )))
  network <- sample(seq_len(n_voxels), n_active)
  data <- .simulate_ar1_errors(n_time, n_voxels, rho)
  data[, network] <- data[, network] + effect * seed_ts

  dataset <- fmridataset::matrix_dataset(
    data, TR = 2, run_length = n_time,
    event_table = data.frame(onset = 0, run = 1L)
  )
  sframe <- sampling_frame(n_time, TR = 2)
  cov_data <- data.frame(seed = seed_ts, null_probe = null_probe)
  event <- event_model(
    onset ~
      covariate(seed, data = cov_data, prefix = "seed") +
      covariate(null_probe, data = cov_data, prefix = "null"),
    data = data.frame(onset = samples(sframe)[1], run = 1L),
    block = ~ run, sampling_frame = sframe
  )
  baseline <- baseline_model(basis = "constant", sframe = sframe)
  model <- fmri_model(event, baseline, dataset)
  fit <- fmri_lm(
    model, dataset = dataset,
    control = fmri_lm_control(
      noise = noise_spec("ar1", pooling = "global")
    )
  )

  pvals <- as.matrix(p_values(fit, type = "estimates"))
  seed_col <- grep("^seed", colnames(pvals), value = TRUE)
  null_col <- grep("^null", colnames(pvals), value = TRUE)
  stopifnot(length(seed_col) == 1L, length(null_col) == 1L)

  seed_discovery <- stats::p.adjust(
    pvals[, seed_col], method = "BH"
  ) < 0.05
  null_discovery <- stats::p.adjust(
    pvals[, null_col], method = "BH"
  ) < 0.05
  n_seed_discoveries <- sum(seed_discovery)

  c(
    sensitivity = mean(seed_discovery[network]),
    background_fpr = mean(seed_discovery[-network]),
    fdp = sum(seed_discovery[-network]) / max(1L, n_seed_discoveries),
    complete_null_rejection = as.numeric(any(null_discovery)),
    fitted_rho = as.numeric(ar_parameters(fit, scope = "raw")[[1]])
  )
}

.one_sided_mean_bounds <- function(x, level = 0.95) {
  critical <- stats::qt(level, df = length(x) - 1L)
  standard_error <- stats::sd(x) / sqrt(length(x))
  c(
    lower = mean(x) - critical * standard_error,
    upper = mean(x) + critical * standard_error
  )
}

.one_sided_binomial_upper <- function(successes, trials, level = 0.95) {
  if (successes == trials) return(1)
  stats::qbeta(level, successes + 1, trials - successes)
}

test_that("seed connectivity fixture matches its declared AR model", {
  skip_on_cran()
  skip_if_not_installed("fmrihrf")

  expect_no_warning(
    result <- .run_connectivity_fixture(42L)
  )
  expect_gte(unname(result[["sensitivity"]]), 0.8)
  expect_lte(unname(result[["background_fpr"]]), 0.05)
  expect_lt(abs(unname(result[["fitted_rho"]]) - 0.3), 0.12)
})

test_that("AR plus BH connectivity is calibrated across independent fixtures", {
  skip_on_cran()
  skip_if_not_installed("fmrihrf")

  n_replicates <- 200L
  captured_warnings <- character()
  results <- withCallingHandlers(
    vapply(
      1000L + seq_len(n_replicates),
      .run_connectivity_fixture,
      numeric(5)
    ),
    warning = function(condition) {
      captured_warnings <<- c(captured_warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(captured_warnings, 0L)

  fdp_upper <- .one_sided_mean_bounds(results["fdp", ])[["upper"]]
  sensitivity_lower <-
    .one_sided_mean_bounds(results["sensitivity", ])[["lower"]]
  background_fpr_upper <-
    .one_sided_mean_bounds(results["background_fpr", ])[["upper"]]
  complete_null_upper <- .one_sided_binomial_upper(
    sum(results["complete_null_rejection", ]), n_replicates
  )

  expect_lte(fdp_upper, 0.08)
  expect_gte(sensitivity_lower, 0.8)
  expect_lte(background_fpr_upper, 0.025)
  expect_lte(complete_null_upper, 0.12)
  expect_lt(abs(mean(results["fitted_rho", ]) - 0.3), 0.08)
})
