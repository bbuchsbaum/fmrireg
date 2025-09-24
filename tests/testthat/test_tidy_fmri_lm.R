test_that("tidy fmri_lm produces expected estimate statistics", {
  set.seed(123)
  TR <- 2
  run_length <- 60
  n_vox <- 3

  event_table <- data.frame(
    run = 1,
    onset = c(10, 20, 30, 40, 50, 55),
    condition = factor(rep(c("condition1", "condition2"), each = 3))
  )

  sframe <- sampling_frame(blocklens = run_length, TR = TR)
  time_points <- samples(sframe, global = TRUE)
  global_onsets <- fmrihrf::global_onsets(sframe, event_table$onset, event_table$run)

  reg1 <- regressor(global_onsets[event_table$condition == "condition1"], fmrihrf::HRF_SPMG1, amplitude = 1)
  reg2 <- regressor(global_onsets[event_table$condition == "condition2"], fmrihrf::HRF_SPMG1, amplitude = 1.5)
  signal <- evaluate(reg1, time_points) + evaluate(reg2, time_points)
  noise <- simulate_noise_vector(length(time_points), TR = TR, ar = 0.4, sd = 0.3, physio = FALSE)

  base_series <- signal + noise
  datamat <- sapply(seq_len(n_vox), function(v) base_series * (1 - (v - 1) * 0.1))

  dset <- fmridataset::matrix_dataset(datamat = datamat, TR = TR,
                                      run_length = run_length, event_table = event_table)

  fit <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset,
                 strategy = "chunkwise", nchunks = 1)

  tidy_est <- tidy(fit, type = "estimates")
  expect_s3_class(tidy_est, "tbl_df")
  expect_true(all(c("voxel", "term", "estimate", "std_error", "statistic", "p_value") %in% names(tidy_est)))
  cond_rows <- dplyr::filter(tidy_est, grepl("condition", term))
  expect_true(all(is.finite(cond_rows$estimate)))
  expect_true(all(is.finite(cond_rows$std_error)))
  expect_true(all(is.finite(cond_rows$statistic)))

  con_spec <- pair_contrast(~ condition == "condition2", ~ condition == "condition1", name = "cond2_minus_cond1")
  fit_con <- fmri_lm(onset ~ hrf(condition, contrasts = con_spec), block = ~ run,
                     dataset = dset, strategy = "chunkwise", nchunks = 1)

  tidy_con <- tidy(fit_con, type = "contrasts")
  expect_s3_class(tidy_con, "tbl_df")
  expect_true("cond2_minus_cond1" %in% tidy_con$term)
  expect_true(all(is.finite(tidy_con$estimate)))
})
