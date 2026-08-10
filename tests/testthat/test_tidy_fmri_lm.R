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

  clean_term <- function(x) {
    x <- gsub("^conditioncondition_condition\\.", "", x)
    x <- gsub("^condition_condition\\.", "", x)
    gsub("\\.", " ", x)
  }
  betas <- t(coef(fit))
  ses <- as.matrix(standard_error(fit))
  tstats <- as.matrix(stats(fit))
  pvals <- as.matrix(p_values(fit))
  colnames(betas) <- clean_term(colnames(betas))
  colnames(ses) <- clean_term(colnames(ses))
  colnames(tstats) <- clean_term(colnames(tstats))
  colnames(pvals) <- clean_term(colnames(pvals))
  cell <- cbind(cond_rows$voxel, match(cond_rows$term, colnames(betas)))

  expect_equal(cond_rows$estimate, unname(betas[cell]), tolerance = 1e-12)
  expect_equal(cond_rows$std_error, unname(ses[cell]), tolerance = 1e-12)
  expect_equal(cond_rows$statistic, unname(tstats[cell]), tolerance = 1e-12)
  expect_equal(cond_rows$p_value, unname(pvals[cell]), tolerance = 1e-12)

  inference_df <- as.numeric(fit$result$df$inference)
  if (length(inference_df) == n_vox) {
    expect_equal(cond_rows$df_inference, inference_df[cond_rows$voxel])
  }

  con_spec <- pair_contrast(~ condition == "condition2", ~ condition == "condition1", name = "cond2_minus_cond1")
  fit_con <- fmri_lm(onset ~ hrf(condition, contrasts = con_spec), block = ~ run,
                     dataset = dset, strategy = "chunkwise", nchunks = 1)

  tidy_con <- tidy(fit_con, type = "contrasts")
  expect_s3_class(tidy_con, "tbl_df")
  expect_true("cond2_minus_cond1" %in% tidy_con$term)
  expect_true(all(is.finite(tidy_con$estimate)))

  con_coef <- as.matrix(coef(fit_con, type = "contrasts"))
  con_se <- as.matrix(standard_error(fit_con, type = "contrasts"))
  con_stat <- as.matrix(stats(fit_con, type = "contrasts"))
  con_p <- as.matrix(p_values(fit_con, type = "contrasts"))
  con_cell <- cbind(tidy_con$voxel, match(tidy_con$term, colnames(con_coef)))
  expect_equal(tidy_con$estimate, unname(con_coef[con_cell]), tolerance = 1e-12)
  expect_equal(tidy_con$std_error, unname(con_se[con_cell]), tolerance = 1e-12)
  expect_equal(tidy_con$statistic, unname(con_stat[con_cell]), tolerance = 1e-12)
  expect_equal(tidy_con$p_value, unname(con_p[con_cell]), tolerance = 1e-12)
})
