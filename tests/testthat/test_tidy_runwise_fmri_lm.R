test_that("tidy fmri_lm normalizes runwise estimate and contrast layouts", {
  set.seed(20260809)
  samples_per_run <- 50L
  n_vox <- 3L
  events <- data.frame(
    onset = c(5, 18, 5, 18),
    condition = factor(rep(c("A", "B"), 2)),
    run = rep(1:2, each = 2)
  )
  Y <- matrix(rnorm(2L * samples_per_run * n_vox), ncol = n_vox)
  dset <- matrix_dataset(
    Y,
    TR = 1,
    run_length = rep(samples_per_run, 2L),
    event_table = events
  )
  con <- pair_contrast(
    ~ condition == "B", ~ condition == "A", name = "B_minus_A"
  )
  fit <- fmri_lm(
    onset ~ hrf(condition, contrasts = con),
    block = ~ run,
    dataset = dset,
    strategy = "runwise",
    nchunks = 1L
  )

  estimates <- tidy(fit, type = "estimates")
  contrasts <- tidy(fit, type = "contrasts")
  expect_s3_class(estimates, "tbl_df")
  expect_s3_class(contrasts, "tbl_df")
  expect_equal(nrow(dplyr::filter(estimates, term %in% c("A", "B"))), 2L * n_vox)
  expect_equal(nrow(contrasts), n_vox)

  beta <- t(coef(fit))
  colnames(beta) <- c("A", "B")
  se_beta <- as.matrix(standard_error(fit))
  stat_beta <- as.matrix(stats(fit))
  p_beta <- as.matrix(p_values(fit))
  colnames(se_beta) <- colnames(stat_beta) <- colnames(p_beta) <- c("A", "B")
  event_rows <- dplyr::filter(estimates, term %in% c("A", "B"))
  event_cell <- cbind(event_rows$voxel, match(event_rows$term, colnames(beta)))
  expect_equal(event_rows$estimate, unname(beta[event_cell]), tolerance = 1e-12)
  expect_equal(event_rows$std_error, unname(se_beta[event_cell]), tolerance = 1e-12)
  expect_equal(event_rows$statistic, unname(stat_beta[event_cell]), tolerance = 1e-12)
  expect_equal(event_rows$p_value, unname(p_beta[event_cell]), tolerance = 1e-12)

  con_coef <- as.matrix(coef(fit, type = "contrasts"))
  con_se <- as.matrix(standard_error(fit, type = "contrasts"))
  con_stat <- as.matrix(stats(fit, type = "contrasts"))
  con_p <- as.matrix(p_values(fit, type = "contrasts"))
  con_cell <- cbind(contrasts$voxel, match(contrasts$term, colnames(con_coef)))
  expect_equal(contrasts$estimate, unname(con_coef[con_cell]), tolerance = 1e-12)
  expect_equal(contrasts$std_error, unname(con_se[con_cell]), tolerance = 1e-12)
  expect_equal(contrasts$statistic, unname(con_stat[con_cell]), tolerance = 1e-12)
  expect_equal(contrasts$p_value, unname(con_p[con_cell]), tolerance = 1e-12)
})
