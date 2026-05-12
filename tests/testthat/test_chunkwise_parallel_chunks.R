test_that("parallel_chunks routes chunkwise fitting through future-compatible path", {
  skip_if_not_installed("future")

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)

  set.seed(185)
  event_data <- data.frame(
    onset = c(2, 10, 18, 26, 34, 42),
    amp = c(1, 0, 1, 0, 1, 0),
    run = rep(1:2, each = 3)
  )
  dset <- fmridataset::matrix_dataset(
    matrix(rnorm(120 * 4), nrow = 120, ncol = 4),
    TR = 1,
    run_length = c(60, 60),
    event_table = event_data
  )

  for (fast_path in c(FALSE, TRUE)) {
    fit_seq <- suppressWarnings(
      fmri_lm(
        onset ~ hrf(amp),
        block = ~run,
        dataset = dset,
        strategy = "chunkwise",
        nchunks = 2,
        use_fast_path = fast_path,
        parallel_chunks = FALSE,
        progress = FALSE
      )
    )
    fit_future <- suppressWarnings(
      fmri_lm(
        onset ~ hrf(amp),
        block = ~run,
        dataset = dset,
        strategy = "chunkwise",
        nchunks = 2,
        use_fast_path = fast_path,
        parallel_chunks = TRUE,
        progress = FALSE
      )
    )

    expect_equal(
      fit_future$result$betas$data[[1]]$estimate[[1]],
      fit_seq$result$betas$data[[1]]$estimate[[1]],
      tolerance = 1e-10
    )
    expect_equal(fit_future$result$cov.unscaled, fit_seq$result$cov.unscaled, tolerance = 1e-10)
    expect_equal(fit_future$result$event_indices, fit_seq$result$event_indices)
    expect_equal(fit_future$result$baseline_indices, fit_seq$result$baseline_indices)
  }
})
