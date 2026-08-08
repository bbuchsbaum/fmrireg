test_that("voxelwise AR preserves degenerate voxels as explicit NA results", {
  skip_if_not_installed("fmridataset")
  set.seed(9251)

  n <- 40L
  events <- data.frame(
    onset = c(2, 10, 18, 26),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  signal <- as.numeric(arima.sim(model = list(ar = 0.25), n = n))
  Y <- cbind(constant = rep(1, n), signal = signal)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = events)

  fit <- fmri_lm(
    onset ~ hrf(condition),
    block = ~run,
    dataset = dset,
    ar_options = list(struct = "ar1", voxelwise = TRUE),
    robust_options = list(type = "huber", max_iter = 3L),
    strategy = "runwise",
    progress = FALSE
  )

  expect_identical(fit$result$voxel_status, c("constant", "ok"))
  expect_true(all(is.na(coef(fit)[, 1L])))
  expect_true(any(is.finite(coef(fit)[, 2L])))
  expect_true(is.list(fit$result$robust_weights))
  expect_null(fit$result$robust_weights[[1L]])
  expect_length(fit$result$robust_weights[[2L]], n)

  beta_data <- fit$result$betas$data[[1L]]
  expect_true(all(is.na(beta_data$se[[1L]][1L, ])))
  expect_true(all(is.na(beta_data$prob[[1L]][1L, ])))
})

test_that("future voxel collection preserves status and robust weights", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("fmridataset")
  set.seed(9252)

  n <- 36L
  events <- data.frame(
    onset = c(2, 9, 16, 23),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- cbind(zero = numeric(n), signal = rnorm(n))
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = events)
  args <- list(
    formula = onset ~ hrf(condition),
    block = ~run,
    dataset = dset,
    ar_options = list(struct = "ar1", voxelwise = TRUE),
    robust_options = list(type = "huber", max_iter = 2L),
    strategy = "runwise",
    progress = FALSE
  )

  sequential <- do.call(fmri_lm, c(args, list(parallel_voxels = FALSE)))
  future_path <- suppressWarnings(
    do.call(fmri_lm, c(args, list(parallel_voxels = TRUE)))
  )

  expect_identical(future_path$result$voxel_status, c("all_zero", "ok"))
  expect_equal(coef(future_path), coef(sequential), tolerance = 1e-10)
  expect_equal(future_path$result$robust_weights,
               sequential$result$robust_weights,
               tolerance = 1e-10)
})

test_that("voxelwise beta statistics use retained robust weights for df", {
  B <- matrix(0.5, nrow = 1L, ncol = 1L)
  weights <- c(rep(1, 10), rep(0.1, 10))
  out <- fmrireg:::beta_stats_matrix_voxelwise(
    Betas = B,
    XtXinv_list = list(matrix(1, 1, 1)),
    sigma = 1,
    dfres = 19,
    varnames = "x",
    robust_weights_list = list(weights),
    ar_order = 0L
  )

  got <- out$data[[1L]]$prob[[1L]][1L, 1L]
  df_effective <- fmrireg:::calculate_effective_df(
    n = 20, p = 1, robust_weights = weights, method = "simple"
  )
  expected <- 2 * stats::pt(-0.5, df = df_effective)
  nominal <- 2 * stats::pt(-0.5, df = 19)

  expect_equal(unname(got), unname(expected), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(got, nominal, tolerance = 1e-12)))
})
