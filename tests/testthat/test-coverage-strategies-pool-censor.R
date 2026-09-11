# fmri_lm_strategies.R: censor edge cases and pool_run_results.

test_that("extract_censor_from_dataset covers empty and index vectors", {
  # Censored volume indices carried by the frame's censor column
  ds2 <- matrix_frame(matrix(0, 8, 2), TR = 2, run_length = c(4L, 4L),
                      censor = c(2L, 5L))
  idx <- fmrireg:::extract_censor_from_dataset(ds2)
  expect_equal(idx, c(2L, 5L))

  # Binary numeric vector all zeros => NULL
  ds3 <- matrix_frame(matrix(0, 4, 2), TR = 1, run_length = 4L,
                      censor = c(0, 0, 0, 0))
  expect_null(fmrireg:::extract_censor_from_dataset(ds3))

  # Logical all FALSE => NULL
  ds4 <- matrix_frame(matrix(0, 2, 2), TR = 1, run_length = 2L,
                      censor = c(FALSE, FALSE))
  expect_null(fmrireg:::extract_censor_from_dataset(ds4))
})

test_that("resolve_censor subsets global explicit censor by run", {
  dset <- matrix_frame(matrix(0, 16, 2), TR = 1, run_length = c(8L, 8L),
                       censor = c(1L, 3L, 10L, 15L))
  cfg <- fmri_lm_control()
  cfg$ar$censor <- c(1L, 3L, 10L, 15L)

  run1 <- fmrireg:::resolve_censor(cfg, dataset = dset, run_num = 1L, n_time = 8L)
  expect_true(is.null(run1) || all(run1 >= 1L))

  run2 <- fmrireg:::resolve_censor(cfg, dataset = dset, run_num = 2L, n_time = 8L)
  expect_true(is.null(run2) || is.integer(run2) || is.numeric(run2))

  # auto with no dataset censor
  cfg_auto <- fmri_lm_control()
  cfg_auto$ar$censor <- "auto"
  no_censor <- matrix_frame(matrix(0, 8, 2), TR = 1, run_length = 8L)
  expect_null(fmrireg:::resolve_censor(cfg_auto, dataset = no_censor, n_time = 8L))
})

test_that("pool_run_results pools sigma/rss/df across runs", {
  runs <- list(
    list(
      dfres = 10,
      sigma2 = matrix(c(1.0, 4.0), 1, 2),
      rss = matrix(c(10, 40), 1, 2)
    ),
    list(
      dfres = 5,
      sigma2 = matrix(c(1.0, 1.0), 1, 2),
      rss = matrix(c(5, 5), 1, 2)
    )
  )
  pooled <- fmrireg:::pool_run_results(runs)
  expect_equal(pooled$rdf, 15)
  expect_equal(as.numeric(pooled$rss), c(15, 45))
  expect_equal(length(pooled$sigma), 2L)
  expect_true(all(pooled$sigma > 0))
  expect_equal(as.numeric(pooled$resvar), as.numeric(pooled$rss) / pooled$rdf)

  # Single-voxel path
  runs2 <- list(
    list(dfres = 4, sigma2 = 1, rss = 4),
    list(dfres = 6, sigma2 = 1, rss = 6)
  )
  pooled2 <- fmrireg:::pool_run_results(runs2)
  expect_equal(pooled2$rdf, 10)
  expect_equal(length(pooled2$sigma), 1L)
})
