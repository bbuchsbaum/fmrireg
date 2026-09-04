# Thirteenth wave: internal LM voxel status, storage helpers, AR config attach.

test_that(".fmri_lm_voxel_status covers finite/NA/Inf series", {
  ok <- fmrireg:::.fmri_lm_voxel_status(rnorm(20))
  bad_na <- fmrireg:::.fmri_lm_voxel_status(c(1, NA, 3))
  bad_inf <- fmrireg:::.fmri_lm_voxel_status(c(1, Inf, 3))
  expect_false(identical(ok, bad_na))
  expect_false(identical(ok, bad_inf))
})

test_that("fmri_lm with robust + AR1 covers combined strategy routing", {
  skip_if_not_installed("fmriAR")
  skip_if_not_installed("robustbase")
  set.seed(321)
  n <- 50L
  etab <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(n * 3), n, 3)
  Y[1, 1] <- Y[1, 1] + 20
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  fit <- fmri_lm(
    onset ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    ar_options = list(struct = "ar1"),
    robust = TRUE,
    robust_options = list(type = "huber", max_iter = 2L),
    strategy = "runwise"
  )
  expect_s3_class(fit, "fmri_lm")
  expect_true(!is.null(coef(fit)))
})

test_that("fmri_lm chunkwise with AR2 and progress=FALSE", {
  skip_if_not_installed("fmriAR")
  set.seed(322)
  n_per <- 30L
  etab <- data.frame(
    onset = c(5, 15, 25, 5, 15, 25),
    condition = factor(rep(c("A", "B"), 3)),
    run = rep(1:2, each = 3)
  )
  Y <- matrix(rnorm((2 * n_per) * 3), 2 * n_per, 3)
  dset <- matrix_dataset(Y, TR = 1, run_length = c(n_per, n_per), event_table = etab)
  fit <- fmri_lm(
    onset ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    strategy = "chunkwise",
    nchunks = 2L,
    ar_options = list(struct = "ar2"),
    progress = FALSE
  )
  expect_s3_class(fit, "fmri_lm")
  cf <- coef(fit, include_baseline = TRUE)
  expect_true(is.matrix(cf) || inherits(cf, "tbl_df") || is.data.frame(cf))
})

test_that(".rrr_normalize_args and extract response matrix edges", {
  args <- fmrireg:::.rrr_normalize_args(list(rank = 3L))
  expect_true(is.list(args))
  expect_true(!is.null(args$rank) || !is.null(args$rank_mode))

  set.seed(323)
  Y <- matrix(rnorm(40 * 4), 40, 4)
  dset <- matrix_dataset(
    Y, TR = 1, run_length = 40L,
    event_table = data.frame(onset = 5, condition = "A", run = 1L)
  )
  mat <- fmrireg:::.rrr_extract_response_matrix(dset)
  expect_equal(dim(mat), dim(Y))
})
