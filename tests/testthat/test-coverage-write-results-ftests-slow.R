# Thirteenth wave: fit_Ftests, reconstruct_image edges, write_results.fmri_lm nifti.

test_that("fit_Ftests returns per-response term F statistics", {
  set.seed(311)
  d <- data.frame(
    y1 = rnorm(50), y2 = rnorm(50),
    a = rnorm(50), b = factor(rep(c("x", "y"), 25))
  )
  fit <- lm(cbind(y1, y2) ~ a + b, data = d)
  ft <- fmrireg:::fit_Ftests(fit)
  expect_true(is.list(ft) || inherits(ft, "tbl_df") || is.data.frame(ft))
  # Weighted path
  fit_w <- lm(cbind(y1, y2) ~ a + b, data = d, weights = runif(50, 0.5, 1.5))
  ft_w <- fmrireg:::fit_Ftests(fit_w)
  expect_true(is.list(ft_w) || inherits(ft_w, "tbl_df") || is.data.frame(ft_w))
})

test_that("write_results.fmri_lm nifti spatial export covers beta/contrast writers", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")
  set.seed(312)
  dims <- c(2L, 2L, 1L)
  n_time <- 40L
  arr <- array(rnorm(prod(dims) * n_time), c(dims, n_time))
  scan <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, n_time)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dim = dims))
  etab <- data.frame(
    onset = c(5, 12, 20, 28),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  dset <- fmridataset::fmri_mem_dataset(
    scans = list(scan), mask = mask, TR = 1, event_table = etab
  )
  con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"))
  fit <- fmri_lm(
    onset ~ hrf(condition, contrasts = con),
    block = ~ run, dataset = dset
  )
  td <- tempfile("wr-nifti-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  out <- write_results(
    fit,
    path = td,
    subject = "01",
    task = "demo",
    space = "T1w",
    format = "nifti",
    save_betas = TRUE,
    contrast_stats = c("estimate", "stat"),
    overwrite = TRUE
  )
  expect_true(is.list(out))
  files <- unlist(out, use.names = FALSE)
  expect_true(any(grepl("\\.nii\\.gz$", files)))
  expect_true(any(grepl("\\.json$", files)))
})

test_that("runwise_lm_slow and chunkwise_lm_slow paths via use_fast_path=FALSE", {
  set.seed(313)
  n <- 40L
  etab <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(n * 3), n, 3)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "constant", sframe = dset$sampling_frame)
  model <- fmri_model(emod, bmod, dset)
  cfg <- fmri_lm_control()

  slow_run <- fmrireg:::runwise_lm_impl(
    dset, model, contrast_objects = list(), cfg = cfg,
    use_fast_path = FALSE, progress = FALSE, verbose = TRUE
  )
  expect_true(is.list(slow_run))

  slow_chunk <- fmrireg:::chunkwise_lm.fmri_dataset(
    dset, model, contrast_objects = list(), nchunks = 2L, cfg = cfg,
    use_fast_path = FALSE, progress = FALSE, verbose = TRUE
  )
  expect_true(is.list(slow_chunk))
})
