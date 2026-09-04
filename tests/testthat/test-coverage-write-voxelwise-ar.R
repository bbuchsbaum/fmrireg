# Spatial write_results NIfTI contrast export + runwise voxelwise AR path.

create_spatial_fit <- function(dims = c(2, 2, 1), n_time = 40L) {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")
  set.seed(61)
  scans <- lapply(1:2, function(run) {
    arr <- array(rnorm(prod(dims) * n_time), c(dims, n_time))
    neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, n_time)))
  })
  mask <- neuroim2::LogicalNeuroVol(
    array(TRUE, dims),
    neuroim2::NeuroSpace(dim = dims)
  )
  event_table <- data.frame(
    onset = c(5, 15, 25, 35, 5, 15, 25, 35),
    condition = factor(rep(c("A", "B"), 4)),
    run = rep(1:2, each = 4)
  )
  dset <- fmridataset::fmri_mem_dataset(
    scans = scans, mask = mask, TR = 1, event_table = event_table
  )
  con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"))
  fit <- fmri_lm(
    onset ~ hrf(condition, contrasts = con),
    block = ~ run, dataset = dset, durations = 0
  )
  list(fit = fit, dataset = dset, dims = dims)
}

test_that("write_results.fmri_lm nifti by_contrast covers spatial export", {
  fx <- create_spatial_fit()
  tmp <- tempfile("wr_nifti_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  out <- write_results(
    fx$fit,
    path = tmp,
    subject = "01",
    task = "test",
    space = "T1w",
    strategy = "by_contrast",
    format = "nifti",
    contrast_stats = c("estimate", "se", "stat"),
    overwrite = TRUE
  )
  expect_true(is.list(out) || is.character(out) || !is.null(out))
  files <- list.files(tmp, recursive = TRUE)
  expect_true(length(files) >= 1L)
  expect_true(any(grepl("nii|json|A_vs_B|contrast", files, ignore.case = TRUE)))
})

test_that(".save_contrasts_by_contrast_nifti succeeds with spatial fmri_lm", {
  fx <- create_spatial_fit(dims = c(2, 2, 1), n_time = 36L)
  tmp <- tempfile("con_nifti_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  cons <- fx$fit$result$contrasts
  cons <- cons[cons$type == "contrast" | cons$type %in% c("contrast", "t"), , drop = FALSE]
  if (nrow(cons) == 0L) {
    # Build from any available contrast rows
    cons <- fx$fit$result$contrasts
  }
  expect_true(nrow(cons) >= 1L)

  entities <- list(sub = "01", task = "test", space = "T1w")
  created <- fmrireg:::.save_contrasts_by_contrast_nifti(
    fx$fit, cons, tmp, entities, desc = "statmap",
    contrast_stats = c("estimate", "stat", "se"),
    overwrite = TRUE, output_formats = "nifti"
  )
  expect_true(is.list(created))
  expect_true(length(created) >= 1L || length(list.files(tmp)) >= 1L)
})

test_that("runwise_lm_voxelwise AR(1) path on multi-run matrix data", {
  skip_if_not_installed("fmriAR")
  set.seed(62)
  n_per <- 36L
  V <- 4L
  etab <- data.frame(
    onset = c(5, 15, 25, 5, 15, 25),
    condition = factor(rep(c("A", "B"), 3)),
    run = rep(1:2, each = 3)
  )
  Y <- matrix(rnorm((2 * n_per) * V), 2 * n_per, V)
  dset <- matrix_dataset(Y, TR = 1, run_length = c(n_per, n_per), event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  model <- fmri_model(emod, bmod, dset)

  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", voxelwise = TRUE))
  # Drive through runwise_lm_impl which dispatches to voxelwise when configured
  res <- fmrireg:::runwise_lm_impl(
    dset, model, contrast_objects = list(), cfg = cfg,
    use_fast_path = TRUE, progress = FALSE, verbose = TRUE
  )
  expect_true(is.list(res))
  expect_true(!is.null(res$betas) || !is.null(res$event_indices) || !is.null(res$sigma))
})

test_that("process_run_ar_robust covers AR+robust run fitting", {
  skip_if_not_installed("fmriAR")
  skip_if_not_installed("robustbase")
  set.seed(63)
  n <- 50L
  V <- 3L
  etab <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(n * V), n, V)
  Y[1:2, 1] <- Y[1:2, 1] + 15
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "constant", sframe = dset$sampling_frame)
  model <- fmri_model(emod, bmod, dset)

  cfg <- fmri_lm_control(
    ar_options = list(struct = "ar1"),
    robust = robust_spec(type = "huber", max_iter = 2L)
  )
  chunk <- list(
    data = Y,
    chunk_num = 1L,
    row_ind = seq_len(n)
  )
  out <- fmrireg:::process_run_ar_robust(chunk, model, cfg)
  expect_true(is.list(out))
  expect_true(all(c("betas", "sigma2", "XtXinv", "robust_weights",
                    "phi_hat", "ar_order") %in% names(out)))
  expect_equal(out$ar_order, 1L)
})
