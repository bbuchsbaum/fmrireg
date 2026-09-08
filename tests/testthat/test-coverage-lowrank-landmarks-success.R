# Successful .run_lowrank_engine landmarks / parcels / latent dense paths.

make_mem_lowrank_fx <- function(dims = c(4L, 4L, 1L), n_time = 60L, seed = 41L) {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")
  skip_if_not_installed("RANN")
  set.seed(seed)
  V <- prod(dims)
  arr <- array(rnorm(V * n_time), c(dims, n_time))
  scan <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, n_time)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dim = dims))
  etab <- data.frame(
    onset = c(5, 15, 25, 35, 45, 55),
    condition = factor(rep(c("A", "B"), 3)),
    run = 1L
  )
  dset <- fmridataset::fmri_mem_dataset(
    scans = list(scan), mask = mask, TR = 1, event_table = etab
  )
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "constant", sframe = dset$sampling_frame)
  list(
    model = fmri_model(emod, bmod, dset),
    dataset = dset,
    V = V
  )
}

test_that(".run_lowrank_engine landmarks succeed for ihs/srht/gaussian/countsketch", {
  skip_if_not_installed("fmriAR")
  fx <- make_mem_lowrank_fx()
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))
  low_base <- list(
    landmarks = 4L,
    k_neighbors = 4L,
    kmeans_iter_max = 15L,
    kmeans_nstart = 1L
  )

  for (method in c("ihs", "srht", "gaussian", "countsketch")) {
    sk <- list(method = method, m = 18L)
    if (identical(method, "ihs")) sk$iters <- 2L
    fit <- fmrireg:::.run_lowrank_engine(
      fx$model, fx$dataset,
      lowrank = c(low_base, list(time_sketch = sk)),
      cfg = cfg
    )
    expect_s3_class(fit, "fmri_lm")
    expect_true(!is.null(fit$result$betas))
    expect_equal(attr(fit, "strategy"), "sketch")
    expect_equal(length(fit$sigma2), fx$V)
  }
})

test_that(".run_lowrank_engine IID landmarks (no AR) succeed", {
  fx <- make_mem_lowrank_fx(n_time = 56L, seed = 42L)
  fit <- fmrireg:::.run_lowrank_engine(
    fx$model, fx$dataset,
    lowrank = list(
      time_sketch = list(method = "gaussian", m = 16L),
      landmarks = 3L,
      k_neighbors = 3L,
      kmeans_iter_max = 10L,
      kmeans_nstart = 1L
    ),
    cfg = fmri_lm_control()
  )
  expect_s3_class(fit, "fmri_lm")
  expect_true(all(is.finite(fit$sigma2)))
})

test_that(".run_lowrank_engine parcels by_cluster with srht/ihs succeed", {
  skip_if_not_installed("fmriAR")
  set.seed(43)
  n <- 72L
  V <- 12L
  etab <- data.frame(
    onset = c(6, 18, 30, 42, 54, 64),
    condition = factor(rep(c("A", "B"), 3)),
    run = 1L
  )
  Y <- matrix(rnorm(n * V), n, V)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  fm <- fmri_model(emod, bmod, dset)
  parcels <- rep(1:3, length.out = V)
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", by_cluster = TRUE))

  for (method in c("srht", "ihs", "gaussian")) {
    sk <- list(method = method, m = 20L)
    if (identical(method, "ihs")) sk$iters <- 2L
    fit <- fmrireg:::.run_lowrank_engine(
      fm, dset,
      lowrank = list(time_sketch = sk, parcels = parcels),
      cfg = cfg
    )
    expect_s3_class(fit, "fmri_lm")
    expect_true(!is.null(fit$ar_coef) || !is.null(fit$result$ar_coef))
  }
})

test_that(".run_lowrank_engine latent dense loadings + cfg=NULL ar_options path", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")

  set.seed(44)
  n_time <- 50L
  n_comp <- 3L
  n_vox <- 12L
  basis <- matrix(rnorm(n_time * n_comp), n_time, n_comp)
  # Dense (non-Matrix) loadings hit the t(lds) branch
  loadings <- matrix(rnorm(n_vox * n_comp), n_vox, n_comp)
  lvec <- suppressWarnings(fmristore::LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = neuroim2::NeuroSpace(c(3, 4, 1, n_time)),
    mask = rep(TRUE, n_vox),
    offset = rep(0, n_vox)
  ))
  lds <- fmridataset::latent_dataset(source = list(lvec), TR = 1, run_length = n_time)
  lds$event_table <- data.frame(
    onset = c(6, 16, 26, 36),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  emod <- event_model(
    onset ~ hrf(condition), data = lds$event_table, block = ~ run,
    sampling_frame = lds$sampling_frame
  )
  bmod <- baseline_model(basis = "constant", sframe = lds$sampling_frame)
  fm <- fmri_model(emod, bmod, lds)

  fit <- fmrireg:::.run_lowrank_engine(
    fm, lds,
    lowrank = list(time_sketch = list(method = "srht", m = 16L)),
    cfg = NULL,
    ar_options = list(struct = "iid")
  )
  expect_s3_class(fit, "fmri_lm")
  expect_equal(length(fit$sigma2), n_vox)
})

test_that(".lowrank_blend_ar and .lowrank_project_voxels helpers", {
  # alpha * local + (1 - alpha) * global
  expect_equal(fmrireg:::.lowrank_blend_ar(0.4, 0.2, 0), 0.2)
  expect_equal(fmrireg:::.lowrank_blend_ar(0.4, 0.2, 1), 0.4)
  blended <- fmrireg:::.lowrank_blend_ar(c(0.5, 0.1), c(0.1, 0.3), 0.5)
  expect_equal(length(blended), 2L)
  expect_equal(
    fmrireg:::.lowrank_blend_ar(list(0.4), list(0.2), 0.5),
    list(0.3)
  )
  expect_error(
    fmrireg:::.lowrank_blend_ar(list(0.1, 0.2), list(0.3), 0.5),
    "different lengths"
  )

  M <- matrix(1:6, 2, 3) # p x r
  expect_equal(fmrireg:::.lowrank_project_voxels(M, NULL, TRUE), M)
  A <- matrix(c(1, 0, 0.5, 0, 1, 0.5), 3, 2) # r x V
  proj <- fmrireg:::.lowrank_project_voxels(M, A, FALSE)
  expect_equal(dim(proj), c(2L, 2L))
})
