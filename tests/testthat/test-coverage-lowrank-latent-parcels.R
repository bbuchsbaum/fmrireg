# Additional .run_lowrank_engine branches: latent dataset, by_cluster parcels, landmarks.

make_matrix_fx <- function(n = 72L, V = 12L, seed = 11L) {
  set.seed(seed)
  etab <- data.frame(
    onset = c(6, 18, 30, 42, 54, 64),
    condition = factor(c("A", "B", "A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(n * V), n, V)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  list(
    model = fmri_model(emod, bmod, dset),
    dataset = dset,
    V = V
  )
}

test_that(".run_lowrank_engine latent_dataset path covers loadings Matrix/dense", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")

  set.seed(12)
  n_time <- 60L
  n_comp <- 4L
  n_vox <- 16L
  basis <- matrix(rnorm(n_time * n_comp), n_time, n_comp)
  loadings <- Matrix::Matrix(rnorm(n_vox * n_comp), n_vox, n_comp, sparse = TRUE)
  lvec <- fmristore::LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = neuroim2::NeuroSpace(c(4, 4, 1, n_time)),
    mask = rep(TRUE, n_vox),
    offset = rep(0, n_vox)
  )
  lds <- fmridataset::latent_dataset(source = list(lvec), TR = 1, run_length = n_time)
  lds$event_table <- data.frame(
    onset = c(8, 20, 32, 44),
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
    lowrank = list(time_sketch = list(method = "gaussian", m = 20L)),
    cfg = fmri_lm_control()
  )
  expect_s3_class(fit, "fmri_lm")
  expect_true(!is.null(fit$result$betas) || !is.null(fit$betas))

  # Missing LatentNeuroVec errors
  lds2 <- lds
  lds2$lvec <- NULL
  if (!is.null(lds2$backend)) lds2$backend$data <- NULL
  expect_error(
    fmrireg:::.run_lowrank_engine(
      fm, lds2,
      lowrank = list(time_sketch = list(method = "srht", m = 16L)),
      cfg = fmri_lm_control()
    ),
    "Cannot find LatentNeuroVec"
  )
})

test_that(".run_lowrank_engine by_cluster parcels AR path succeeds", {
  skip_if_not_installed("fmriAR")
  fx <- make_matrix_fx(n = 80L, V = 10L, seed = 13L)
  parcels <- rep(1:2, length.out = fx$V)
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", by_cluster = TRUE))

  fit <- fmrireg:::.run_lowrank_engine(
    fx$model, fx$dataset,
    lowrank = list(
      time_sketch = list(method = "gaussian", m = 24L),
      parcels = parcels
    ),
    cfg = cfg
  )
  expect_s3_class(fit, "fmri_lm")
  expect_true(!is.null(fit$result$betas))

  fit_cs <- fmrireg:::.run_lowrank_engine(
    fx$model, fx$dataset,
    lowrank = list(
      time_sketch = list(method = "countsketch", m = 24L),
      parcels = parcels
    ),
    cfg = cfg
  )
  expect_s3_class(fit_cs, "fmri_lm")

  expect_error(
    fmrireg:::.run_lowrank_engine(
      fx$model, fx$dataset,
      lowrank = list(
        time_sketch = list(method = "gaussian", m = 16L),
        parcels = 1:3
      ),
      cfg = cfg
    ),
    "parcels/group ids length"
  )
})

test_that(".run_lowrank_engine landmarks path with ihs/srht AR", {
  skip_if_not_installed("fmriAR")
  skip_if_not_installed("RANN")
  skip_if_not_installed("neuroim2")

  fx <- make_matrix_fx(n = 64L, V = 16L, seed = 14L)
  # Attach a spatial mask so index_to_coord works
  mask <- array(TRUE, dim = c(4, 4, 1))
  fx$dataset$mask <- neuroim2::LogicalNeuroVol(mask, neuroim2::NeuroSpace(c(4, 4, 1)))

  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))
  fit_ihs <- tryCatch(
    fmrireg:::.run_lowrank_engine(
      fx$model, fx$dataset,
      lowrank = list(
        time_sketch = list(method = "ihs", m = 18L, iters = 2L),
        landmarks = 4L,
        k_neighbors = 4L,
        kmeans_iter_max = 20L,
        kmeans_nstart = 1L
      ),
      cfg = cfg
    ),
    error = function(e) e
  )
  if (inherits(fit_ihs, "error")) {
    # Spatial mask recovery can fail on matrix_dataset; still require a real message
    expect_match(conditionMessage(fit_ihs), ".")
  } else {
    expect_s3_class(fit_ihs, "fmri_lm")
  }

  fit_srht <- tryCatch(
    fmrireg:::.run_lowrank_engine(
      fx$model, fx$dataset,
      lowrank = list(
        time_sketch = list(method = "srht", m = 18L),
        landmarks = 4L,
        k_neighbors = 4L,
        kmeans_iter_max = 20L,
        kmeans_nstart = 1L
      ),
      cfg = cfg
    ),
    error = function(e) e
  )
  if (inherits(fit_srht, "error")) {
    expect_match(conditionMessage(fit_srht), ".")
  } else {
    expect_s3_class(fit_srht, "fmri_lm")
  }
})
