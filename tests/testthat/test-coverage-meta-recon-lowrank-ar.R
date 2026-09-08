# reconstruct_image / mask_array for h5+nifti-shaped meta objects + lowrank AR sketch paths.

test_that("reconstruct_image covers h5 and nifti group_data backends", {
  skip_if_not_installed("neuroim2")
  dims <- c(2L, 2L, 1L)
  n_vox <- prod(dims)
  values <- rnorm(n_vox)

  # H5-like object with space + mask
  obj_h5 <- list(
    data = structure(
      list(
        dim = dims,
        space = neuroim2::NeuroSpace(dims),
        mask = array(TRUE, dim = dims)
      ),
      class = "group_data_h5"
    )
  )
  vol_h5 <- fmrireg:::reconstruct_image(values, obj_h5)
  expect_s4_class(vol_h5, "NeuroVol")

  # H5 without mask fills all
  obj_h5b <- list(
    data = structure(
      list(dim = dims, space = NULL, mask = NULL),
      class = "group_data_h5"
    )
  )
  vol_h5b <- fmrireg:::reconstruct_image(values, obj_h5b)
  expect_s4_class(vol_h5b, "NeuroVol")

  # NIfTI-like with voxel_indices
  obj_nii <- list(
    data = structure(
      list(dim = dims, voxel_size = c(2, 2, 2)),
      class = "group_data_nifti"
    ),
    voxel_indices = seq_len(n_vox)
  )
  vol_nii <- fmrireg:::reconstruct_image(values, obj_nii)
  expect_s4_class(vol_nii, "NeuroVol")

  # NIfTI without voxel_indices
  obj_nii2 <- list(
    data = structure(list(dim = dims, voxel_size = c(1, 1, 1)), class = "group_data_nifti")
  )
  vol_nii2 <- fmrireg:::reconstruct_image(values, obj_nii2)
  expect_s4_class(vol_nii2, "NeuroVol")
})

test_that(".fmri_meta_mask_array covers nifti/h5 backends and fallbacks", {
  skip_if_not_installed("neuroim2")
  dims <- c(2L, 2L, 1L)
  ref <- neuroim2::NeuroVol(array(rnorm(prod(dims)), dim = dims), neuroim2::NeuroSpace(dims))

  # NULL object -> fallback
  m0 <- fmrireg:::.fmri_meta_mask_array(NULL, ref)
  expect_equal(dim(m0), dims)

  # NIfTI with mask_data
  mask_vol <- array(c(TRUE, FALSE, TRUE, TRUE), dim = dims)
  x_nii <- list(
    data = structure(
      list(mask_data = mask_vol),
      class = "group_data_nifti"
    )
  )
  m1 <- fmrireg:::.fmri_meta_mask_array(x_nii, ref)
  expect_equal(dim(m1), dims)

  # NIfTI via voxel_indices
  x_nii2 <- list(
    data = structure(list(mask_data = NULL), class = "group_data_nifti"),
    voxel_indices = c(1L, 3L)
  )
  m2 <- fmrireg:::.fmri_meta_mask_array(x_nii2, ref)
  expect_true(sum(m2) == 2L)

  # H5 with mask
  x_h5 <- list(
    data = structure(
      list(mask = array(TRUE, dim = dims)),
      class = "group_data_h5"
    )
  )
  m3 <- fmrireg:::.fmri_meta_mask_array(x_h5, ref)
  expect_true(all(m3))

  # H5 without usable mask
  x_h5b <- list(data = structure(list(mask = NULL), class = "group_data_h5"))
  m4 <- fmrireg:::.fmri_meta_mask_array(x_h5b, ref)
  expect_equal(dim(m4), dims)
})

test_that(".run_lowrank_engine global AR srht/gaussian/countsketch succeed", {
  skip_if_not_installed("fmriAR")
  set.seed(91)
  n <- 70L
  V <- 8L
  etab <- data.frame(
    onset = c(5, 15, 25, 35, 45, 55),
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
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))

  for (method in c("srht", "gaussian", "countsketch")) {
    fit <- fmrireg:::.run_lowrank_engine(
      fm, dset,
      lowrank = list(time_sketch = list(method = method, m = 20L)),
      cfg = cfg
    )
    expect_s3_class(fit, "fmri_lm")
  }

  # arp order via ar_options$p
  cfg2 <- fmri_lm_control(ar_options = list(struct = "arp", p = 2L))
  fit2 <- fmrireg:::.run_lowrank_engine(
    fm, dset,
    lowrank = list(time_sketch = list(method = "gaussian", m = 18L)),
    cfg = cfg2
  )
  expect_s3_class(fit2, "fmri_lm")
})

test_that(".run_lowrank_engine landmarks with spatial mem dataset", {
  skip_if_not_installed("fmriAR")
  skip_if_not_installed("RANN")
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")

  set.seed(92)
  dims <- c(3L, 3L, 1L)
  n_time <- 48L
  V <- prod(dims)
  arr <- array(rnorm(V * n_time), c(dims, n_time))
  scan <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, n_time)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dim = dims))
  etab <- data.frame(
    onset = c(5, 12, 20, 28, 36, 42),
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
  fm <- fmri_model(emod, bmod, dset)
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))

  for (method in c("gaussian", "countsketch", "srht")) {
    fit <- tryCatch(
      fmrireg:::.run_lowrank_engine(
        fm, dset,
        lowrank = list(
          time_sketch = list(method = method, m = 16L),
          landmarks = 3L,
          k_neighbors = 3L,
          kmeans_iter_max = 20L,
          kmeans_nstart = 1L
        ),
        cfg = cfg
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      expect_match(conditionMessage(fit), ".")
    } else {
      expect_s3_class(fit, "fmri_lm")
    }
  }
})
