# Frame constructors: fmrireg builds fmridataset::fmri_frame objects from
# matrices, NeuroVec objects, NIfTI files, and latent decompositions.

test_that("matrix_frame builds a frame with the run structure and events", {
  Y <- matrix(seq_len(40), 20, 2)
  events <- data.frame(onset = c(2, 12), condition = factor(c("A", "B")), run = c(1, 2))
  frame <- matrix_frame(Y, TR = 2, run_length = c(10, 10), event_table = events)

  expect_s3_class(frame, "fmri_frame")
  expect_equal(dim(frame), c(20L, 2L))
  schema <- fmridataset::temporal_schema(frame)
  expect_equal(unname(schema$run_lengths), c(10L, 10L))
  expect_equal(unname(schema$TR), c(2, 2))
  expect_equal(schema$block_ids, rep(1:2, each = 10))
  expect_null(schema$censor)

  sf <- fmridataset::as_sampling_frame(frame)
  expect_equal(fmrihrf::blocklens(sf), c(10, 10))
  expect_equal(fmrihrf::blockids(sf), rep(1:2, each = 10))

  expect_equal(unname(fmridataset::collect_assay(frame)), Y)
  ev <- frame_events(frame)
  expect_equal(ev$onset, events$onset)
  expect_equal(as.character(ev$condition), c("A", "B"))
  expect_true("event_id" %in% names(ev))

  expect_equal(fmridataset::observation_ids(frame)[c(1, 11)],
               c("run-1-vol-00001", "run-2-vol-00001"))
  expect_equal(fmridataset::feature_ids(frame), c("feature-1", "feature-2"))
})

test_that("matrix_frame validates run lengths and stores censor as a logical column", {
  Y <- matrix(0, 10, 1)
  expect_error(matrix_frame(Y, TR = 1, run_length = c(4, 4)), "sum\\(run_length\\)")

  f01 <- matrix_frame(Y, TR = 1, run_length = 10, censor = c(0, 0, 1, 0, 0, 0, 1, 0, 0, 0))
  expect_equal(which(fmridataset::temporal_schema(f01)$censor), c(3L, 7L))
  fidx <- matrix_frame(Y, TR = 1, run_length = 10, censor = c(3L, 7L))
  expect_equal(which(fmridataset::temporal_schema(fidx)$censor), c(3L, 7L))
  flog <- matrix_frame(Y, TR = 1, run_length = 10, censor = seq_len(10) %in% c(3, 7))
  expect_equal(fmridataset::temporal_schema(flog)$censor,
               fmridataset::temporal_schema(f01)$censor)
})

test_that("matrix_frame uses unique column names as feature IDs and keeps a per-run TR", {
  Y <- matrix(rnorm(12), 6, 2, dimnames = list(NULL, c("roi_a", "roi_b")))
  frame <- matrix_frame(Y, TR = c(1, 2), run_length = c(3, 3))
  expect_equal(fmridataset::feature_ids(frame), c("roi_a", "roi_b"))
  expect_equal(unname(fmridataset::temporal_schema(frame)$TR), c(1, 2))
})

test_that("neurovec_frame extracts masked voxel series onto a volume_space", {
  facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header = TRUE)
  facedes$repnum <- factor(facedes$rep_num)
  n_runs <- length(unique(facedes$run))

  set.seed(1)
  scans <- lapply(seq_len(n_runs), function(i) {
    arr <- array(rnorm(4 * 4 * 3 * 20), c(4, 4, 3, 20))
    neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(4, 4, 3, 20)))
  })
  mask <- neuroim2::LogicalNeuroVol(array(rnorm(4 * 4 * 3), c(4, 4, 3)) > 0,
                                    neuroim2::NeuroSpace(dim = c(4, 4, 3)))

  frame <- neurovec_frame(scans = scans, mask = mask, TR = 1.5,
                          event_table = tibble::as_tibble(facedes))

  expect_s3_class(frame, "fmri_frame")
  expect_s3_class(fmridataset::space(frame), "volume_space")
  expect_equal(nrow(frame), 20L * n_runs)
  expect_equal(ncol(frame), sum(mask))
  expect_equal(unname(fmridataset::temporal_schema(frame)$run_lengths), rep(20L, n_runs))

  # values agree with neuroim2::series() on the masked voxels
  voxel_ind <- which(as.vector(mask) != 0)
  expected <- do.call(rbind, lapply(scans, function(s) neuroim2::series(s, voxel_ind)))
  expect_equal(unname(fmridataset::collect_assay(frame)), unname(expected))

  # the mask and NeuroSpace round-trip through the volume_space
  spatial <- fmrireg:::.fmri_dataset_mask_space(frame, "test")
  expect_equal(spatial$mask_array, array(as.vector(mask) != 0, dim(mask)))
  expect_equal(dim(spatial$space), dim(mask))
})

test_that("nifti_frame builds a lazy file-backed frame", {
  skip_if_not_installed("RNifti")
  dims <- c(4, 4, 3)
  nii_dir <- tempfile("nifti-frame-")
  dir.create(nii_dir)
  on.exit(unlink(nii_dir, recursive = TRUE))

  set.seed(2)
  scans <- lapply(1:2, function(run) {
    neuroim2::NeuroVec(array(rnorm(prod(dims) * 6), c(dims, 6)),
                       neuroim2::NeuroSpace(dim = c(dims, 6)))
  })
  files <- vapply(seq_along(scans), function(i) {
    f <- file.path(nii_dir, sprintf("run-%d.nii", i))
    neuroim2::write_vec(scans[[i]], f)
    f
  }, character(1))
  mask_array <- array(FALSE, dims)
  mask_array[1:10] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(dim = dims))
  mask_file <- file.path(nii_dir, "mask.nii")
  neuroim2::write_vol(mask, mask_file)

  frame <- nifti_frame(basename(files), basename(mask_file), TR = 2,
                       base_path = nii_dir,
                       event_table = data.frame(onset = c(1, 3), run = c(1, 2)))
  expect_s3_class(frame, "fmri_frame")
  expect_equal(dim(frame), c(12L, 10L))
  expect_equal(unname(fmridataset::temporal_schema(frame)$run_lengths), c(6L, 6L))

  expected <- do.call(rbind, lapply(scans, function(s) neuroim2::series(s, 1:10)))
  expect_equal(unname(fmridataset::collect_assay(frame)), unname(expected), tolerance = 1e-6)

  # a NeuroVol mask is accepted too and yields the same space
  frame2 <- nifti_frame(files, mask, TR = 2, run_length = c(6, 6))
  expect_true(fmridataset::same_space(fmridataset::space(frame), fmridataset::space(frame2))$same)

  expect_error(nifti_frame(files, mask_file, TR = 2, run_length = c(5, 5)), "sum\\(run_length\\)")
})

test_that("latent_frame stores scores as the assay and loadings on a basis_space", {
  skip_if_not_installed("fmristore")
  n_time <- 12; n_comp <- 3; n_vox <- 24
  basis <- matrix(rnorm(n_time * n_comp), n_time, n_comp)
  loadings <- matrix(rnorm(n_vox * n_comp), n_vox, n_comp)
  lvec <- fmristore::LatentNeuroVec(
    basis = basis, loadings = loadings,
    space = neuroim2::NeuroSpace(c(4, 3, 2, n_time)),
    mask = rep(TRUE, n_vox), offset = rep(0, n_vox)
  )
  frame <- latent_frame(lvec, TR = 2, run_length = c(6, 6))
  expect_s3_class(frame, "fmri_frame")
  expect_true(fmrireg:::.dset_is_latent(frame))
  expect_equal(unname(fmridataset::collect_assay(frame)), unname(basis))
  expect_equal(unname(as.matrix(fmrireg:::.dset_loadings(frame))), unname(loadings))
  expect_s3_class(fmridataset::parent_space(fmridataset::space(frame)), "volume_space")
  expect_equal(dim(fmrireg:::.fmri_dataset_mask_space(frame, "test")$mask_array), c(4L, 3L, 2L))
})
