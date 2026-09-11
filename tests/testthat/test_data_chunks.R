# Chunk iteration over frames: fmrireg's runwise and chunkwise engines consume
# data_chunk objects (data, voxel_ind, row_ind, chunk_num) produced lazily
# from an fmri_frame.

test_that("runwise chunks cover each run in acquisition order", {
  n_time <- 100
  n_vox <- 10
  n_runs <- 2
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  run_length <- rep(n_time / n_runs, n_runs)
  dset <- matrix_frame(Y, TR = 1, run_length = run_length)

  chunks <- fmrireg:::.dset_run_chunks(dset)
  expect_length(chunks, n_runs)

  chunk1 <- chunks[[1]]
  expect_s3_class(chunk1, "data_chunk")
  expect_true(all(c("data", "voxel_ind", "row_ind", "chunk_num") %in% names(chunk1)))
  expect_equal(nrow(chunk1$data), n_time / n_runs)
  expect_equal(ncol(chunk1$data), n_vox)
  expect_equal(chunk1$chunk_num, 1)
  expect_equal(chunk1$row_ind, 1:(n_time / n_runs))
  expect_equal(chunk1$voxel_ind, seq_len(n_vox))
  expect_equal(unname(chunk1$data), Y[1:50, ])
  expect_equal(unname(chunks[[2]]$data), Y[51:100, ])
  expect_equal(chunks[[2]]$row_ind, 51:100)
})

test_that("a single chunk holds the whole matrix", {
  n_time <- 50
  n_vox <- 5
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  dset <- matrix_frame(Y, TR = 1, run_length = n_time)

  chunks <- fmrireg:::.dset_chunks(dset, nchunks = 1)
  expect_length(chunks, 1)
  chunk <- chunks[[1]]
  expect_s3_class(chunk, "data_chunk")
  expect_equal(dim(chunk$data), dim(Y))
  expect_equal(unname(chunk$data), Y)
  expect_equal(chunk$chunk_num, 1)
  expect_equal(chunk$voxel_ind, 1:n_vox)
  expect_equal(chunk$row_ind, 1:n_time)
})

test_that("feature-wise chunks partition the voxels", {
  n_time <- 50
  n_vox <- 20
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  dset <- matrix_frame(Y, TR = 1, run_length = n_time)

  chunks <- fmrireg:::.dset_chunks(dset, nchunks = 4)
  expect_length(chunks, 4)
  all_vox_ind <- unlist(lapply(chunks, function(ch) ch$voxel_ind))
  expect_equal(sort(all_vox_ind), 1:n_vox)
  for (i in 1:4) {
    expect_equal(nrow(chunks[[i]]$data), n_time)
    expect_true(ncol(chunks[[i]]$data) > 0)
    expect_equal(chunks[[i]]$chunk_num, i)
    expect_equal(unname(chunks[[i]]$data), Y[, chunks[[i]]$voxel_ind, drop = FALSE])
  }
})

test_that("requesting more chunks than voxels warns and caps the count", {
  Y <- matrix(rnorm(30), 10, 3)
  dset <- matrix_frame(Y, TR = 1, run_length = 10)
  expect_warning(chunks <- fmrireg:::.dset_chunks(dset, nchunks = 10),
                 "greater than number of voxels")
  expect_length(chunks, 3)
})

test_that("chunks read lazily from a volumetric frame and honour a feature view", {
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, c(5, 5, 4)), neuroim2::NeuroSpace(c(5, 5, 4)))
  ntp <- 30
  scans <- replicate(3, {
    neuroim2::SparseNeuroVec(matrix(rnorm(ntp * 100), ntp, 100),
                             space = neuroim2::NeuroSpace(c(5, 5, 4, ntp)), mask = mask)
  }, simplify = FALSE)
  dset <- neurovec_frame(scans, mask, TR = 2)

  rchunks <- fmrireg:::.dset_run_chunks(dset)
  expect_length(rchunks, 3)
  expect_equal(vapply(rchunks, function(ch) nrow(ch$data), integer(1)), rep(ntp, 3))
  expect_equal(range(unlist(lapply(rchunks, `[[`, "row_ind"))), c(1, 3 * ntp))
  expect_equal(unname(rchunks[[2]]$data), neuroim2::series(scans[[2]], 1:100))

  # the same iteration on a feature view only touches the selected voxels
  view <- dset[, 1:10]
  vchunks <- fmrireg:::.dset_run_chunks(view)
  expect_equal(ncol(vchunks[[1]]$data), 10)
  expect_equal(unname(vchunks[[3]]$data), neuroim2::series(scans[[3]], 1:10))

  fchunks <- fmrireg:::.dset_chunks(dset, nchunks = 7)
  expect_length(fchunks, 7)
  expect_equal(sort(unlist(lapply(fchunks, `[[`, "voxel_ind"))), 1:100)
})

test_that("data_chunk objects have the expected structure", {
  test_mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  chunk <- fmrireg:::.data_chunk(test_mat, voxel_ind = 1:10, row_ind = 1:10, chunk_num = 1)
  expect_s3_class(chunk, "data_chunk")
  expect_equal(chunk$data, test_mat)
  expect_equal(chunk$voxel_ind, 1:10)
  expect_equal(chunk$row_ind, 1:10)
  expect_equal(chunk$chunk_num, 1)
})
