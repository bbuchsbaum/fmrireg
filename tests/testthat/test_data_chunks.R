# Test data chunking functionality

library(fmrireg)

test_that("matrix_dataset chunking works correctly", {
  # Create test data
  n_time <- 100
  n_vox <- 10
  n_runs <- 2
  
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  run_length <- rep(n_time/n_runs, n_runs)
  
  dset <- matrix_dataset(Y, TR = 1, run_length = run_length)
  
  # Test runwise chunking
  chunks <- data_chunks(dset, runwise = TRUE)
  expect_s3_class(chunks, "chunkiter")
  
  # Should have 2 chunks (one per run)
  expect_equal(chunks$nchunks, 2)
  
  # Collect all chunks
  chunk_list <- list()
  for (i in 1:chunks$nchunks) {
    chunk_list[[i]] <- chunks$nextElem()
  }
  
  expect_equal(length(chunk_list), n_runs)
  
  # Check first chunk structure
  chunk1 <- chunk_list[[1]]
  expect_s3_class(chunk1, "data_chunk")
  expect_true(all(c("data", "voxel_ind", "row_ind", "chunk_num") %in% names(chunk1)))
  
  # Check dimensions
  expect_equal(nrow(chunk1$data), n_time/n_runs)
  expect_equal(ncol(chunk1$data), n_vox)
  expect_equal(chunk1$chunk_num, 1)
  expect_equal(chunk1$row_ind, 1:(n_time/n_runs))
})

test_that("matrix_dataset single chunk works", {
  n_time <- 50
  n_vox <- 5
  
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  dset <- matrix_dataset(Y, TR = 1, run_length = n_time)
  
  chunks <- data_chunks(dset, nchunks = 1)
  chunk <- chunks$nextElem()
  
  expect_s3_class(chunk, "data_chunk")
  expect_equal(dim(chunk$data), dim(Y))
  expect_equal(chunk$chunk_num, 1)
  expect_equal(chunk$voxel_ind, 1:n_vox)
})

test_that("matrix_dataset voxel chunking works", {
  n_time <- 50
  n_vox = 20
  
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  dset <- matrix_dataset(Y, TR = 1, run_length = n_time)
  
  # Split into 4 chunks
  chunks <- data_chunks(dset, nchunks = 4)
  expect_equal(chunks$nchunks, 4)
  
  chunk_list <- list()
  for (i in 1:chunks$nchunks) {
    chunk_list[[i]] <- chunks$nextElem()
  }
  
  expect_equal(length(chunk_list), 4)
  
  # Check that all voxels are covered
  all_vox_ind <- unlist(lapply(chunk_list, function(ch) ch$voxel_ind))
  expect_equal(sort(all_vox_ind), 1:n_vox)
  
  # Check chunk dimensions
  for (i in 1:4) {
    expect_equal(nrow(chunk_list[[i]]$data), n_time)
    expect_true(ncol(chunk_list[[i]]$data) > 0)
    expect_equal(chunk_list[[i]]$chunk_num, i)
  }
})

test_that("data_chunk object has correct structure", {
  skip("data_chunk is an internal function in fmridataset, not exported")
})