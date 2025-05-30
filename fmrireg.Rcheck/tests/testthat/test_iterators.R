options(mc.cores=2)

library(foreach)

gen_dataset <- function(nruns, ntp, nvox=1000) {
  mask <- neuroim2::LogicalNeuroVol(array(1, c(10,10,10)), space=neuroim2::NeuroSpace(c(10,10,10)))
  vec <- neuroim2::SparseNeuroVec(matrix(rnorm(ntp*nvox), ntp, nvox), space=neuroim2::NeuroSpace(c(10,10,10,ntp)), mask=mask)
  scans <- replicate(nruns, vec, simplify=FALSE)
  
  dset <- fmri_mem_dataset(scans, mask, TR=2)
}

test_that("runwise iterator maintains data integrity", {
  dset <- gen_dataset(3, 100, nvox=1000)
  rchunks <- data_chunks(dset, runwise=TRUE)
  
  # Test sequential access
  all_data <- list()
  foreach::foreach(chunk=rchunks) %do% {
    all_data[[length(all_data) + 1]] <- chunk$data
  }
  
  # Verify dimensions
  expect_equal(length(all_data), 3)
  expect_equal(nrow(all_data[[1]]), 100)
  expect_equal(ncol(all_data[[1]]), 1000)
  
  # Test parallel access with foreach
  res1 <- foreach::foreach(chunk = rchunks) %do% {
    colMeans(chunk$data)
  }
  
  # Test parallel access again to ensure iterator resets correctly
  res2 <- foreach::foreach(chunk = rchunks) %do% {
    colMeans(chunk$data)
  }
  
  # Results should be identical across iterations
  expect_equal(res1, res2)
})

test_that("arbitrary chunking works correctly", {
  dset <- gen_dataset(2, 100, nvox=1000)
  
  # Test different chunk sizes
  chunk_sizes <- c(2, 5, 10, 20)
  
  for (nchunks in chunk_sizes) {
    rchunks <- data_chunks(dset, nchunks=nchunks)
    
    # Collect all chunk data
    res <- foreach::foreach(chunk = rchunks) %do% {
      list(
        data_dim = dim(chunk$data),
        voxel_indices = range(chunk$voxel_ind),
        row_indices = range(chunk$row_ind),
        chunk_num = chunk$chunk_num
      )
    }
    
    # Verify chunk count
    expect_equal(length(res), nchunks)
    
    # Verify chunk numbers are sequential
    expect_equal(sapply(res, function(x) x$chunk_num), 1:nchunks)
    
    # Verify all voxels are covered
    all_voxels <- sort(unique(unlist(lapply(res, function(x) x$voxel_indices[1]:x$voxel_indices[2]))))
    expect_equal(all_voxels, 1:1000)
  }
})

test_that("iterator handles edge cases correctly", {
  # Test with single run
  dset_single <- gen_dataset(1, 100, nvox=1000)
  chunks_single <- data_chunks(dset_single, runwise=TRUE)
  res_single <- foreach::foreach(chunk = chunks_single) %do% dim(chunk$data)
  expect_equal(length(res_single), 1)
  
  # Test with single voxel
  #dset_small <- gen_dataset(2, 100, nvox=1)
  #chunks_small <- data_chunks(dset_small, nchunks=1)
  #res_small <- foreach::foreach(chunk = chunks_small) %do% dim(chunk$data)
  #expect_equal(res_small[[1]][2], 1)
  
  # Test with more chunks than voxels
  #dset_over <- gen_dataset(2, 100, nvox=5)
  #chunks_over <- data_chunks(dset_over, nchunks=10)
  #res_over <- foreach::foreach(chunk = chunks_over) %do% dim(chunk$data)
  #expect_equal(length(res_over), 5)  # Should limit to number of voxels
})

test_that("matrix_dataset chunking works correctly", {
  # Create a matrix dataset
  data_mat <- matrix(rnorm(1000 * 100), 100, 1000)
  sframe <- sampling_frame(c(50, 50), TR=2)
  mset <- matrix_dataset(data_mat, TR=2, run_length=c(50,50))
  
  # Test runwise chunking
  rchunks <- data_chunks.matrix_dataset(mset, runwise=TRUE)
  res_run <- foreach(chunk = rchunks) %do% {
    list(dim=dim(chunk$data), 
         row_indices=range(chunk$row_ind),
         chunk_num=chunk$chunk_num)
  }
  
  expect_equal(length(res_run), 2)  # Two runs
  expect_equal(res_run[[1]]$dim[1], 50)  # First run length
  expect_equal(res_run[[2]]$dim[1], 50)  # Second run length
  
  # Test arbitrary chunking
  chunks_arb <- data_chunks.matrix_dataset(mset, nchunks=4)
  res_arb <- foreach(chunk = chunks_arb) %do% {
    list(dim=dim(chunk$data),
         voxel_indices=range(chunk$voxel_ind),
         chunk_num=chunk$chunk_num)
  }
  
  expect_equal(length(res_arb), 4)
  # Verify voxel coverage
  all_voxels <- sort(unique(unlist(lapply(res_arb, function(x) x$voxel_indices[1]:x$voxel_indices[2]))))
  expect_equal(length(all_voxels), 1000)
})

test_that("parallel processing with foreach works correctly", {
  dset <- gen_dataset(4, 100, nvox=1000)
  
  # Test parallel processing with different chunk sizes
  chunk_sizes <- c(2, 4, 8)
  
  for (nchunks in chunk_sizes) {
    chunks <- data_chunks(dset, nchunks=nchunks)
    
    # Run parallel computation
    res_par <- foreach::foreach(chunk = chunks) %dopar% {
      colMeans(chunk$data)
    }
    
    chunks <- data_chunks(dset, nchunks=nchunks)
    # Run sequential computation
    res_seq <- foreach::foreach(chunk = chunks) %do% {
      colMeans(chunk$data)
    }
    
    # Results should be the same
    expect_equal(res_par, res_seq)
    expect_equal(length(res_par), nchunks)
  }
})

test_that("iterator reset functionality works", {
  dset <- gen_dataset(3, 100, nvox=1000)
  chunks <- data_chunks(dset, nchunks=5)
  
  # First iteration
  res1 <- foreach::foreach(chunk = chunks) %do% {
    sum(chunk$data)
  }
  
  chunks =  data_chunks(dset, nchunks=5)
  # Second iteration
  res2 <- foreach::foreach(chunk = chunks) %do% {
    sum(chunk$data)
  }
  
  # Results should be identical
  expect_equal(res1, res2)
  expect_equal(length(res1), 5)
  
  chunks =  data_chunks(dset, nchunks=5)
  
  # Test mixed foreach operators
  res3 <- foreach(chunk = chunks) %do% {
    sum(chunk$data)
  }
  
  chunks =  data_chunks(dset, nchunks=5)
  
  res4 <- foreach(chunk = chunks) %dopar% {
    sum(chunk$data)
  }
  
  expect_equal(res3, res4)
})

test_that("chunk indices are correct and complete", {
  dset <- gen_dataset(2, 100, nvox=1000)
  
  # Test runwise chunks
  rchunks <- data_chunks(dset, runwise=TRUE)
  run_indices <- list()
  
  foreach(chunk=rchunks) %do% {
    run_indices[[length(run_indices) + 1]] <- chunk$row_ind
  }
  
  # Verify run indices
  expect_equal(length(unlist(run_indices)), 200)  # Total timepoints
  expect_equal(range(unlist(run_indices)), c(1, 200))
  expect_equal(length(unique(unlist(run_indices))), 200)  # No duplicates
  
  # Test arbitrary chunks
  chunks <- data_chunks(dset, nchunks=3)
  voxel_indices <- list()
  
  foreach(chunk = chunks) %do% {
    voxel_indices[[length(voxel_indices) + 1]] <- chunk$voxel_ind
  }
  
  # Verify voxel indices
  all_voxels <- sort(unique(unlist(voxel_indices)))
  expect_equal(range(all_voxels), c(1, 1000))
  expect_equal(length(all_voxels), 1000)  # All voxels covered
})

test_that("can construct and iterate over a runwise iterator", {
  dset <- gen_dataset(5, 100, nvox=1000)
  rchunks <- data_chunks(dset, runwise=TRUE)
  res <- foreach (chunk = rchunks) %do% {
    list(ncol=ncol(chunk$data), nrow=nrow(chunk$data))
  }
  
  expect_equal(length(res), 5)
  expect_true(all(sapply(res, "[[", "ncol") == 1000))
  expect_true(all(sapply(res, "[[", "nrow") == 100))
})

test_that("can construct and iterate over a runwise-chunked iterator", {
  dset <- gen_dataset(5, 100, nvox=1000)
  rchunks <- data_chunks(dset, nchunks=5, runwise=TRUE)
  res <- foreach (chunk = rchunks) %do% {
    list(ncol=ncol(chunk$data), nrow=nrow(chunk$data))
  }
  
  expect_equal(length(res), 5)
  expect_true(all(sapply(res, "[[", "ncol") == 1000))
  expect_true(all(sapply(res, "[[", "nrow") == 100))
})