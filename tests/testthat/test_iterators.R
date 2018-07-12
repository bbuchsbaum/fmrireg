

library(foreach)

gen_dataset <- function(nruns, ntp, nvox=1000) {
  mask <- neuroim2::LogicalNeuroVol(array(1, c(10,10,10)), space=neuroim2::NeuroSpace(c(10,10,10)))
  vec <- neuroim2::SparseNeuroVec(matrix(rnorm(ntp*nvox), ntp, nvox), space=neuroim2::NeuroSpace(c(10,10,10,ntp)), mask=mask)
  scans <- replicate(nruns, vec, simplify=FALSE)
  
  dset <- fmri_mem_dataset(scans, mask, TR=2)
}

test_that("can construct a runwise iterator", {
  dset <- gen_dataset(5, 100, nvox=1000)
  rchunks <- data_chunks(dset, runwise=TRUE)
  res <- foreach (chunk = rchunks) %do% {
    list(ncol=ncol(chunk$data), nrow=nrow(chunk$data))
  }
  
  expect_equal(length(res), 5)
  expect_true(all(sapply(res, "[[", "ncol") == 1000))
  expect_true(all(sapply(res, "[[", "nrow") == 100))
})