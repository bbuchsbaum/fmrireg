context("ls_svd engine with benchmark dataset")

library(fmrireg)

# This test loads the BM_Canonical_HighSNR benchmark dataset and
# tests event model creation with benchmark data. The actual LS+SVD+1ALS 
# engine (fmrireg_cfals) is not yet implemented, so we test the infrastructure.

test_that("benchmark data can be used with event_model infrastructure", {
  skip_if_not_installed("fmrireg")

  bm <- load_benchmark_dataset("BM_Canonical_HighSNR")
  dset <- bm$core_data

  # Use only the first 5 voxels for speed
  Y <- dset$datamat[, 1:5]

  # Create the missing sampling_frame from benchmark data
  # The benchmark has TR and run_length information
  TR <- dset$TR
  run_length <- dset$run_length
  sampling_frame <- sampling_frame(blocklens = run_length, TR = TR)

  evtab <- dset$event_table
  evtab$block <- 1
  emod <- event_model(onset ~ hrf(condition), data = evtab,
                      block = ~ block, sampling_frame = sampling_frame)

  # Test that event model was created successfully
  expect_s3_class(emod, "event_model")
  expect_equal(nrow(design_matrix(emod)), run_length)
  expect_true(ncol(design_matrix(emod)) > 0)
  
  # Test that conditions are correctly extracted
  conds <- conditions(emod)
  expect_true(length(conds) > 0)
  expect_true(all(c("condition.Cond1", "condition.Cond2", "condition.Cond3") %in% conds))
})


