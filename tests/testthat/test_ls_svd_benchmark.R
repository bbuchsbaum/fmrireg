context("ls_svd engine with benchmark dataset")

library(fmrireg)

# This test loads the BM_Canonical_HighSNR benchmark dataset and
# runs the LS+SVD+1ALS engine via fmrireg_cfals. Only a subset of voxels
# are used to keep the test lightweight.

test_that("ls_svd_1als_engine works on benchmark data", {
  skip_if_not_installed("fmrireg")

  bm <- load_benchmark_dataset("BM_Canonical_HighSNR")
  dset <- bm$core_data

  # Use only the first 5 voxels for speed
  Y <- dset$datamat[, 1:5]

  evtab <- dset$event_table
  evtab$block <- 1
  emod <- event_model(onset ~ hrf(condition), data = evtab,
                      block = ~ block, sampling_frame = dset$sampling_frame)

  fit <- fmrireg_cfals(Y, emod, HRF_SPMG1,
                       method = "ls_svd_1als",
                       lambda_init = 0.1,
                       lambda_b = 0.1,
                       lambda_h = 0.1)

  expect_equal(nrow(fit$h_coeffs), nbasis(HRF_SPMG1))
  expect_equal(ncol(fit$h_coeffs), ncol(Y))
  expect_equal(dim(fit$beta_amps), c(length(unique(evtab$condition)), ncol(Y)))
  expect_true(all(is.finite(fit$h_coeffs)))
})
