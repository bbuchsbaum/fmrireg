context("Comprehensive voxelwise AR contrasts")

set.seed(123)

# Helper data setup
etab <- data.frame(onset = c(1, 10), cond = factor(c("A", "B")), run = c(1, 1))
con <- contrast_set(pair_contrast(~ cond == "A", ~ cond == "B", name = "AvB"))


test_that("single voxel voxelwise matches global AR", {
  y <- as.numeric(arima.sim(model = list(ar = 0.4), n = 20))
  dset <- matrix_frame(matrix(y, ncol = 1), TR = 1, run_length = 20,
                         event_table = etab)

  mod_global <- fmri_lm(onset ~ hrf(cond, contrasts = con), block = ~ run,
                        dataset = dset, use_fast_path = FALSE,
                        cor_struct = "ar1", ar_voxelwise = FALSE)
  mod_vox <- fmri_lm(onset ~ hrf(cond, contrasts = con), block = ~ run,
                     dataset = dset, use_fast_path = FALSE,
                     cor_struct = "ar1", ar_voxelwise = TRUE)
  expect_equal(coef(mod_vox), coef(mod_global), tolerance = 1e-1)
  expect_equal(stats(mod_vox, "contrasts"), stats(mod_global, "contrasts"),
               tolerance = 5e-1)
})
