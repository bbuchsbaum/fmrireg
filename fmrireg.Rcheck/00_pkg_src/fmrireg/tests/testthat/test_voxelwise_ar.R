context("voxelwise AR slow path")

options(mc.cores=1)

etab <- data.frame(onset = c(1,10), repnum = factor(c("A","B")), run = c(1,1))
Y <- matrix(rnorm(20*3), 20, 3)
dset <- matrix_dataset(Y, TR=1, run_length=20, event_table=etab)

con <- contrast_set(pair_contrast(~ repnum == "A", ~ repnum == "B", name="AvB"))

test_that("fmri_lm voxelwise AR1 fits and computes contrast", {
  mod <- fmri_lm(onset ~ hrf(repnum, contrasts = con), block = ~ run,
                 dataset = dset, use_fast_path = FALSE,
                 ar_voxelwise = TRUE,
                 cor_struct = "ar1")
  expect_equal(dim(coef(mod)), c(2, 3))
  ctab <- stats(mod, "contrasts")
  expect_equal(ncol(ctab), 1)
})

