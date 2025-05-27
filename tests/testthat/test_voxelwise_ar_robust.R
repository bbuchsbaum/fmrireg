context("voxelwise AR+robust")

options(mc.cores=1)

etab <- data.frame(onset = c(1,10), repnum = factor(c("A","B")), run = c(1,1))
Y <- matrix(rnorm(20*4), 20, 4)
dset <- matrix_dataset(Y, TR=1, run_length=20, event_table=etab)

con <- contrast_set(pair_contrast(~ repnum == "A", ~ repnum == "B", name="AvB"))

rob_opts <- list(type = "huber", max_iter = 3)

# Simple check that model runs and produces expected dimensions

test_that("fmri_lm voxelwise AR1 robust fits", {
  mod <- fmri_lm(onset ~ hrf(repnum, contrasts = con), block = ~ run,
                 dataset = dset, use_fast_path = FALSE,
                 ar_voxelwise = TRUE,
                 cor_struct = "ar1",
                 robust = "huber",
                 robust_options = rob_opts)
  expect_equal(dim(coef(mod)), c(2, 4))
  ctab <- stats(mod, "contrasts")
  expect_equal(ncol(ctab), 1)
})
