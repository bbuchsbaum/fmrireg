context("Robust argument passing")

# simple synthetic dataset
etab <- data.frame(onset=c(1,10), repnum=factor(c("A","B")), run=c(1,1))
Y <- matrix(rnorm(20*5), 20, 5)
dset <- matrix_dataset(Y, TR=1, run_length=20, event_table=etab)

# accept valid robust arguments
test_that("fmri_lm accepts robust arguments", {
  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            robust = TRUE, robust_psi = "huber", robust_max_iter = 2,
            robust_scale_scope = "run"),
    NA
  )

  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            robust = TRUE, robust_psi = "bisquare", robust_max_iter = 1,
            robust_scale_scope = "global"),
    NA
  )
})

# validation of arguments
test_that("robust argument validation works", {
  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            robust = TRUE, robust_psi = "bogus"),
    "robust_psi"
  )

  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            robust = TRUE, robust_max_iter = 0),
    "robust_max_iter"
  )
})
