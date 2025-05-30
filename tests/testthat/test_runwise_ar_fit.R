context("runwise AR prewhitening")

options(mc.cores=1)

etab <- data.frame(onset=c(1,10), repnum=factor(c("A","B")), run=c(1,1))
Y <- matrix(rnorm(20*3), 20, 3)
dset <- matrix_dataset(Y, TR=1, run_length=20, event_table=etab)

test_that("fmri_lm runwise AR1 fits without error", {
  # Use ar_options parameter
  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            use_fast_path = TRUE,
            ar_options = list(
              struct = "ar1",
              iter_gls = 2,
              exact_first = TRUE
            )),
    NA
  )
})

test_that("fmri_lm runwise AR1 global fits without error", {
  # Use ar_options parameter with global pooling
  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            use_fast_path = TRUE,
            ar_options = list(
              struct = "ar1",
              global = TRUE
            )),
    NA
  )
})
