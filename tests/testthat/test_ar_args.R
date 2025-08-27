context("AR argument passing")

options(mc.cores=1)

library(testthat)

# simple synthetic dataset
etab <- data.frame(onset=c(1,10), repnum=factor(c("A","B")), run=c(1,1))
Y <- matrix(rnorm(20*5), 20, 5)
dset <- fmridataset::matrix_dataset(Y, TR=1, run_length=20, event_table=etab)

test_that("fmri_lm accepts AR arguments", {
  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            cor_struct = "ar1", cor_iter = 2, cor_global = FALSE,
            ar1_exact_first = TRUE),
    NA
  )

  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            strategy = "chunkwise", nchunks = 2,
            cor_struct = "ar1", cor_iter = 1, cor_global = TRUE,
            ar1_exact_first = FALSE),
    NA
  )
})

test_that("ar_p validation works", {
  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            cor_struct = "arp"),
    "p must be specified in ar_options when struct is 'arp'"
  )

  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            cor_struct = "arp", ar_p = 2),
    NA
  )
})
