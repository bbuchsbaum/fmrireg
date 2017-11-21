test_that("can construct an fmri_dataset", {
  dset <- fmri_dataset(
    scans=c("scan1.nii", "scan2.nii", "scan3.nii"),
    mask="mask.nii",
    run_length=c(100,100,100),
    TR=2
  )
  expect_true(!is.null(dset))
  
})

test_that("can read a config file to create fmri_dataset", {
  fname <- system.file("extdata", "config.R", package = "fmrireg")
  base_path=dirname(fname)
  
  config <- read_fmri_config(fname, base_path)
  expect_true(!is.null(config))
})