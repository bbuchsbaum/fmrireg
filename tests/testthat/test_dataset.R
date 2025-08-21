test_that("can construct an fmri_dataset", {
  # Use dummy_mode to test without actual files
  dset <- fmridataset::fmri_dataset(
    scans = c("dummy1.nii", "dummy2.nii"),
    mask = "dummy_mask.nii",
    TR = 2,
    run_length = c(100, 100),
    dummy_mode = TRUE
  )
  
  expect_s3_class(dset, "fmri_dataset")
  expect_equal(length(dset$sampling_frame$blocklens), 2)
  expect_equal(dset$sampling_frame$blocklens, c(100, 100))
  # TR is stored per run, so check all are equal to 2
  expect_equal(dset$sampling_frame$TR, c(2, 2))
})


## design file not found during testing
#test_that("can read a config file to create fmri_dataset", {
  #fname <- system.file("extdata", "config.R", package = "fmrireg")
  #base_path=dirname(fname)
  
  #config <- read_fmri_config(fname, base_path)
  #expect_true(!is.null(config))
#})

test_that("can construct an fmri_mem_dataset", {
  
  facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
  facedes$repnum <- factor(facedes$rep_num)
  
  scans <- lapply(1:length(unique(facedes$run)), function(i) {
    arr <- array(rnorm(10*10*10*244), c(10,10,10, 244))
    bspace <- neuroim2::NeuroSpace(dim=c(10,10,10,244))
    neuroim2::NeuroVec(arr, bspace)
  })
  
  mask <- neuroim2::LogicalNeuroVol(array(rnorm(10*10*10), c(10,10,10)) > 0, neuroim2::NeuroSpace(dim=c(10,10,10)))
  
  #scans <- list.files("test_data/images_study/epi/", "rscan0.*nii", full.names=TRUE)
  dset <- fmridataset::fmri_mem_dataset(scans=scans, 
                           mask=mask, 
                           TR=1.5, 
                           event_table=tibble::as_tibble(facedes))
  
  expect_true(!is.null(dset))
  
  
})

