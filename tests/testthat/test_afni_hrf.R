facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
library(assertthat)
options(mc.cores=2)

# Helper function to assign a dummy mask to a dataset object
# This prevents errors when tests use dummy file paths for masks/scans
# with fmri_file_dataset, as some operations might try to access dset$mask.
.assign_dummy_mask_to_dataset <- function(dset) {
  dummy_mask_space <- neuroim2::NeuroSpace(dim = c(5L, 5L, 5L)) # Arbitrary small dimensions
  dset$mask <- neuroim2::NeuroVol(array(1L, dim = dim(dummy_mask_space)), dummy_mask_space)
  return(dset) # Return the modified dataset
}

context("afni")

test_that("can construct an simple afni native stimulus model", {
  
  facedes$repnum <- factor(facedes$rep_num)
  scans <- paste0("rscan0", 1:6, ".nii")
  dset <- fmri_dataset(scans=scans,
                       mask="mask.nii",
                       TR=1.5,
                       run_length=rep(436,6),
                       event_table=facedes)
  
  dset <- .assign_dummy_mask_to_dataset(dset)
  
  espec <- event_model(onset ~ afni_hrf(repnum, basis="csplin", nbasis=8, start=0, stop=18), data=facedes, block=~run, sampling_frame=dset$sampling_frame)
  bspec <- baseline_model(basis="bs", degree=5, sframe=dset$sampling_frame)
  fmod <- fmri_model(espec, bspec)
  alm <- afni_lm(fmod, dset)
  
  expect_true(!is.null(fmod))
  expect_equal(length(terms(fmod)), 3)
  expect_error(design_matrix(fmod))
  expect_equal(2, length(baseline_terms(fmod)))
  expect_null(contrast_weights(fmod)$repnum)
 
})

test_that("can construct an an afni model with trialwise regressor", {
  
  facedes$repnum <- factor(facedes$rep_num)
  facedes$constant <- rep(1, nrow(facedes))
  scans <- paste0("rscan0", 1:6, ".nii")
  dset <- fmri_dataset(scans=scans,
                       mask="mask.nii",
                       TR=1.5,
                       run_length=rep(436,6),
                       event_table=facedes)
  
  dset <- .assign_dummy_mask_to_dataset(dset)
  
  espec <- event_model(onset ~ hrf(constant, basis="spmg1") + afni_trialwise("trial", basis="spmg1"), data=facedes, block=~run, sampling_frame=dset$sampling_frame)
  bspec <- baseline_model(basis="bs", degree=5, sframe=dset$sampling_frame)
  fmod <- fmri_model(espec, bspec)
  alm <- afni_lm(fmod, dset)
  
  expect_true(!is.null(fmod))
  expect_equal(length(terms(fmod)), 4)
  expect_error(design_matrix(fmod))
  expect_equal(2, length(baseline_terms(fmod)))
  expect_null(contrast_weights(fmod)$repnum)
  
})

test_that("can construct an an afni model with a constant", {
  
  facedes$repnum <- factor(facedes$rep_num)
  facedes$constant <- factor(rep(1, nrow(facedes)))
  scans <- paste0("rscan0", 1:6, ".nii")
  dset <- fmri_dataset(scans=scans,
                       mask="mask.nii",
                       TR=1.5,
                       run_length=rep(436,6),
                       event_table=facedes)
  
  dset <- .assign_dummy_mask_to_dataset(dset)
  
  espec <- event_model(onset ~ afni_hrf(constant, basis="spmg1"), data=facedes, 
                       block=~run, sampling_frame=dset$sampling_frame)
  bspec <- baseline_model(basis="bs", degree=5, sframe=dset$sampling_frame)
  fmod <- fmri_model(espec, bspec)
  alm <- afni_lm(fmod, dset)
  
  expect_true(!is.null(fmod))
 
  
})

test_that("can construct an an afni model with trialwise regressor and a Polynomial modulator", {
  
  facedes$repnum <- as.numeric(as.character(facedes$rep_num))
  scans <- paste0("rscan0", 1:6, ".nii")
  dset <- fmri_dataset(scans=scans,
                       mask="mask.nii",
                       TR=1.5,
                       run_length=rep(436,6),
                       event_table=facedes)
  
  dset <- .assign_dummy_mask_to_dataset(dset)
  
  espec <- event_model(onset ~ afni_trialwise("trial") + hrf(fmrireg::Poly(repnum,2)), 
                      data=facedes, 
                      block=~run, 
                      sampling_frame=dset$sampling_frame)
  bspec <- baseline_model(basis="bs", degree=5, sframe=dset$sampling_frame)
  fmod <- fmri_model(espec, bspec)
  alm <- afni_lm(fmod, dset)
  
  expect_true(!is.null(fmod))
  expect_equal(length(terms(fmod)), 4)
  expect_error(design_matrix(fmod))
  expect_equal(2, length(baseline_terms(fmod)))
  expect_null(contrast_weights(fmod)$repnum)
})

test_that("afni_hrf accepts precision, summate, and lag", {

  facedes$repnum <- factor(facedes$rep_num)
  scans <- paste0("rscan0", 1:6, ".nii")
  dset <- fmri_dataset(scans=scans,
                       mask="mask.nii",
                       TR=1.5,
                       run_length=rep(436,6),
                       event_table=facedes)

  dset <- .assign_dummy_mask_to_dataset(dset)

  espec <- event_model(onset ~ afni_hrf(repnum, basis="spmg1", precision=0.2, summate=FALSE, lag=2),
                       data=facedes,
                       block=~run,
                       sampling_frame=dset$sampling_frame)

  term <- terms(espec)[[1]]
  expect_equal(term$hrfspec$precision, 0.2)
  expect_false(term$hrfspec$summate)
  expect_equal(term$hrfspec$lag, 2)
})

test_that("afni_trialwise accepts precision and summate", {

  scans <- paste0("rscan0", 1:6, ".nii")
  dset <- fmri_dataset(scans=scans,
                       mask="mask.nii",
                       TR=1.5,
                       run_length=rep(436,6),
                       event_table=facedes)

  dset <- .assign_dummy_mask_to_dataset(dset)

  espec <- event_model(onset ~ afni_trialwise("trial", precision=0.2, summate=FALSE),
                       data=facedes,
                       block=~run,
                       sampling_frame=dset$sampling_frame)

  term <- terms(espec)[[1]]
  expect_equal(term$hrfspec$precision, 0.2)
  expect_false(term$hrfspec$summate)
})

