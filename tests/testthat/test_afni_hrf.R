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
  
  # Check that design_matrix(fmod) does not error, but gives warnings due to AFNI term
  # and that the resulting matrix effectively only contains the baseline model columns
  dm_baseline <- design_matrix(bspec)
  dm_fmod <- design_matrix(fmod)
  expect_equal(ncol(dm_fmod), ncol(dm_baseline))
  expect_equal(colnames(dm_fmod), colnames(dm_baseline))
  
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
  
  # Check design matrix behavior for mixed model
  dm_baseline <- design_matrix(bspec)
  # Expect columns from hrf(constant) + baseline. AFNI term should warn and not add columns.
  # Find the non-AFNI term (hrf(constant)) specifically
  event_terms <- terms(espec)
  hrf_constant_term <- NULL
  for (term_name in names(event_terms)) {
    term <- event_terms[[term_name]]
    if (!inherits(term, c("afni_hrf_convolved_term", "afni_trialwise_convolved_term"))) {
      hrf_constant_term <- term
      break
    }
  }
  
  expect_true(!is.null(hrf_constant_term), "Should find at least one non-AFNI term")
  dm_hrf_constant <- design_matrix(hrf_constant_term, sampling_frame=dset$sampling_frame)
  expected_event_cols <- ncol(dm_hrf_constant)

  dm_fmod_mixed <- design_matrix(fmod)
  expect_equal(ncol(dm_fmod_mixed), expected_event_cols + ncol(dm_baseline))
  
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
  
  dm_baseline <- design_matrix(bspec)
  # Find the non-AFNI term (hrf(Poly...)) specifically
  event_terms <- terms(espec)
  hrf_poly_term <- NULL
  for (term_name in names(event_terms)) {
    term <- event_terms[[term_name]]
    if (!inherits(term, c("afni_hrf_convolved_term", "afni_trialwise_convolved_term"))) {
      hrf_poly_term <- term
      break
    }
  }
  
  expect_true(!is.null(hrf_poly_term), "Should find at least one non-AFNI term")
  dm_hrf_poly <- design_matrix(hrf_poly_term, sampling_frame=dset$sampling_frame)
  expected_event_cols_poly <- ncol(dm_hrf_poly)
  
  dm_fmod_poly <- design_matrix(fmod)
  expect_equal(ncol(dm_fmod_poly), expected_event_cols_poly + ncol(dm_baseline))
  
  expect_equal(2, length(baseline_terms(fmod)))
  expect_null(contrast_weights(fmod)$repnum)
})





