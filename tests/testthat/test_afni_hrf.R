facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)


test_that("can construct an simple afni native stimulus model", {
  
  facedes$repnum <- factor(facedes$rep_num)
  scans <- paste0("rscan0", 1:6, ".nii")
  dset <- fmri_dataset(scans=scans,
                       mask="mask.nii",
                       TR=1.5,
                       run_length=rep(436,6),
                       event_table=facedes)
  
  
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
