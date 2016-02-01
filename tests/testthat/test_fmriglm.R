facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)

test_that("can construct a simple fmri glm", {
 
  dset <- fmri_dataset(scans=paste0(letters[1:6], ".nii"), mask=NULL, TR=2, blocklens=rep(226,6), blockids=facedes$run, event_table=facedes)
  mod <- fmri_glm(onset ~ hrf(rep_num), dataset=dset, durations=0)
  cmod <- construct(mod$model_spec)
  
})

test_that("can construct an fmri_dataset from test_data", {
  scans <- "test_data/images_study/epi/rscan01.nii"
  mask <- "test_data/images_study/epi/global_mask.nii"
  dset <- fmri_dataset(scans=scans, mask=mask, TR=2, blocklens=rep(348,1), blockids=1)
  iter <- data_chunks(dset)
  
})

