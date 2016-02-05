facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)

imagedes <- read.table("test_data/images_study/behavior/design.txt", header=TRUE)

test_that("can construct a simple fmri glm", {
 
  dset <- fmri_dataset(scans=paste0(letters[1:6], ".nii"), mask=NULL, TR=2, blocklens=rep(226,6), blockids=facedes$run, event_table=facedes)
  mod <- fmri_glm(onset ~ hrf(rep_num), dataset=dset, durations=0)
  cmod <- construct(mod$model_spec)
  
})

test_that("can construct an fmri_dataset from test_data", {
  df1 <- subset(imagedes,run==1)
  df1 <- subset(df1, !is.na(onsetTime))
  scans <- "test_data/images_study/epi/rscan01.nii"
  mask <- "test_data/images_study/epi/global_mask.nii"
  dset <- fmri_dataset(scans=scans, mask=mask, TR=2, event_table=df1, blocklens=rep(348,1), blockids=rep(1,nrow(df1)))
  #iter <- data_chunks(dset, nchunks=5)
  con <- contrast(imageName=="Shark", imageName == "Tools")
  mod <- fmri_glm(onsetTime ~ hrf(imageName, subset=!is.na(imageName), contrasts=con), dataset=dset, durations=0)
  
})

