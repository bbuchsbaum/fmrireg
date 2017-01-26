facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)

imagedes <- read.table("test_data/images_study/behavior/design.txt", header=TRUE)
image_mask_file = "test_data/images_study/epi/global_mask.nii"

test_that("can construct a simple fmri glm", {
 

  scans <- list.files("test_data/images_study/epi/", "rscan0.*nii", full.names=TRUE)
  dset <- fmri_dataset(scans=scans, 
                       mask=image_mask_file, 
                       TR=1.5, 
                       run_length=rep(348,6), 
                       event_table=imagedes)
  
  con <- contrast_set(contrast( ~ Thorns - Massage, name="Thorns_Massage"))
  mod <- fmri_glm(onsetTime ~ hrf(imageName, subset = !is.na(onsetTime), contrasts=con), ~ run, dataset=dset, durations=0)
 
  
})

test_that("can construct an fmri_dataset from test_data", {
  df1 <- subset(imagedes,run==1)
  df1 <- subset(df1, !is.na(onsetTime))
  
  df1$sdur <- scale(df1$duration)[,1]
  
  scans <- "test_data/images_study/epi/rscan01.nii"
  mask <- "test_data/images_study/epi/global_mask.nii"
  dset <- fmri_dataset(scans=scans, mask=mask, TR=2, event_table=df1, blocklens=rep(348,1), blockids=rep(1,nrow(df1)))
  #iter <- data_chunks(dset, nchunks=5)
  con <- contrast(imageName=="Shark", imageName == "Tools")
  con2 <- contrast(imageName == "Candy", imageName == "Massage")
  con3 <- contrast(sdur == "sdur")
  
  mod <- fmri_glm(onsetTime ~ hrf(imageName, subset=!is.na(imageName), 
                                  contrasts=contrast_set(con=con, con2=con2)) + 
                              hrf(sdur, contrasts=contrast(sdur== "sdur")) + block(run),
                              dataset=dset, durations=0)
  
})

test_that("can load and run a simple config file", {
  config <- read_fmri_config("test_data/images_study/config.R")
  dset <- fmri_dataset(config$scans, config$mask, config$TR, 
                           config$run_length,
                           config$event_table,
                           config$aux_table, 
                           base_path=config$base_path)
  
  frame <- sampling_frame(dset$run_length, config$TR)
  mod <- fmri_model(config$event_model, config$baseline_model, config$design, dset$aux_table,
                    basis=HRF_SPMG1, dset$runids, dset$run_length, config$TR, drop_empty=TRUE)
                    
  
  mod <- fmri_glm(config$event_model, 
                  dataset=dset, durations=0)
})



