facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
facedes$repnum <- factor(facedes$rep_num)

test_that("can construct and run a simple fmri glm from in memory dataset", {
  
   scans <- lapply(1:length(unique(facedes$run)), function(i) {
     arr <- array(rnorm(10*10*10*244), c(10,10,10, 244))
     bspace <- neuroim2::NeuroSpace(dim=c(10,10,10,244))
     neuroim2::NeuroVec(arr, bspace)
   })
   
   mask <- neuroim2::LogicalNeuroVol(array(rnorm(10*10*10), c(10,10,10)) > 0, neuroim2::NeuroSpace(dim=c(10,10,10)))
   
   #scans <- list.files("test_data/images_study/epi/", "rscan0.*nii", full.names=TRUE)
   dset <- fmri_mem_dataset(scans=scans, 
                        mask=mask, 
                        TR=1.5, 
                        event_table=facedes)
   
   
   mod <- fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset=dset, durations=0)
   expect_true(!is.null(mod))
  
})

test_that("can construct and run a simple fmri glm from in memory dataset and one contrast", {
  
  scans <- lapply(1:length(unique(facedes$run)), function(i) {
    arr <- array(rnorm(10*10*10*244), c(10,10,10, 244))
    bspace <- neuroim2::NeuroSpace(dim=c(10,10,10,244))
    neuroim2::NeuroVec(arr, bspace)
  })
  
  mask <- neuroim2::LogicalNeuroVol(array(rnorm(10*10*10), c(10,10,10)) > 0, neuroim2::NeuroSpace(dim=c(10,10,10)))
  
  #scans <- list.files("test_data/images_study/epi/", "rscan0.*nii", full.names=TRUE)
  dset <- fmri_mem_dataset(scans=scans, 
                           mask=mask, 
                           TR=1.5, 
                           event_table=facedes)
  
  con <- contrast_set(pair_contrast( ~ repnum == 1, ~ repnum == 2, name="rep2_rep1"))
  
  mod <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0)
  expect_true(!is.null(mod))
  expect_equal(ncol(mod$result$contrasts$estimate()), 1)
  
})

# test_that("a one-run, one-contrast linear model analysis", {
#   df1 <- subset(imagedes,run==1)
#   df1 <- subset(df1, !is.na(onsetTime))
#   
#   df1$sdur <- scale(df1$duration)[,1]
#   
#   dmat <- matrix(rnorm(400*100), 400, 100)
#   md <- matrix_dataset(dmat, TR=1.5, run_length=400, event_table=df1)
#   con <- contrast_set(contrast( ~ Thorns - Massage, name="Thorns_Massage"))
#   mod <- fmri_lm(onsetTime ~ hrf(imageName, subset = !is.na(onsetTime), contrasts=con), ~ run, dataset=md, durations=sdur)
#  
# })


# test_that("a two-run, one contrast linear model analysis", {
#   df1 <- subset(imagedes,run %in% c(1,2))
#   df1 <- subset(df1, !is.na(onsetTime))
#   
#   df1$sdur <- scale(df1$duration)[,1]
#   
#   dmat <- matrix(rnorm(800*100), 800, 100)
#   md <- matrix_dataset(dmat, TR=1.5, run_length=c(400,400), event_table=df1)
#   con <- contrast_set(contrast( ~ Thorns - Massage, name="Thorns_Massage"))
#   mod <- fmri_lm(onsetTime ~ hrf(imageName, contrasts=con), ~ run, dataset=md, durations=sdur)
#   
#   
# })

# test_that("can load and run a simple config file", {
#   config <- read_fmri_config("test_data/images_study/config.R")
#   dset <- fmri_dataset(config$scans, config$mask, config$TR, 
#                            config$run_length,
#                            config$event_table,
#                            config$aux_table, 
#                            base_path=config$base_path)
#   
#   frame <- sampling_frame(dset$run_length, config$TR)
#   mod <- fmri_model(config$event_model, config$baseline_model, config$design, dset$aux_table,
#                     basis=HRF_SPMG1, dset$runids, dset$run_length, config$TR, drop_empty=TRUE)
#                     
#   
#   mod <- fmri_glm(config$event_model, 
#                   dataset=dset, durations=0)
# })



