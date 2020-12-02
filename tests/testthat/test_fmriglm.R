facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
facedes$repnum <- factor(facedes$rep_num)


gen_mask_file <- function(d, perc) {
  arr = array(0,d)
  vals <- ifelse(runif(prod(d)) > .5, 1, 0)
  vol <- NeuroVol(vals, NeuroSpace(d))
  fname <- paste0(tempfile(), ".nii")
  write_vol(vol, fname)
  fname
}

gen_fake_dataset <- function(d, nscans) {
  
  onames <- vector(length=nscans, mode="list")
  for (i in 1:nscans) {
    arr <- array(rnorm(prod(d)), d)
    bspace <- neuroim2::NeuroSpace(dim=d)
    vec <- neuroim2::NeuroVec(arr, bspace)
    fname <- paste0(tempfile(), ".nii")
    write_vec(vec, fname)
    onames[i] <- fname
  }
  onames
}

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
  
  con <<- contrast_set(pair_contrast( ~ repnum == 1, ~ repnum == 2, name="rep2_rep1"))
  
  mod1 <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0)
  mod2 <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0, 
                  strategy="chunkwise", nchunks=10)
  expect_true(!is.null(mod1))
  expect_true(!is.null(mod2))
  expect_equal(ncol(mod1$result$contrasts$estimate()), 1)
  expect_equal(ncol(mod2$result$contrasts$estimate()), 1)
  
})

test_that("can construct and run a simple fmri glm from a matrix_dataset with 1 column", {
  
  vals <- rep(rnorm(244),6)
  dset <- matrix_dataset(as.matrix(vals),TR=1.5, run_length=rep(244,6), event_table=facedes)
  
  c1 <- pair_contrast( ~ repnum == 1, ~ repnum == 2, name="rep2_rep1")
  c2 <- pair_contrast( ~ repnum == 2, ~ repnum == 3, name="rep3_rep2")
  con <<- contrast_set(c1,c2)
  
  mod1 <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0)
  mod2 <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0, 
                  strategy="chunkwise", nchunks=1)
  
  expect_true(!is.null(mod1))
  expect_equal(ncol(mod1$result$contrasts$estimate()), 2)
  expect_equal(ncol(mod2$result$contrasts$estimate()), 2)
  
})

test_that("can construct and run a simple fmri glm from a matrix_dataset with 2 columns", {
  
  vals <- cbind(rep(rnorm(244),6), rep(rnorm(244),6))
  dset <- matrix_dataset(as.matrix(vals),TR=1.5, run_length=rep(244,6), event_table=facedes)
  
  c1 <- pair_contrast( ~ repnum == 1, ~ repnum == 2, name="rep2_rep1")
  c2 <- pair_contrast( ~ repnum == 2, ~ repnum == 3, name="rep3_rep2")
  con <<- contrast_set(c1,c2)
  
  mod1 <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0)
  mod2 <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0, 
                  strategy="chunkwise", nchunks=1)
  
  expect_true(!is.null(mod1))
  expect_true(!is.null(mod2))
  expect_equal(ncol(mod1$result$contrasts$estimate()), 2)
  expect_equal(ncol(mod2$result$contrasts$estimate()), 2)
  
})


test_that("can construct and run a simple fmri glm two terms and prefix args", {
  
  vals <- cbind(rep(rnorm(244),6), rep(rnorm(244),6))
  dset <- matrix_dataset(as.matrix(vals),TR=1.5, run_length=rep(244,6), event_table=facedes)
  

  mod1 <- fmri_lm(onset ~ hrf(repnum, subset=repnum %in% c(1,2), prefix="r12")+ 
                          hrf(repnum, subset=repnum %in% c(3,4), prefix="r34"),
                  block = ~ run, dataset=dset, durations=0)
 
  
  expect_true(!is.null(mod1))
  #expect_true(!is.null(mod2))
  expect_equal(ncol(mod1$result$betas$estimate()), 4)
  #expect_equal(ncol(mod2$result$contrasts$estimate()), 2)
  
})


test_that("can run video fmri design with matrix_dataset", {
  des <- read.table(system.file("extdata", "video_design.txt", package = "fmrireg"), header=TRUE)
  events <- rep(320,7)
  sframe <- sampling_frame(rep(320, length(events)), TR=1.5)
  
  evmod <- event_model(Onset ~ hrf(Video, Condition, basis="spmg1"), 
                       block = ~ run, sampling_frame=sframe, data=des)
  bmod <- baseline_model(basis="bs", degree=4, sframe=sframe)
  fmod <- fmri_model(evmod, bmod)
  
  dset <- matrix_dataset(matrix(rnorm(320*7*100), 320*7, 100),TR=1.5, run_length=rep(320,7), event_table=des)

  #conset <- fmrireg::one_against_all_contrast(levels(des$Video), "Video")
  
  conset <<- do.call("contrast_set", lapply(levels(factor(des$Video)), function(v) {
    f1 <- as.formula(paste("~ Video == ", paste0('"', v, '"')))
    f2 <- as.formula(paste("~ Video != ", paste0('"', v, '"')))
    pair_contrast(f1, f2, name=paste0(v, "_vsall"))
  }))
  
  res1 <- fmrireg:::fmri_lm(Onset ~ hrf(Video, subset=Condition=="Encod", contrasts=conset) + 
                                   hrf(Video, subset=Condition=="Recall", prefix="rec"), block= ~ run, dataset=dset, 
                                  strategy="runwise")
  res2 <- fmrireg:::fmri_lm(Onset ~ hrf(Video, subset=Condition=="Encod", contrasts=conset) + 
                                    hrf(Video, subset=Condition=="Recall", prefix="rec"), block= ~ run, dataset=dset, 
                            strategy="chunkwise", nchunks=12)
  
  res3 <- fmrireg:::fmri_lm(Onset ~ hrf(Video, subset=Condition=="Encod", contrasts=conset) + 
                              hrf(Video, subset=Condition=="Recall", prefix="rec"), block= ~ run, dataset=dset, 
                            strategy="chunkwise", nchunks=1)
  
  expect_true(!is.null(coef(res1)))
  expect_true(!is.null(coef(res2)))
  expect_true(!is.null(coef(res3)))
  
  expect_true(!is.null(coef(res1, "contrasts")))
  expect_true(!is.null(coef(res2, "contrasts")))
  expect_true(!is.null(coef(res3, "contrasts")))
  

})

test_that("can run video fmri design with fmri_file_dataset", {
  library(neuroim2)
  des <- read.table(system.file("extdata", "video_design.txt", package = "fmrireg"), header=TRUE)
  events <- rep(320,7)
  sframe <- sampling_frame(rep(320, length(events)), TR=1.5)
  
  scans <- gen_fake_dataset(c(10,10,10,320), 7)
  maskfile <- gen_mask_file(c(10,10,10))
  
  dset <- fmri_dataset(scans, maskfile,TR=1.5, rep(320,7), base_path="/", mode="bigvec", event_table=des)
  evmod <- event_model(Onset ~ hrf(Video, Condition, basis="spmg1"), 
                       block = ~ run, sampling_frame=sframe, data=des)
  bmod <- baseline_model(basis="bs", degree=4, sframe=sframe)
  fmod <- fmri_model(evmod, bmod)

  conset <<- NULL
  conset <<- do.call("contrast_set", lapply(levels(factor(des$Video)), function(v) {
    f1 <- as.formula(paste("~ Video == ", paste0('"', v, '"')))
    f2 <- as.formula(paste("~ Video != ", paste0('"', v, '"')))
    pair_contrast(f1, f2, name=paste0(v, "_vsall"))
  }))
  
 
  res2 <- fmrireg:::fmri_lm(Onset ~ hrf(Video, subset=Condition=="Encod", contrasts=conset) + 
                              hrf(Video, subset=Condition=="Recall", prefix="rec"), block= ~ run, 
                            dataset=dset, 
                            strategy="chunkwise", nchunks=22)
  
  res3 <- fmrireg:::fmri_lm(Onset ~ hrf(Video, subset=Condition=="Encod", contrasts=conset) + 
                              hrf(Video, subset=Condition=="Recall", prefix="rec"), block= ~ run, dataset=dset, 
                            strategy="chunkwise", nchunks=2)
  
  expect_true(!is.null(coef(res2)))
  expect_true(!is.null(coef(res3)))
  
  expect_true(!is.null(coef(res2, "contrasts")))
  expect_true(!is.null(coef(res3, "contrasts")))
  
  
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



