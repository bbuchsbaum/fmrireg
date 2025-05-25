options(mc.cores=1)
facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
facedes$repnum <- factor(facedes$rep_num)

library(testthat)
library(foreach)
library(ggrepel)

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


## test that latent and fmri_mem_dataset of same underlying latent dataset produce the same betas

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
   
   
   mod <- fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset=dset, durations=0, strategy="chunkwise", nchunks=4)
   expect_true(!is.null(mod))
  
   # Fast path version
   mod_fast <- fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset=dset, durations=0, strategy="chunkwise", nchunks=4, use_fast_path=TRUE)
   expect_true(!is.null(mod_fast))
   
   # Compare results
   expect_equal(coef(mod), coef(mod_fast), tolerance=1e-8)
   expect_equal(standard_error(mod), standard_error(mod_fast), tolerance=1e-8)
   expect_equal(stats(mod), stats(mod_fast), tolerance=1e-8)

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
  
  mod1 <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0, use_fast_path=TRUE)
  mod1a <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0)
  mod2 <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0, 
                  strategy="chunkwise", nchunks=10, verbose=FALSE)
  
  expect_true(!is.null(mod1))
  expect_true(!is.null(mod1a))
  expect_true(!is.null(mod2))
  
  expect_equal(ncol(stats(mod1, "contrasts")), 1)
  expect_equal(ncol(stats(mod1a, "contrasts")), 1)
  expect_equal(ncol(stats(mod2, "contrasts")), 1)
 
  # Fast path versions
  mod1_fast <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0, use_fast_path=TRUE)
  mod2_fast <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0, 
                       strategy="chunkwise", nchunks=10, verbose=FALSE, use_fast_path=TRUE)

  # Compare original runwise (mod1) vs fast runwise (mod1_fast)
  expect_equal(coef(mod1), coef(mod1_fast), tolerance=1e-8)
  expect_equal(standard_error(mod1), standard_error(mod1_fast), tolerance=1e-8)
  expect_equal(stats(mod1), stats(mod1_fast), tolerance=1e-8)
  expect_equal(coef(mod1, "contrasts"), coef(mod1_fast, "contrasts"), tolerance=1e-8)
  expect_equal(standard_error(mod1, "contrasts"), standard_error(mod1_fast, "contrasts"), tolerance=1e-8)
  expect_equal(stats(mod1, "contrasts"), stats(mod1_fast, "contrasts"), tolerance=1e-8)
  
  # Compare original chunkwise (mod2) vs fast chunkwise (mod2_fast)
  expect_equal(coef(mod2), coef(mod2_fast), tolerance=1e-8)
  expect_equal(standard_error(mod2), standard_error(mod2_fast), tolerance=1e-8)
  expect_equal(stats(mod2), stats(mod2_fast), tolerance=1e-8)
  expect_equal(coef(mod2, "contrasts"), coef(mod2_fast, "contrasts"), tolerance=1e-8)
  expect_equal(standard_error(mod2, "contrasts"), standard_error(mod2_fast, "contrasts"), tolerance=1e-8)
  expect_equal(stats(mod2, "contrasts"), stats(mod2_fast, "contrasts"), tolerance=1e-8)
 
})

test_that("can construct and run a simple fmri glm from a matrix_dataset with 1 column", {
  
  vals <- rep(rnorm(244),6)
  dset <- matrix_dataset(as.matrix(vals),TR=1.5, run_length=rep(244,6), event_table=facedes)
  
  c1 <- pair_contrast( ~ repnum == 1, ~ repnum == 2, name="rep2_rep1")
  c2 <- pair_contrast( ~ repnum == 3, ~ repnum == 4, name="rep3_rep4")
  con <<- contrast_set(c1,c2)
  
  mod1 <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0)
  mod2 <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0, 
                  strategy="chunkwise", nchunks=1)
  
  expect_true(!is.null(mod1))
  expect_equal(ncol(stats(mod1, "contrasts")), 2)
  expect_equal(ncol(stats(mod2, "contrasts")), 2)
  
  # Fast path versions
  mod1_fast <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0, use_fast_path=TRUE)
  mod2_fast <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0, 
                       strategy="chunkwise", nchunks=1, use_fast_path=TRUE)

  # Compare original runwise (mod1) vs fast runwise (mod1_fast)
  # Coefficients (betas) might differ slightly if only 1 voxel
  expect_equal(as.numeric(coef(mod1)), as.numeric(coef(mod1_fast)), tolerance=1e-8) 
  expect_equal(stats(mod1, "contrasts"), stats(mod1_fast, "contrasts"), tolerance=1e-8)
  expect_equal(standard_error(mod1, "contrasts"), standard_error(mod1_fast, "contrasts"), tolerance=1e-8)
  expect_equal(coef(mod1, "contrasts"), coef(mod1_fast, "contrasts"), tolerance=1e-8)

  # Compare original chunkwise (mod2) vs fast chunkwise (mod2_fast)
  expect_equal(as.numeric(coef(mod2)), as.numeric(coef(mod2_fast)), tolerance=1e-8)
  expect_equal(stats(mod2, "contrasts"), stats(mod2_fast, "contrasts"), tolerance=1e-8)
  expect_equal(standard_error(mod2, "contrasts"), standard_error(mod2_fast, "contrasts"), tolerance=1e-8)
  expect_equal(coef(mod2, "contrasts"), coef(mod2_fast, "contrasts"), tolerance=1e-8)
}) 

test_that("fmri glm for multivariate matrix and complex contrast ", {
  
  vals <- do.call(cbind, lapply(1:100, function(i) rnorm(244*6)))
  fd <- subset(facedes, null == 0 & rt < 2)
  fd$letter <- sample(factor(rep(letters[1:4], length.out=nrow(fd))))
  dset <- matrix_dataset(vals,TR=1.5, run_length=rep(244,6), event_table=fd)
  
  cset <<- contrast_set(pair_contrast( ~ letter %in% c("a", "b"), 
                       ~ letter %in% c("c", "d"),
                       name="abcd_efgh"),
                     pair_contrast( ~ letter %in% c("a", "c"), 
                       ~ letter %in% c("b", "d"),
                       name="ijkl_mnop"),
                     unit_contrast(~ letter, "letter"))
  
  #c3 <- unit_contrast(~ letter, "letter")

 
 # bmod <- baseline_model(basis="constant", degree=1, intercept="none", sframe=dset$sampling_frame)
  mod1 <- fmri_lm(onset ~ hrf(letter,  contrasts=cset), 
                  #baseline_model=bmod,
                  block = ~ run, dataset=dset, durations=0, nchunks=1,strategy="chunkwise")
  
  zz <- stats(mod1, "contrasts")
  expect_equal(ncol(zz),3)
  
  # Fast path version
  mod1_fast <- fmri_lm(onset ~ hrf(letter,  contrasts=cset), 
                       block = ~ run, dataset=dset, durations=0, nchunks=1, strategy="chunkwise", use_fast_path=TRUE)
  
  # Compare original chunkwise (mod1) vs fast chunkwise (mod1_fast)
  expect_equal(coef(mod1), coef(mod1_fast), tolerance=1e-8)
  expect_equal(stats(mod1, "contrasts"), stats(mod1_fast, "contrasts"), tolerance=1e-8)
  expect_equal(standard_error(mod1, "contrasts"), standard_error(mod1_fast, "contrasts"), tolerance=1e-8)
  expect_equal(coef(mod1, "contrasts"), coef(mod1_fast, "contrasts"), tolerance=1e-8)
  expect_true(!is.null(mod1))
  expect_true(!is.null(mod1_fast))
  
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
  expect_equal(ncol(stats(mod1, "contrasts")), 2)
  expect_equal(ncol(stats(mod2, "contrasts")), 2)

  # Fast path versions
  mod1_fast <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0, use_fast_path=TRUE)
  mod2_fast <- fmri_lm(onset ~ hrf(repnum,  contrasts=con), block = ~ run, dataset=dset, durations=0, 
                       strategy="chunkwise", nchunks=1, use_fast_path=TRUE)

  # Compare original runwise (mod1) vs fast runwise (mod1_fast)
  expect_equal(coef(mod1), coef(mod1_fast), tolerance=1e-8)
  expect_equal(stats(mod1, "contrasts"), stats(mod1_fast, "contrasts"), tolerance=1e-8)
  expect_equal(standard_error(mod1, "contrasts"), standard_error(mod1_fast, "contrasts"), tolerance=1e-8)
  expect_equal(coef(mod1, "contrasts"), coef(mod1_fast, "contrasts"), tolerance=1e-8)

  # Compare original chunkwise (mod2) vs fast chunkwise (mod2_fast)
  expect_equal(coef(mod2), coef(mod2_fast), tolerance=1e-8)
  expect_equal(stats(mod2, "contrasts"), stats(mod2_fast, "contrasts"), tolerance=1e-8)
  expect_equal(standard_error(mod2, "contrasts"), standard_error(mod2_fast, "contrasts"), tolerance=1e-8)
  expect_equal(coef(mod2, "contrasts"), coef(mod2_fast, "contrasts"), tolerance=1e-8)
  
})


test_that("can construct and run a simple fmri glm two terms and prefix args", {
  
  vals <- cbind(rep(rnorm(244),6), rep(rnorm(244),6))
  dset <- matrix_dataset(as.matrix(vals),TR=1.5, run_length=rep(244,6), event_table=facedes)
  
  
  mod1 <- fmri_lm(onset ~ hrf(repnum, subset=repnum %in% c(1,2), prefix="r12")+ 
                    hrf(repnum, subset=repnum %in% c(3,4), prefix="r34"),
                  block = ~ run, dataset=dset, durations=0)
 
  
  expect_true(!is.null(mod1))
  #expect_true(!is.null(mod2))
  expect_equal(ncol(stats(mod1)), 4)
  #expect_equal(ncol(mod2$result$contrasts$estimate()), 2)
  
  # Fast path version
  mod1_fast <- fmri_lm(onset ~ hrf(repnum, subset=repnum %in% c(1,2), prefix="r12")+ 
                         hrf(repnum, subset=repnum %in% c(3,4), prefix="r34"),
                       block = ~ run, dataset=dset, durations=0, use_fast_path=TRUE)
  
  expect_true(!is.null(mod1_fast))
  expect_equal(ncol(stats(mod1_fast)), 4)
  
  # Compare
  expect_equal(coef(mod1), coef(mod1_fast), tolerance=1e-8)
  expect_equal(stats(mod1), stats(mod1_fast), tolerance=1e-8)
  expect_equal(standard_error(mod1), standard_error(mod1_fast), tolerance=1e-8)
  
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
  

  # Fast path versions
  res1_fast <- fmrireg:::fmri_lm(Onset ~ hrf(Video, subset=Condition=="Encod", contrasts=conset) + 
                                   hrf(Video, subset=Condition=="Recall", prefix="rec"), block= ~ run, dataset=dset, 
                                  strategy="runwise", use_fast_path=TRUE)
  res2_fast <- fmrireg:::fmri_lm(Onset ~ hrf(Video, subset=Condition=="Encod", contrasts=conset) + 
                                    hrf(Video, subset=Condition=="Recall", prefix="rec"), block= ~ run, dataset=dset, 
                                  strategy="chunkwise", nchunks=12, use_fast_path=TRUE)
  res3_fast <- fmrireg:::fmri_lm(Onset ~ hrf(Video, subset=Condition=="Encod", contrasts=conset) + 
                              hrf(Video, subset=Condition=="Recall", prefix="rec"), block= ~ run, dataset=dset, 
                            strategy="chunkwise", nchunks=1, use_fast_path=TRUE)
  
  # Compare res1 vs res1_fast (runwise)
  expect_equal(coef(res1), coef(res1_fast), tolerance=1e-8)
  expect_equal(coef(res1, "contrasts"), coef(res1_fast, "contrasts"), tolerance=1e-8)
  expect_equal(stats(res1), stats(res1_fast), tolerance=1e-8)
  expect_equal(stats(res1, "contrasts"), stats(res1_fast, "contrasts"), tolerance=1e-8)
  expect_equal(standard_error(res1), standard_error(res1_fast), tolerance=1e-8)
  expect_equal(standard_error(res1, "contrasts"), standard_error(res1_fast, "contrasts"), tolerance=1e-8)
  
  # Compare res2 vs res2_fast (chunkwise, many chunks)
  expect_equal(coef(res2), coef(res2_fast), tolerance=1e-8)
  expect_equal(coef(res2, "contrasts"), coef(res2_fast, "contrasts"), tolerance=1e-8)
  expect_equal(stats(res2), stats(res2_fast), tolerance=1e-8)
  expect_equal(stats(res2, "contrasts"), stats(res2_fast, "contrasts"), tolerance=1e-8)
  expect_equal(standard_error(res2), standard_error(res2_fast), tolerance=1e-8)
  expect_equal(standard_error(res2, "contrasts"), standard_error(res2_fast, "contrasts"), tolerance=1e-8)

  # Compare res3 vs res3_fast (chunkwise, 1 chunk)
  expect_equal(coef(res3), coef(res3_fast), tolerance=1e-8)
  expect_equal(coef(res3, "contrasts"), coef(res3_fast, "contrasts"), tolerance=1e-8)
  expect_equal(stats(res3), stats(res3_fast), tolerance=1e-8)
  expect_equal(stats(res3, "contrasts"), stats(res3_fast, "contrasts"), tolerance=1e-8)
  expect_equal(standard_error(res3), standard_error(res3_fast), tolerance=1e-8)
  expect_equal(standard_error(res3, "contrasts"), standard_error(res3_fast, "contrasts"), tolerance=1e-8)

})

test_that("can run video fmri design with fmri_file_dataset", {
  library(neuroim2)
  des <- read.table(system.file("extdata", "video_design.txt", package = "fmrireg"), header=TRUE)
  events <- rep(320,7)
  sframe <- sampling_frame(rep(320, length(events)), TR=1.5)
  
  scans <- gen_fake_dataset(c(10,10,10,320), 7)
  maskfile <- gen_mask_file(c(10,10,10))
  
  dset <- fmri_dataset(scans, maskfile,TR=1.5, rep(320,7), base_path="/",mode="normal",  event_table=as_tibble(des))
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
                            strategy="chunkwise", nchunks=1)
  
  # Fast path versions
  res2_fast <- fmrireg:::fmri_lm(Onset ~ hrf(Video, subset=Condition=="Encod", contrasts=conset) + 
                              hrf(Video, subset=Condition=="Recall", prefix="rec"), block= ~ run, 
                            dataset=dset, 
                            strategy="chunkwise", nchunks=22, use_fast_path=TRUE)
  res3_fast <- fmrireg:::fmri_lm(Onset ~ hrf(Video, subset=Condition=="Encod", contrasts=conset) + 
                              hrf(Video, subset=Condition=="Recall", prefix="rec"), block= ~ run, dataset=dset, 
                            strategy="chunkwise", nchunks=1, use_fast_path=TRUE)
                            
  expect_true(!is.null(coef(res2)))
  expect_true(!is.null(coef(res3)))
  
  expect_true(!is.null(coef(res2_fast)))
  expect_true(!is.null(coef(res3_fast)))
  
  expect_true(!is.null(coef(res2, "contrasts")))
  expect_true(!is.null(coef(res3, "contrasts")))
  
  # Compare res2 vs res2_fast (chunkwise, many chunks)
  expect_equal(coef(res2), coef(res2_fast), tolerance=1e-8)
  expect_equal(coef(res2, "contrasts"), coef(res2_fast, "contrasts"), tolerance=1e-8)
  expect_equal(stats(res2), stats(res2_fast), tolerance=1e-8)
  expect_equal(stats(res2, "contrasts"), stats(res2_fast, "contrasts"), tolerance=1e-8)
  expect_equal(standard_error(res2), standard_error(res2_fast), tolerance=1e-8)
  expect_equal(standard_error(res2, "contrasts"), standard_error(res2_fast, "contrasts"), tolerance=1e-8)

  # Compare res3 vs res3_fast (chunkwise, 1 chunk)
  expect_equal(coef(res3), coef(res3_fast), tolerance=1e-8)
  expect_equal(coef(res3, "contrasts"), coef(res3_fast, "contrasts"), tolerance=1e-8)
  expect_equal(stats(res3), stats(res3_fast), tolerance=1e-8)
  expect_equal(stats(res3, "contrasts"), stats(res3_fast, "contrasts"), tolerance=1e-8)
  expect_equal(standard_error(res3), standard_error(res3_fast), tolerance=1e-8)
  expect_equal(standard_error(res3, "contrasts"), standard_error(res3_fast, "contrasts"), tolerance=1e-8)
  
  # Clean up temporary files
  unlink(scans)
  unlink(maskfile)
})

test_that("can run video fmri design with latent_dataset", {
  #library(multivarious)
  des <- read.table(system.file("extdata", "video_design.txt", package = "fmrireg"), header=TRUE)
  events <- rep(320,7)
  sframe <- sampling_frame(rep(320, length(events)), TR=1.5)
  
  scans <- gen_fake_dataset(c(10,10,10,320), 7)
  vecs <- lapply(scans, read_vec)
  maskfile <- gen_mask_file(c(10,10,10))
  mask <- read_vol(maskfile)
  
  mats <- lapply(vecs, function(v) series(v, mask!=0))
  mat <- do.call(rbind, mats)
  pres <- multivarious::pca(mat, ncomp=488, preproc=multivarious::pass())
  lvec <- fmristore::LatentNeuroVec(pres$s, pres$v, add_dim(space(mask), nrow(mat)), 
                                   mask=mask)
  ldset <- latent_dataset(lvec, 1.5, run_length=rep(320,7), des)
  
  evmod <- event_model(Onset ~ hrf(Video, Condition, basis="spmg1"), 
                       block = ~ run, sampling_frame=sframe, data=des)
   
  conset <<- NULL
  conset <<- do.call("contrast_set", lapply(levels(factor(des$Video)), function(v) {
    f1 <- as.formula(paste("~ Video == ", paste0('"', v, '"')))
    f2 <- as.formula(paste("~ Video != ", paste0('"', v, '"')))
    pair_contrast(f1, f2, name=paste0(v, "_vsall"))
  }))
  
  conset2 <<- do.call("contrast_set", lapply(levels(factor(des$Video)), function(v) {
    f1 <- as.formula(paste("~ rec_Video == ", paste0('"', v, '"')))
    f2 <- as.formula(paste("~ rec_Video != ", paste0('"', v, '"')))
    pair_contrast(f1, f2, name=paste0("rec_", v, "_vsall"))
  }))
  
  
  # Note: fmri_latent_lm does not currently have the fast path implementation
  res2 <- fmrireg:::fmri_latent_lm(Onset ~ hrf(Video, subset=Condition=="Encod", contrasts=conset) + 
                              hrf(Video, subset=Condition=="Recall", prefix="rec", contrasts=conset2), 
                              block= ~ run, 
                              autocor="none", dataset=ldset)
  
  se1 <- standard_error(res2, "contrasts", recon=TRUE)
  con1 <- stats.fmri_latent_lm(res2, "contrasts", recon=TRUE)
  
  # Run standard fmri_lm for comparison (if needed, but latent path is different)
  # dset <- fmri_dataset(scans, maskfile,TR=1.5, rep(320,7), base_path="/", event_table=des)
  # res3 <- fmrireg:::fmri_lm(Onset ~ hrf(Video, subset=Condition=="Encod", contrasts=conset) + 
  #                                    hrf(Video, subset=Condition=="Recall", prefix="rec"), block= ~ run, 
  #                                  strategy="chunkwise", nchunks=1, dataset=dset)
  # se2 <- standard_error(res3, "contrasts")
  # con2 <- stats(res3, "contrasts")
  
  expect_true(!is.null(se1))
  expect_true(!is.null(con1))
  
  # Clean up temporary files
  unlink(scans)
  unlink(maskfile)
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



