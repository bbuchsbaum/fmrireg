

facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
facedes$repnum <- factor(facedes$rep_num)

test_that("can run a beta estimation", {
  
  facedes$frun <- factor(facedes$run)
  facedes$constant <- factor(rep(1, nrow(facedes)))
  
  scans <- lapply(1:length(unique(facedes$run)), function(i) {
    arr <- array(rnorm(10*10*10*300), c(10,10,10, 300))
    bspace <- neuroim2::NeuroSpace(dim=c(10,10,10,300))
    neuroim2::NeuroVec(arr, bspace)
  })
  
  mask <- neuroim2::LogicalNeuroVol(array(rnorm(10*10*10), c(10,10,10)) > 0, neuroim2::NeuroSpace(dim=c(10,10,10)))
  
  #scans <- list.files("test_data/images_study/epi/", "rscan0.*nii", full.names=TRUE)
  dset <- fmri_mem_dataset(scans=scans, 
                           mask=mask, 
                           TR=1.5, 
                           event_table=facedes)
  
  
 ret1 <- estimate_betas(dset, fixed = onset ~ hrf(constant), ran = onset ~ trialwise(), block = ~ run, 
                       method="pls", ncomp=1)
 ret2 <- estimate_betas(dset, fixed = onset ~ hrf(constant), ran = onset ~ trialwise(), block = ~ run, 
                         method="pls", ncomp=3)
 ret3 <- estimate_betas(dset, fixed = onset ~ hrf(constant), ran = onset ~ trialwise(), block = ~ run, 
                        method="slm")
 ret4 <- estimate_betas(dset, fixed = onset ~ hrf(constant), ran = onset ~ trialwise(), block = ~ run, 
                        method="mixed")
 ret5 <- estimate_betas(dset, fixed = onset ~ hrf(constant), ran = onset ~ trialwise(), block = ~ run, 
                        method="pls_searchlight", niter=3)
 
 expect_true(!is.null(ret1))
 expect_true(!is.null(ret2))
 expect_true(!is.null(ret3))
 expect_true(!is.null(ret4))
  
})