Sys.unsetenv("R_TESTS")
options(mc.cores=1)

facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
facedes$repnum <- factor(facedes$rep_num)
library(dplyr)


gen_dset <- function(D=5, des=facedes) {
  D <- 5
  scans <- lapply(1:length(unique(des$run)), function(i) {
    arr <- array(rnorm(D*D*D*300), c(D,D,D, 300))
    bspace <- neuroim2::NeuroSpace(dim=c(D,D,D,300))
    neuroim2::NeuroVec(arr, bspace)
  })
  
  mask <- neuroim2::LogicalNeuroVol(array(rnorm(D*D*D), c(D,D,D)) > 0, neuroim2::NeuroSpace(dim=c(D,D,D)))
  
  
  #scans <- list.files("test_data/images_study/epi/", "rscan0.*nii", full.names=TRUE)
  fmri_mem_dataset(scans=scans, 
                   mask=mask, 
                   TR=1.5, 
                   event_table=des)
  
  
}


test_that("can run a beta estimation", {
  
  facedes$frun <- factor(facedes$run)
  facedes$constant <- factor(rep(1, nrow(facedes)))
  

  facedes <- facedes %>% dplyr::filter(run==1)
  
  dset <- gen_dset(5, facedes)
  
  basis <- gen_hrf(HRF_SPMG3)
  
  ret1 <- estimate_betas(dset, fixed = onset ~ hrf(constant), ran = onset ~ trialwise(), block = ~ run, 
                       method="pls", ncomp=1)
  ret2 <- estimate_betas(dset, fixed = onset ~ hrf(constant), ran = onset ~ trialwise(), block = ~ run, 
                         method="pls", ncomp=3)
  ret3 <- estimate_betas(dset, fixed = onset ~ hrf(constant), ran = onset ~ trialwise(), block = ~ run, 
                        method="mixed")
  
  ret4 <- estimate_betas(dset, ran = onset ~ trialwise(), block = ~ run, 
                         method="ols")
  ret5 <- estimate_betas(dset, ran = onset ~ trialwise(), block = ~ run, 
                         method="lss")
  ret6 <- estimate_betas(dset, ran = onset ~ trialwise(), block = ~ run, 
                         method="lss_cpp")
  
  #ret7 <- estimate_betas(dset, ran = onset ~ trialwise(), block = ~ run, 
  #                       method="lss_naive")
  
  hrf_basis <- basis(0:25)
  hrf_reference <- HRF_SPMG1(0:25)
  ret8 <- estimate_betas(dset, ran = onset ~ trialwise(basis=basis), block = ~ run, 
                         method="r1", hrf_basis=hrf_basis, hrf_ref=hrf_reference)
  
  ret8 <- estimate_betas(dset, ran = onset ~ trialwise(), block = ~ run, 
                         method="lowrank_hrf")
  

  
 
  expect_true(!is.null(ret1))
  expect_true(!is.null(ret2))
  expect_true(!is.null(ret3))
  expect_true(!is.null(ret4))
  expect_true(!is.null(ret5))
  expect_true(!is.null(ret6))
  expect_true(!is.null(ret7))
  expect_true(!is.null(ret8))
  
})

test_that("can accurate estimate an hrf shape with appropriate methods", {
  amps <- runif(20)*2
  hrf <- gen_hrf(hrf_half_cosine,h1=2, h2=7, h4=11 )
  ret <- sim_ts(ncond=1, hrf,nreps=20, amps=amps,isi=16)
  
  etab <- data.frame(onset=ret$onset, fac=rep("a", length(ret$onset)), run=factor(rep(1, length(ret$onset))))
  matrix_dataset(as.matrix(ret$mat[,2]), TR=1.5, run_length=143, event_table=etab)
  
})

test_that("can run a beta estimation with different durations", {

  facedes$frun <- factor(facedes$run)
  facedes$constant <- factor(rep(1, nrow(facedes)))

  facedes <- facedes %>% dplyr::filter(run==1)
  dset <- gen_dset(5, facedes)

  hf <- gen_hrf(hrf_spmg1, width=10)
  ret1 <- estimate_betas(dset, fixed = onset ~ hrf(constant, durations=3), ran = onset ~ trialwise(durations=5),
                         block = ~ run,
                         method="pls", ncomp=1)


  expect_true(!is.null(ret1))



})


test_that("can run a beta estimation with multiple basis functions", {

  facedes$frun <- factor(facedes$run)
  facedes$constant <- factor(rep(1, nrow(facedes)))

  dset <- gen_dset(5, facedes)

  #hrfbasis = NULL
  hrfbasis <- do.call(gen_hrf_set, lapply(0:12, function(i) {
    gen_hrf(hrf_gaussian, lag=i,width=.01)
  }))


  # est <- estimate_betas(dset, fixed = onset ~ hrf(constant, basis=hrfbasis, durations=0),
  #                                 ran = onset ~ trialwise(basis=hrfbasis, durations=0), block = ~ run,
  #                                 method="pls",ncomp=3)
  
  est <- estimate_betas(dset, fixed = onset ~ hrf(constant, durations=0),
                                   ran = onset ~ trialwise(durations=0), block = ~ run,
                                   method="pls",ncomp=3)


  expect_true(!is.null(est))

})

test_that("can run a beta estimation with custom basis", {

  facedes$frun <- factor(facedes$run)
  facedes$constant <- factor(rep(1, nrow(facedes)))
  dset <- gen_dset(5, facedes)

  b1 <<- gen_hrf(hrf_spmg1, lag=1, width=3, normalize=TRUE)


  # est <- estimate_betas(dset, fixed = onset ~ hrf(constant),
  #                       ran = onset ~ trialwise(basis=b1, durations=0), block = ~ run,
  #                       method="pls_global",ncomp=30)
  
  est <- estimate_betas(dset, fixed = onset ~ hrf(constant,durations=0),
                                   ran = onset ~ trialwise(durations=0), block = ~ run,
                                   method="pls",ncomp=3)

  expect_true(!is.null(est))


})

# 
# test_that("can run a beta estimation with fixed duration", {
#   
#   facedes$frun <- factor(facedes$run)
#   facedes$constant <- factor(rep(1, nrow(facedes)))
#   
#   dset <- gen_dset(5,facedes)
#   #scans <- list.files("test_data/images_study/epi/", "rscan0.*nii", full.names=TRUE)
#   
#   dset <- fmri_mem_dataset(scans=dset$scans, 
#                            mask=dset$mask, 
#                            TR=1.5, 
#                            event_table=facedes)
#   
#   b1 <<- gen_hrf(hrf_spmg1, lag=1, width=3, normalize=TRUE)
#   
#   
#   est <- estimate_betas(dset, fixed = onset ~ hrf(constant),
#                         ran = onset ~ trialwise(durations=4), block = ~ run, 
#                         method="pls_global",ncomp=20)
#   
#   expect_true(!is.null(est))
#   
#   
# })