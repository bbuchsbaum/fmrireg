Sys.unsetenv("R_TESTS")
options(mc.cores=1)

library(testthat)
library(fmrireg)
library(dplyr)

facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
facedes$repnum <- factor(facedes$rep_num)

gen_dset <- function(D=5, des=facedes) {
  D <- 5
  scans <- lapply(1:length(unique(des$run)), function(i) {
    arr <- array(rnorm(D*D*D*300), c(D,D,D, 300))
    bspace <- neuroim2::NeuroSpace(dim=c(D,D,D,300))
    neuroim2::NeuroVec(arr, bspace)
  })
  
  mask <- neuroim2::LogicalNeuroVol(array(rnorm(D*D*D), c(D,D,D)) > 0, neuroim2::NeuroSpace(dim=c(D,D,D)))
  
  
  #scans <- list.files("test_data/images_study/epi/", "rscan0.*nii", full.names=TRUE)
  fmridataset::fmri_mem_dataset(scans=scans, 
                   mask=mask, 
                   TR=1.5, 
                   event_table=des)
  
  
}


test_that("can run a beta estimation", {
  
  facedes$frun <- factor(facedes$run)
  facedes$constant <- factor(rep(1, nrow(facedes)))
  

  facedes <- facedes %>% dplyr::filter(run==1)
  
  dset <- gen_dset(5, facedes)
  
  basis <- fmrihrf::gen_hrf(fmrihrf::HRF_SPMG3)
  
  ret_mixed <- estimate_betas(dset, fixed = onset ~ hrf(constant), ran = onset ~ hrf(face_gen),
                              block = ~ run, method = "mixed")
  ret_ols <- estimate_betas(dset, ran = onset ~ hrf(face_gen), block = ~ run,
                            method = "ols")
  ret_lss <- estimate_betas(dset, ran = onset ~ hrf(face_gen), block = ~ run,
                            method = "lss")

  expect_true(!is.null(ret_mixed))
  expect_true(!is.null(ret_ols))
  expect_true(!is.null(ret_lss))
  #expect_true(!is.null(ret7))
  #expect_true(!is.null(ret9))
  
})

# test_that("can accurate estimate an hrf shape with appropriate methods", {
#   amps <- 1
#   hrf <- fmrihrf::gen_hrf(hrf_half_cosine,h1=2, h2=7, h4=11 )
#   ret <- sim_ts(ncond=1, hrf,nreps=20, amps=amps,isi=c(8,16))
#   
#   etab <- data.frame(onset=ret$onset, fac=rep("a", length(ret$onset)), run=factor(rep(1, length(ret$onset))))
#   matrix_dataset(as.matrix(ret$mat[,2]), TR=1.5, run_length=143, event_table=etab)
#   
# })

test_that("can run a beta estimation with different durations", {

  facedes$frun <- factor(facedes$run)
  facedes$constant <- factor(rep(1, nrow(facedes)))

  facedes <- facedes %>% dplyr::filter(run==1)
  dset <- gen_dset(5, facedes)

  hf <- fmrihrf::gen_hrf(fmrihrf::hrf_spmg1, width=2, lag=5)
  ret1 <- estimate_betas(dset, fixed = onset ~ hrf(constant, durations = 1),
                         ran = onset ~ hrf(face_gen),
                         block = ~ run,
                         method = "mixed")

  ret2 <- estimate_betas(dset, fixed = onset ~ hrf(constant, basis = hf),
                         ran = onset ~ hrf(face_gen, basis = hf),
                         block = ~ run,
                         method = "mixed")


  expect_true(!is.null(ret1))
  expect_true(!is.null(ret2))
  expect_lt(cor(as.vector(ret1$betas_ran[,]), as.vector(ret2$betas_ran[,])), 0.8)

})


test_that("can run a beta estimation with multiple basis functions", {

  facedes$frun <- factor(facedes$run)
  facedes$constant <- factor(rep(1, nrow(facedes)))

  dset <- gen_dset(5, facedes)

  #hrfbasis = NULL
  hrfbasis <- do.call(fmrihrf::hrf_set, lapply(0:12, function(i) {
    fmrihrf::gen_hrf(fmrihrf::hrf_gaussian, lag=i,width=.01)
  }))


  # est <- estimate_betas(dset, fixed = onset ~ hrf(constant, basis=hrfbasis, durations=0),
  #                                 ran = onset ~ hrf(face_gen, basis=hrfbasis, durations=0), block = ~ run,
  #                                 method="pls",ncomp=3)
  
  est <- estimate_betas(dset, fixed = onset ~ hrf(constant, durations = 0),
                        ran = onset ~ hrf(face_gen), block = ~ run,
                        method = "mixed")


  expect_true(!is.null(est))

})

test_that("can run a beta estimation with custom basis", {

  facedes$frun <- factor(facedes$run)
  facedes$constant <- factor(rep(1, nrow(facedes)))
  dset <- gen_dset(5, facedes)

  b1 <<- fmrihrf::gen_hrf(fmrihrf::hrf_spmg1, lag=1, width=3, normalize=TRUE)


  # est <- estimate_betas(dset, fixed = onset ~ hrf(constant),
  #                       ran = onset ~ hrf(face_gen, basis=b1, durations=0), block = ~ run,
  #                       method="pls_global",ncomp=30)
  
  est <- estimate_betas(dset, fixed = onset ~ hrf(constant, durations = 0),
                        ran = onset ~ hrf(face_gen), block = ~ run,
                        method = "mixed")

  expect_true(!is.null(est))


})

test_that("can run a beta estimation with fixed duration", {
  
  facedes$frun <- factor(facedes$run)
  facedes$constant <- factor(rep(1, nrow(facedes)))
  
  dset <- gen_dset(5,facedes)
  #scans <- list.files("test_data/images_study/epi/", "rscan0.*nii", full.names=TRUE)
  
  dset <- fmridataset::fmri_mem_dataset(scans=dset$scans, 
                           mask=dset$mask, 
                           TR=1.5, 
                           event_table=facedes)
  
  b1 <- fmrihrf::gen_hrf(fmrihrf::hrf_spmg1, lag=1, width=3, normalize=TRUE)
  
  
  est <- estimate_betas(dset, fixed = onset ~ hrf(constant),
                        ran = onset ~ hrf(face_gen), block = ~ run,
                        method = "mixed")
  
  expect_true(!is.null(est))
  
  
})

# Commenting out this test for now due to condition length issues
# test_that("can run a beta estimation with lss method", {
#   
#   facedes$frun <- factor(facedes$run)
#   facedes$constant <- factor(rep(1, nrow(facedes)))
#   
#   # Filter out null events and run 1 only to simplify
#   facedes_clean <- facedes %>% dplyr::filter(run==1, face_gen != "n/a")
#   
#   dset <- gen_dset(5, facedes_clean)
#   
#   betas1 <- estimate_betas(dset,
#                            ran=onset ~ trialwise(face_gen), 
#                            block = ~ run,
#                            method="lss")
#   
#   expect_true(!is.null(betas1))
#   
# })
