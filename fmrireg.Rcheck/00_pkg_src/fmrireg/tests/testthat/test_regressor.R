options(mc.cores=1)
library(testthat)
library(assertthat)

facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
lopdes <- read.table(system.file("extdata", "design_aug.txt", package = "fmrireg"), header=TRUE)


test_that("regressor generates correct outputs", {
  onsets <- seq(1,100, by=5)
  durations=1
  
  reg1 <- regressor(onsets, HRF_GAMMA)

  expect_equal(class(reg1)[1],"Reg")
  expect_equal(length(reg1$onsets), length(onsets))
  expect_equal(length(reg1$duration), length(onsets))
  expect_true(durations(reg1)[1] == 0)
  expect_true(amplitudes(reg1)[1] == 1)
})

test_that("generate an event model with one observation per level", {
  sframe <- sampling_frame(blocklens=rep(401,5), TR=1.5)
  lopdes$onset <- lopdes$WordPresentationOnset/1000
  lopdes$Target <- factor(lopdes$Target)
  ev <- event_model(onset ~ hrf(Target), data=lopdes, block= ~ Run, sampling_frame=sframe)
  expect_true(!is.null(ev))
})

test_that("can construct a convolved term from an hrfspec with one factor and one run", {
  N <- 100
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  sframe <- sampling_frame(blocklens=100, TR=1)
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), Run=rep(1,4))
  ev <- event_model(onsets ~ hrf(fac), data=etab, block= ~ Run, sampling_frame=sframe)

  expect_equal(ncol(design_matrix(ev)), length(levels(etab$fac)))
})

test_that("can construct a convolved term from an hrfspec with two factors and one run", {
  N <- 100
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  sframe <- sampling_frame(blocklens=100, TR=1)
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), fac2=factor(letters[1:2]), Run=1)
  espec <- event_model(onsets ~ hrf(fac) + hrf(fac2), block= ~ Run, data=etab, sampling_frame=sframe)
  expect_equal(dim(design_matrix(espec)), c(N,length(levels(interaction(etab$fac, etab$fac2)))))
})

test_that("can construct a convolved term from an hrfspec with one factor and two runs", {
  N <- 100
  onsets1 <- seq(1,N,by=5)
  onsets2 <- seq(1,N,by=5)
  onsets <- c(onsets1,onsets2)
  durations <- 0
  
  sframe <- sampling_frame(blocklens=c(100,100), TR=1)
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), block=rep(1:2, c(length(onsets1), length(onsets2))))
  espec <- event_model(onsets ~ hrf(fac), data=etab, block = ~ block, sampling_frame=sframe)

  expect_equal(dim(design_matrix(espec)), c(N*2,length(levels(etab$fac))))
})


test_that("can construct a convolved trialwise term from an hrfspec with one factor and one run", {
  N <- 500
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  sframe <- sampling_frame(blocklens=500, TR=1)
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), block=1)
  espec <- event_model(onsets ~ trialwise(), data=etab, block=~block, sampling_frame=sframe)
 
  expect_equal(dim(design_matrix(espec)), c(N,length(onsets)))
})

test_that("can construct a convolved trialwise term with bspline basis from an hrfspec with one factor, one run", {
  N <- 500
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  sframe <- sampling_frame(blocklens=500, TR=1)
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), block=1)
  
  hrfb <- gen_hrf(hrf_bspline, N=5)
  espec <- event_model(onsets ~ trialwise(basis=hrfb), data=etab, block=~block, 
  sampling_frame=sframe)
  
  expect_equal(dim(design_matrix(espec)), c(N,length(onsets)*5))
})




test_that("can extract a design matrix from an fmri_model with two terms and one continuous variable", {
  N <- 100
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  sframe <- sampling_frame(blocklens=100, TR=1)
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), fac2=factor(letters[1:2]), z=rnorm(length(onsets)), run=1)
  
  espec <- event_model(onsets ~ hrf(fac,fac2) + hrf(z), data=etab, block = ~ run, sampling_frame=sframe)
  dmat <- design_matrix(espec)
  expect_equal(dim(dmat), c(N, 5))
})

test_that("can extract a design matrix from an fmri_model with one factor crossed with one continuous variable", {
  N <- 100
  onsets <- seq(1,N,by=10)
  sframe <- sampling_frame(blocklens=100, TR=1)
  
  etab <- data.frame(onsets=onsets, fac=factor(rep(c(1,2),5)), z=1:10, run=1)
  espec <- event_model(onsets ~ hrf(fac) + hrf(fac,z), etab, block=~run, sampling_frame=sframe)
  dmat <- design_matrix(espec)
  expect_equal(dim(dmat), c(N, length(levels(etab$fac)) * 2))
})

test_that("can extract a design matrix from an fmri_model with one factor and a 3rd degree bspline with 5 basis functions", {
  N <- 100
  onsets <- seq(1,N,by=10)
  durations <- 0
  sframe <- sampling_frame(blocklens=100, TR=1)

  etab <- data.frame(onsets=onsets, fac=factor(rep(c(1,2),5)), run=1)
  espec <- event_model(onsets ~ hrf(fac, basis="bs", nbasis=5), data=etab, block=~run, sampling_frame=sframe)
  dmat <- design_matrix(espec)
  expect_equal(dim(dmat), c(N, length(levels(etab$fac)) * 5))
})

test_that("can extract a design matrix from an fmri_model with one factor and SPMG3 basis", {
  N <- 100
  onsets <- seq(1,N,by=10)
  durations <- 0
  sframe <- sampling_frame(blocklens=100, TR=1)
  
  etab <- data.frame(onsets=onsets, fac=factor(rep(c(1,2),5)), run=1)
  espec <- event_model(onsets ~ hrf(fac, basis="spmg3"), data=etab, block=~run, sampling_frame=sframe)
  dmat <- design_matrix(espec)
  expect_equal(dim(dmat), c(N, length(levels(etab$fac)) * 3))
})


test_that("can extract a design matrix from an fmri_model with one trialwise factor and a 3rd degree bspline with 5 basis functions", {
  N <- 100
  onsets <- seq(1,N,by=10)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(rep(c(1,2),5)), run=1)
  sframe <- sampling_frame(blocklens=100, TR=1)
  
  espec <- event_model(onsets ~ trialwise(basis="bspline"), 
                       data=etab, block=~run, sampling_frame = sframe)
  dmat <- design_matrix(espec)
  expect_equal(dim(dmat), c(N, 4 * length(onsets)))
})

test_that("event_model with duplicate terms at different lags", {
  N <- 100
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  sframe <- sampling_frame(blocklens=100, TR=1)
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), Run=rep(1,4))
  ev <- event_model(onsets ~ hrf(fac) + hrf(fac,lag=5, prefix="phase2") , data=etab, block= ~ Run, sampling_frame=sframe)
  expect_true(!is.null(ev))
  
})

test_that("can construct a model with a raw covariate term", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  cdat <- data.frame(x = rnorm(sum(sframe$blocklens)), y = rnorm(sum(sframe$blocklens)))
  
  espec <- event_model(onset ~ hrf(rt) + covariate(x, y, data=cdat) , data=facedes, block=~run, sampling_frame=sframe)
  dmat <- design_matrix(espec)

  expect_equal(dim(dmat), c(sum(sframe$blocklens), 3))
  expect_equal(dmat$x, cdat$x)
  expect_equal(dmat$y, cdat$y)
  
  bmod <- baseline_model(basis="constant", sframe=sframe)
  fmod <- fmri_model(espec,bmod)
  #alm <- afni_lm(fmod)
  expect_true(!is.null(fmod))
})

test_that("covariate captures subset expression", {
  df <- data.frame(x = 1:4, y = 5:8, flag = c(TRUE, FALSE, TRUE, FALSE))
  cspec <- covariate(x, y, data = df, subset = flag)
  expect_true(identical(cspec$subset, rlang::expr(flag)))
})

test_that("facedes model with rep_num", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  espec <- event_model(onset ~ hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  dmat <- design_matrix(espec)
  expect_equal(dim(dmat), c(sum(sframe$blocklens), length(levels(facedes$repnum))))
})

test_that("facedes model with rep_num, subsetting rep_num == -1 ", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  
  espec <- event_model(onset ~ hrf(repnum, subset=repnum != "-1"), data=facedes, block=~run, sampling_frame=sframe)
  dmat <- design_matrix(espec)
  expect_equal(dim(dmat), c(sum(rep(436/2,max(facedes$run))), length(levels(facedes$repnum))-1))
})

test_that("facedes model with rep_num, and rep_num by rt ", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  
  espec <- event_model(onset ~ hrf(repnum, subset=repnum != "-1") + hrf(repnum,rt,subset=rep_num!= "-1"), 
                       data=facedes, block=~run,sampling_frame = sframe)
  
  dmat <- design_matrix(espec)
  expect_equal(dim(dmat), c(sum(rep(436/2,max(facedes$run))), (length(levels(facedes$repnum))-1)*2))
})

test_that("facedes model with RT and RT^2 as separate regressors ", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  
  facedes$rt_squared <- facedes$rt^2
  facedes$rt_cubed <- facedes$rt^3
  
  #conset <- contrast_set(contrast(~ rt - rt_sqared, "rt_diff"))
  espec <- event_model(onset ~ hrf(repnum) + hrf(Ident(rt, rt_squared)) + hrf(rt_cubed), 
                       data=facedes, block=~run,sampling_frame = sframe)
  
  dmat <- design_matrix(espec)
  expect_equal(dim(dmat), c(sum(rep(436/2,max(facedes$run))), (length(levels(facedes$repnum))-1)*2))
})




test_that("facedes model with polynomial parametric basis", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  
  espec <- event_model(onset ~  hrf(Poly(rt,3)), data=facedes, block=~run, sampling_frame=sframe)
  
  dmat <- design_matrix(espec)
  expect_equal(dim(dmat), c(sum(sframe$blocklens), 3))
})

test_that("facedes model with polynomial parametric basis and overriding onsets", {
  facedes$repnum <- factor(facedes$rep_num)
  facedes$onset2 <- facedes$onset+4
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  
  espec <- event_model(onset ~  hrf(Poly(rt,2), onsets=onset2), 
                       data=facedes, block=~run, sampling_frame=sframe)
  
  dmat <- design_matrix(espec)
  expect_equal(dim(dmat), c(sum(sframe$blocklens), 2))
})

test_that("facedes model with bspline parametric basis", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  
  espec <- event_model(onset ~  hrf(BSpline(rt,3)), data=facedes, block=~run, sampling_frame=sframe)
  
  dmat <- design_matrix(espec)
  expect_equal(dim(dmat), c(sum(sframe$blocklens), 3))
})

# test_that("rcan create an identity regressor", {
#   onsets <- seq(0,9,by=1)
#   durations=1
#   
#   reg1 <- regressor(onsets, HRF_IDENT, duration=0, amplitude=1:10)
#   
#   expect_equal(class(reg1)[1],"regressor")
#   expect_equal(length(reg1$onsets), length(onsets))
#   expect_equal(length(reg1$duration), length(onsets))
#   expect_true(durations(reg1)[1] == 0)
#   expect_true(amplitudes(reg1)[1] == 1)
# })

# Test null_regressor
test_that("null_regressor creates a valid object", {
  nr <- null_regressor()
  expect_is(nr, "Reg")
})

# Test single_trial_regressor
test_that("single_trial_regressor creates a valid object", {
  str <- single_trial_regressor(onsets = 10)
  expect_is(str, "Reg")
})

# Test regressor
test_that("regressor creates a valid object", {
  reg <- regressor(onsets = c(10, 12, 14, 16, 18, 40))
  expect_is(reg, "Reg")

})

# Test evaluate.single_trial_regressor
test_that("evaluate.single_trial_regressor produces correct output", {
  str <- single_trial_regressor(onsets = 10)
  grid <- seq(0, 100, by = 2)
  eval_str <- evaluate(str, grid)
  expect_is(eval_str, "numeric")
})

# Test evaluate.null_regressor
test_that("evaluate.null_regressor produces correct output", {
  nr <- null_regressor()
  grid <- seq(0, 100, by = 2)
  eval_nr <- evaluate(nr, grid)
  expect_is(eval_nr, "matrix")
})






# cset = list(
#   repnum == 1 ~ repnum == 2 | fac == 1,
#   repnum == 3 ~ repnum == 2
# )
# 
# cset <- contrast_set (
#    c1=contrast(
#      repnum == 1,
#      repnum == 2
#    ),
#    c2=contrast(
#      repnum==1
#    ),
#    c3=poly_contrast(
#      repnum,
#      value_map=list("-1" =0)
#    )
#  )





