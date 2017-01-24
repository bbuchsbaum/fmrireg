library(testthat)
facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)

test_that("regressor generates correct outputs", {
  onsets <- seq(1,100, by=5)
  durations=1
  
  reg1 <- regressor(onsets, HRF_GAMMA)

  expect_equal(class(reg1)[1],"regressor")
  expect_equal(length(reg1$onsets), length(onsets))
  expect_equal(length(reg1$duration), length(onsets))
})


test_that("can construct a convolved term from an hrfspec with one factor and one run", {
  N <- 100
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)))
  mspec <- fmri_model(onsets ~ hrf(fac), etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=100, TR=1)
  
  expect_equal(dim(mspec$terms[[1]]$design_matrix), c(N,length(levels(etab$fac))))
})

test_that("can construct a convolved term from an hrfspec with two factors and one run", {
  N <- 100
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), fac2=letters[1:2])
  mspec <- fmri_model(onsets ~ hrf(fac) + hrf(fac2), etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=100, TR=1)
 
  expect_equal(dim(design_matrix(mspec)), c(N,length(levels(interaction(etab$fac, etab$fac2)))))
})

test_that("can construct a convolved term from an hrfspec with one factor and two runs", {
  N <- 100
  onsets1 <- seq(1,N,by=5)
  onsets2 <- seq(1,N,by=5)
  onsets <- c(onsets1,onsets2)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), block=rep(1:2, c(length(onsets1), length(onsets2))))
  mspec <- fmri_model(onsets ~ hrf(fac), etab, durations=1, blockids=etab$block, blocklens=c(100,100), TR=1)

  
  expect_equal(dim(design_matrix(mspec)), c(N*2,length(levels(etab$fac))))
})


test_that("can construct a convolved trialwise term from an hrfspec with one factor and one run", {
  N <- 500
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)))
  mspec <- fmri_model(onsets ~ fac, etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=N, TR=1)
  hspec <- trialwise(fac)
  
  cterm <- construct(hspec, mspec)
  expect_equal(dim(cterm$design_matrix), c(N,length(onsets)))
})


test_that("can extract a design matrix from an fmri_model with one term", {
  N <- 100
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)))
  mspec <- fmri_model(onsets ~ hrf(fac), etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=N, TR=1)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(N, length(levels(etab$fac))))
})

test_that("can extract a design matrix from an fmri_model with two terms", {
  N <- 100
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), fac2=factor(letters[1:2]))
  mspec <- fmri_model(onsets ~ hrf(fac,fac2), etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=N, TR=1)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(N, length(levels(interaction(etab$fac, etab$fac2)))))
})

test_that("can extract a design matrix from an fmri_model with two terms and one continuous variable", {
  N <- 100
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), fac2=factor(letters[1:2]), z=rnorm(length(onsets)))
  mspec <- fmri_model(onsets ~ hrf(fac,fac2) + hrf(z), etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=N, TR=1)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(N, 5))
})

test_that("can extract a design matrix from an fmri_model with one factor crossed with one continuous variable", {
  N <- 100
  onsets <- seq(1,N,by=10)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(rep(c(1,2),5)), z=1:10)
  mspec <- fmri_model(onsets ~ hrf(fac) + hrf(fac,z), etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=N, TR=1)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(N, length(levels(etab$fac)) * 2))
})

test_that("can extract a design matrix from an fmri_model with one factor and a 3rd degree bspline with 5 basis functions", {
  N <- 100
  onsets <- seq(1,N,by=10)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(rep(c(1,2),5)))
  mspec <- fmri_model(onsets ~ hrf(fac, basis="bspline", nbasis=5), etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=N, TR=1)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(N, length(levels(etab$fac)) * 5))
})

test_that("can extract a design matrix from an fmri_model with one factor and SPMG3 basis", {
  N <- 100
  onsets <- seq(1,N,by=10)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(rep(c(1,2),5)))
  mspec <- fmri_model(onsets ~ hrf(fac, basis="spmg3"), etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=N, TR=1)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(N, length(levels(etab$fac)) * 3))
})


test_that("can extract a design matrix from an fmri_model with one trialwise factor and a 3rd degree bspline with 5 basis functions", {
  N <- 100
  onsets <- seq(1,N,by=10)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(rep(c(1,2),5)))
  mspec <- fmri_model(onsets ~ trialwise(fac, basis="bspline", nbasis=5), etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=N, TR=1)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(N, 5 * length(onsets)))
})

test_that("facedes model with rep_num", {
  facedes$repnum <- factor(facedes$rep_num)
  mspec <- fmri_model(onset ~ hrf(repnum), facedes, durations=0, blockids=facedes$run, blocklens=rep(436/2,max(facedes$run)), TR=2)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(sum(rep(436/2,max(facedes$run))), length(levels(facedes$repnum))))
})

test_that("facedes model with rep_num, subsetting rep_num == -1 ", {
  facedes$repnum <- factor(facedes$rep_num)
  mspec <- fmri_model(onset ~ hrf(repnum, subset=repnum != "-1"), facedes, durations=0, blockids=facedes$run, blocklens=rep(436/2,max(facedes$run)), TR=2)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(sum(rep(436/2,max(facedes$run))), length(levels(facedes$repnum))-1))
})

test_that("facedes model with rep_num, and rep_num by rt ", {
  facedes$repnum <- factor(facedes$rep_num)
  mspec <- fmri_model(onset ~ hrf(repnum, subset=repnum != "-1") + hrf(repnum,rt,subset=rep_num!= "-1"), facedes, durations=0, blockids=facedes$run, blocklens=rep(436/2,max(facedes$run)), TR=2)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(sum(rep(436/2,max(facedes$run))), (length(levels(facedes$repnum))-1)*2))
})

test_that("facedes model block variable", {
  facedes$repnum <- factor(facedes$rep_num)
  
  aux_data <- data.frame(run=rep(1:6, each=218))
  mspec <- fmri_model(onset ~  hrf(repnum) + block(run), facedes, durations=0, blockids=facedes$run, 
                      blocklens=rep(436/2,max(facedes$run)), TR=2, aux_data=aux_data)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(sum(rep(436/2,max(facedes$run))), 11))
})

test_that("facedes model with polynomial parametric basis", {
  facedes$repnum <- factor(facedes$rep_num)
  
  aux_table <- data.frame(run=rep(1:6, each=218))
  mspec <- fmri_model(onset ~  hrf(Poly(rt,3)) + block(run), facedes, durations=0, blockids=facedes$run, 
                      blocklens=rep(436/2,max(facedes$run)), TR=2, aux_data=aux_table)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(sum(rep(436/2,max(facedes$run))), 9))
})

test_that("facedes model with bspline baseline term and repnum regressor", {
  facedes$repnum <- factor(facedes$rep_num)
  
  aux_table <- data.frame(run=rep(1:6, each=218))
  mspec <- fmri_model(onset ~  hrf(repnum) + baseline(degree=3, basis="bs"), facedes, durations=0, blockids=facedes$run, 
                      blocklens=rep(436/2,max(facedes$run)), TR=2)
  dmat <- design_matrix(mspec)
  expect_equal(dim(dmat), c(sum(rep(436/2,max(facedes$run))), 9))
})

test_that("can build a simple contrast from a convolved term", {
  facedes$repnum <- factor(facedes$rep_num)
  aux_table <- data.frame(run=rep(1:6, each=218))
  mspec <- fmri_model(onset ~  hrf(repnum) + baseline(degree=3, basis="bs"), facedes, durations=0, blockids=facedes$run, 
                      blocklens=rep(436/2,max(facedes$run)), TR=2)
  
  term <- construct(mspec$varspec[[1]], mspec)
  con <- contrast(A=repnum==-1, B=repnum==1)
  expect_equal(as.vector(contrast_weights(con, term)), c(1,-1,0,0,0))
})

test_that("can build a contrast versus the intercept from a convolved term", {
  facedes$repnum <- factor(facedes$rep_num)
  aux_table <- data.frame(run=rep(1:6, each=218))
  mspec <- fmri_model(onset ~  hrf(repnum) + baseline(degree=3, basis="bs"), facedes, durations=0, blockids=facedes$run, 
                      blocklens=rep(436/2,max(facedes$run)), TR=2)
  
  term <- construct(mspec$varspec[[1]], mspec)
  con <- contrast(A=repnum==-1)
  expect_equal(as.vector(contrast_weights(con, term)), c(1,0,0,0,0))
})

test_that("can build a contrast versus the intercept and add to hrfspec", {
  facedes$repnum <- factor(facedes$rep_num)
  aux_table <- data.frame(run=rep(1:6, each=218))
  
  
  
  conf <- contrast_formula(~ `2` - !`1`, id="repnum")
  
  sframe <- sampling_frame(rep(436/2,max(facedes$run)), TR=2)
  
  nuisance <- matrix(rnorm(2*length(sframe$blockids)), length(sframe$blockids), 2)
  
  bm <- baseline_model(basis="bs", degree=3, sampling_frame=sframe, nuisance_matrix=nuisance)
  em <- event_model(onset ~ hrf(repnum, contrasts=con), block = ~ run, data=facedes, sampling_frame=sframe)
  mod <- fmri_model(em, bm)
  
  term <- construct(mspec$varspec[[1]], mspec)
  expect_equal(as.vector(contrast_weights(con, term)), c(1,0,0,0,0))
})

test_that("can build a linear contrast from repnum and value_map", {
  facedes$repnum <- factor(facedes$rep_num)
  aux_table <- data.frame(run=rep(1:6, each=218))
  con <- poly_contrast(A=repnum, value_map=list("-1"=0, "1"=1, "2"=2, "3"=3, "4"=4))
  
  sframe <- sampling_frame(rep(436/2,max(facedes$run)), TR=2)
  mspec <- fmri_model(onset ~  hrf(repnum, contrasts=con), ~ baseline(degree=3, basis="bs"), 
                      event_table=facedes, event_block_ids=facedes$run, sampling_frame=sframe)
  
  term <- construct(mspec$varspec[[1]], mspec)
  expect_equal(as.vector(contrast_weights(con, term)), as.vector(poly(c(0,1,2,3,4))))
})




# cset = list(
#   repnum == 1 ~ repnum == 2 | fac == 1,
#   repnum == 3 ~ repnum == 2
# )
# 
cset <- contrast_set (
   c1=contrast(
     repnum == 1,
     repnum == 2
   ),
   c2=contrast(
     repnum==1
   ),
   c3=poly_contrast(
     repnum,
     value_map=list("-1" =0)
   )
 )





