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
  mspec <- fmri_model(onsets ~ fac, etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=100, TR=1)
  hspec <- hrf(fac)
  
  cterm <- construct(hspec, mspec)
  expect_equal(dim(cterm$convolved_term), c(N,length(levels(etab$fac))))
})

test_that("can construct a convolved term from an hrfspec with two factors and one run", {
  N <- 100
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), fac2=letters[1:2])
  mspec <- fmri_model(onsets ~ fac + fac2, etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=100, TR=1)
  hspec <- hrf(fac,fac2)
  
  cterm <- construct(hspec, mspec)
  expect_equal(dim(cterm$convolved_term), c(N,length(levels(interaction(etab$fac, etab$fac2)))))
})

test_that("can construct a convolved term from an hrfspec with one factor and two runs", {
  N <- 100
  onsets1 <- seq(1,N,by=5)
  onsets2 <- seq(1,N,by=5)
  onsets <- c(onsets1,onsets2)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)), block=rep(1:2, c(length(onsets1), length(onsets2))))
  mspec <- fmri_model(onsets ~ fac, etab, durations=1, blockids=etab$block, blocklens=c(100,100), TR=1)
  hspec <- hrf(fac)
  
  cterm <- construct(hspec, mspec)
  expect_equal(dim(cterm$convolved_term), c(N*2,length(levels(etab$fac))))
})


test_that("can construct a convolved trialwise term from an hrfspec with one factor and one run", {
  N <- 100
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)))
  mspec <- fmri_model(onsets ~ fac, etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=N, TR=1)
  hspec <- trialwise(fac)
  
  cterm <- construct(hspec, mspec)
  expect_equal(dim(cterm$convolved_term), c(N,length(onsets)))
})


test_that("can extract a design matrix from an fmri_model with one term", {
  N <- 100
  onsets <- seq(1,N,by=5)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)))
  mspec <- fmri_model(onsets ~ hrf(fac), etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=N, TR=1)
  dmat <- design_matrix(mspec)

})


