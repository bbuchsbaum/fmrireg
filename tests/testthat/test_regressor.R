facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)

test_that("regressor generates correct outputs", {
  onsets <- seq(1,100, by=5)
  durations=1
  
  reg1 <- regressor(onsets, HRF_GAMMA)

  expect_equal(class(reg1)[1],"regressor")
  expect_equal(length(reg1$onsets), length(onsets))
  expect_equal(length(reg1$duration), length(onsets))
})


test_that("can constucrt a convolved term from an hrfspec with one factor and one run", {
  onsets <- seq(1,100,by=5)
  durations <- 0
  
  etab <- data.frame(onsets=onsets, fac=factor(c(1,1,2,2)))
  mspec <- fmri_model(onsets ~ fac, etab, durations=1, blockids=rep(1,nrow(etab)), blocklens=100, TR=1)
  hspec <- hrf(fac)
  
  cterm <- construct(hspec, mspec)
  cterm$convolved_term
})