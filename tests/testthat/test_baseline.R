
test_that("can construct a baseline_model with 1 block and 5th order bspline basis", {
  sframe <- sampling_frame(blocklens=100, TR=2)
  bmodel <- baseline_model(basis="bs", degree=5, sframe=sframe)
  
  ## 5 + 1 for constant term
  expect_equal(ncol(design_matrix(bmodel)),6)
  expect_equal(length(terms(bmodel)), 2)
})

test_that("can construct a baseline_model with 1 block and 5th order bspline basis and nuisance_list", {
  sframe <- sampling_frame(blocklens=100, TR=2)
  nlist <- list(
    poly(samples(sframe), 3)
  )
  
  bmodel <- baseline_model(basis="bs", degree=5, sframe=sframe, nuisance_list=nlist)
  
  ## 5 + 1 for constant term
  expect_equal(ncol(design_matrix(bmodel)),9)
  expect_equal(length(terms(bmodel)), 3)
})

test_that("can construct a baseline_model with 2 blocks and 5th order bspline basis", {
  sframe <- sampling_frame(blocklens=c(100,100), TR=2)
  bmodel <- baseline_model(basis="bs", degree=5, sframe=sframe)
  
  expect_equal(ncol(design_matrix(bmodel)),6*2)
  expect_equal(length(terms(bmodel)), 2)
  #p <- print(bmodel)
  #expect_true(TRUE)
  
})