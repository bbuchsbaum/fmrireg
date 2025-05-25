options(mc.cores=2)
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

test_that("can construct a baseline_model with 2 blocks and global intercept", {
  sframe <- sampling_frame(blocklens=c(100,100), TR=2)
  # Use poly basis to avoid degree restrictions of bs/ns
  bmodel <- baseline_model(basis="poly", degree=3, sframe=sframe, intercept="global")
  
  # Expect 3 poly terms per block + 1 global constant term
  expect_equal(ncol(design_matrix(bmodel)), 3*2 + 1)
  # Drift term + Block term (global constant)
  expect_equal(length(terms(bmodel)), 2) 
  # Check term names
  expect_equal(names(terms(bmodel)), c("drift", "block"))
})

test_that("can construct a baseline_model with 2 blocks and nuisance_list", {
  sframe <- sampling_frame(blocklens=c(100,100), TR=2)
  nlist <- list(
    matrix(rnorm(100*2), 100, 2), # 2 nuisance regressors block 1
    matrix(rnorm(100*3), 100, 3)  # 3 nuisance regressors block 2
  )
  
  bmodel <- baseline_model(basis="bs", degree=5, sframe=sframe, nuisance_list=nlist)
  
  # Expect 5 bspline terms per block + 1 constant term per block + 5 nuisance regressors
  # (5*2) + (1*2) + (2+3) = 10 + 2 + 5 = 17
  expect_equal(ncol(design_matrix(bmodel)), 17)
  # Drift term + Block term + Nuisance term
  expect_equal(length(terms(bmodel)), 3) 
  # Check that the expected terms are present, regardless of order
  expect_true(setequal(names(terms(bmodel)), c("block", "drift", "nuisance")))
})

test_that("can construct a baseline_model with basis='constant'", {
  sframe <- sampling_frame(blocklens=c(100,100), TR=2)
  # Intercept option should be ignored/overridden when basis is constant
  bmodel_runwise <- baseline_model(basis="constant", sframe=sframe, intercept="runwise")
  bmodel_global <- baseline_model(basis="constant", sframe=sframe, intercept="global")
  bmodel_none <- baseline_model(basis="constant", sframe=sframe, intercept="none")
  
  # Expect only 1 column (the constant baseline itself) regardless of intercept 
  # because the baseline IS the constant. Intercept logic is skipped.
  # Update: The current logic *does* add a separate block term if intercept != "none".
  # Let's test the current behavior. If basis='constant', drift term is 1 col per block.
  # If intercept='runwise', block term is 1 col per block. Total = 2*2 = 4? No, drift is global constant if intercept='global'.
  # If basis='constant' & intercept='runwise': drift term is 1 col/block (but value is 1?), block term is 1 col/block? -> This seems redundant.
  # Let's re-read baseline_model carefully for basis='constant'.
  # baseline(): sets degree=1. bfun returns matrix(rep(1, length(x))). 
  # construct.baselinespec(): calculates ret_list. Each element is matrix(1, block_len, 1). 
  # If intercept='global', returns a single column matrix of 1s. Term name 'baseline_constant_1'.
  # If intercept!='global', builds block-diagonal matrix of 1s (1 col per block). Term name 'baseline_constant_1'.
  # baseline_model(): 
  #   drift = construct(...) based on above.
  #   block = NULL because basis='constant'. 
  # So, only the 'drift' term is created.
  
  # Therefore: 
  # basis='constant', intercept='runwise' -> drift term = block-diagonal 1s (1 col per block). Total = 2 cols.
  expect_equal(ncol(design_matrix(bmodel_runwise)), 2)
  expect_equal(length(terms(bmodel_runwise)), 1)
  expect_equal(names(terms(bmodel_runwise)), c("drift"))
  
  # basis='constant', intercept='global' -> drift term = single column of 1s. Total = 1 col.
  expect_equal(ncol(design_matrix(bmodel_global)), 1)
  expect_equal(length(terms(bmodel_global)), 1)
  expect_equal(names(terms(bmodel_global)), c("drift"))
  
  # basis='constant', intercept='none' -> treated like 'runwise' by construct.baselinespec? No, intercept is passed.
  # construct.baselinespec doesn't explicitly use intercept='none' for constant basis differently than 'runwise'.
  # Let's assume it behaves like runwise. -> drift term = block-diagonal 1s (1 col per block). Total = 2 cols.
  expect_equal(ncol(design_matrix(bmodel_none)), 2)
  expect_equal(length(terms(bmodel_none)), 1)
  expect_equal(names(terms(bmodel_none)), c("drift"))

})