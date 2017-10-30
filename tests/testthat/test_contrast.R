library(testthat)


test_that("a 2-by-2 Fcontrast", {
  F_a <- factor(rep(letters[1:2], 8))
  F_b <- factor(rep(c("V1", "V2"), each=8))
  onsets <- seq(1, length(F_a))
  blockids <- rep(1, length(onsets))
  
  et <- event_term(list(Fa=F_a, Fb=F_b), onsets, blockids)
  Fcontrasts(et)
})

test_that("a 3-by-2 Fcontrast", {
  F_a <- factor(rep(letters[1:2], 8))
  F_b <- factor(rep(c("V1", "V2"), each=8))
  F_c <- factor(rep(c("B3", "B3", "B4", "B4"), 4))
  onsets <- seq(1, length(F_a))
  blockids <- rep(1, length(onsets))
  
  et <- event_term(list(Fa=F_a, Fb=F_b, Fc=F_c), onsets, blockids)
  Fcontrasts(et)
})

test_that("a 3-by-3 Fcontrast", {
  F_a <- factor(sample(letters[1:3], 100, replace=TRUE))
  F_b <- factor(sample(c("V1", "V2", "V3"), 100, replace=TRUE))
  F_c <- factor(sample(c("B3", "B3", "B4", "B4", "B5", "B5"),100, replace=TRUE))
  onsets <- seq(1, length(F_a))
  blockids <- rep(1, length(onsets))
  
  et <- event_term(list(Fa=F_a, Fb=F_b, Fc=F_c), onsets, blockids)
  Fcontrasts(et)
})


# test_that("can build a simple contrast from a convolved term", {
#   facedes$repnum <- factor(facedes$rep_num)
#   aux_table <- data.frame(run=rep(1:6, each=218))
#   mspec <- fmri_model(onset ~  hrf(repnum) + baseline(degree=3, basis="bs"), facedes, durations=0, blockids=facedes$run, 
#                       blocklens=rep(436/2,max(facedes$run)), TR=2)
#   
#   term <- construct(mspec$varspec[[1]], mspec)
#   con <- contrast(A=repnum==-1, B=repnum==1)
#   expect_equal(as.vector(contrast_weights(con, term)), c(1,-1,0,0,0))
# })
# 
# test_that("can build a contrast versus the intercept from a convolved term", {
#   facedes$repnum <- factor(facedes$rep_num)
#   aux_table <- data.frame(run=rep(1:6, each=218))
#   mspec <- fmri_model(onset ~  hrf(repnum) + baseline(degree=3, basis="bs"), facedes, durations=0, blockids=facedes$run, 
#                       blocklens=rep(436/2,max(facedes$run)), TR=2)
#   
#   term <- construct(mspec$varspec[[1]], mspec)
#   con <- contrast(A=repnum==-1)
#   expect_equal(as.vector(contrast_weights(con, term)), c(1,0,0,0,0))
# })
# 
# test_that("can build a contrast versus the intercept and add to hrfspec", {
#   facedes$repnum <- factor(facedes$rep_num)
#   aux_table <- data.frame(run=rep(1:6, each=218))
#   
#   
#   conf <- contrast_formula(~ `2` - !`1`, id="repnum")
#   sframe <- sampling_frame(rep(436/2,max(facedes$run)), TR=2)
#   nuisance <- matrix(rnorm(2*length(sframe$blockids)), length(sframe$blockids), 2)
#   
#   bm <- baseline_model(basis="bs", degree=3, sampling_frame=sframe, nuisance_matrix=nuisance)
#   em <- event_model(onset ~ hrf(repnum, contrasts=con), block = ~ run, data=facedes, sampling_frame=sframe)
#   mod <- fmri_model(em, bm)
#   
#   term <- construct(mspec$varspec[[1]], mspec)
#   expect_equal(as.vector(contrast_weights(con, term)), c(1,0,0,0,0))
# })
# 
# test_that("can build a linear contrast from repnum and value_map", {
#   facedes$repnum <- factor(facedes$rep_num)
#   aux_table <- data.frame(run=rep(1:6, each=218))
#   con <- poly_contrast(A=repnum, value_map=list("-1"=0, "1"=1, "2"=2, "3"=3, "4"=4))
#   
#   sframe <- sampling_frame(rep(436/2,max(facedes$run)), TR=2)
#   mspec <- fmri_model(onset ~  hrf(repnum, contrasts=con), ~ baseline(degree=3, basis="bs"), 
#                       event_table=facedes, event_block_ids=facedes$run, sampling_frame=sframe)
#   
#   term <- construct(mspec$varspec[[1]], mspec)
#   expect_equal(as.vector(contrast_weights(con, term)), as.vector(poly(c(0,1,2,3,4))))
# })

