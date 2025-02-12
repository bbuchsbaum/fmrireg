options(mc.cores=2)

library(testthat)
library(assertthat)
facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)


test_that("a 2-by-2 Fcontrast", {
  F_a <- factor(rep(letters[1:2], 8))
  F_b <- factor(rep(c("V1", "V2"), each=8))
  onsets <- seq(1, length(F_a))
  blockids <- rep(1, length(onsets))
  
  et <- event_term(list(Fa=F_a, Fb=F_b), onsets, blockids)
  expect_true(!is.null(Fcontrasts(et)))
})

test_that("a 2-by-3 Fcontrast", {
  F_a <- factor(rep(letters[1:2], 9))
  F_b <- factor(rep(c("V1", "V2", "V3"), each=6))
 
  onsets <- seq(1, length(F_a))
  blockids <- rep(1, length(onsets))
  
  et <- event_term(list(Fa=F_a, Fb=F_b), onsets, blockids)
  expect_true(!is.null(Fcontrasts(et)))
})


test_that("a 3-by-2 Fcontrast", {
  F_a <- factor(rep(letters[1:2], 8))
  F_b <- factor(rep(c("V1", "V2"), each=8))
  F_c <- factor(rep(c("B3", "B3", "B4", "B4"), 4))
  onsets <- seq(1, length(F_a))
  blockids <- rep(1, length(onsets))
  
  et <- event_term(list(Fa=F_a, Fb=F_b, Fc=F_c), onsets, blockids)
  expect_true(!is.null(Fcontrasts(et)))
})

test_that("a 3-by-3 Fcontrast", {
  F_a <- factor(sample(letters[1:3], 200, replace=TRUE))
  F_b <- factor(sample(c("V1", "V2", "V3"), 200, replace=TRUE))
  F_c <- factor(sample(c("B3", "B3", "B4", "B4", "B5", "B5"),200, replace=TRUE))
  onsets <- seq(1, length(F_a))
  blockids <- rep(1, length(onsets))
  
  et <- event_term(list(Fa=F_a, Fb=F_b, Fc=F_c), onsets, blockids)
  expect_true(!is.null(Fcontrasts(et)))
})




test_that("can build a simple contrast from a convolved term", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  con <- pair_contrast(~ repnum==-1, ~ repnum==1, name="A_B")
  
  expect_equal(as.vector(contrast_weights(con, terms(espec)[[1]])$weights), c(1,-1,0,0,0))
})

test_that("can build a simple contrast from a convolved term and convert to glt", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  con <- pair_contrast(~ repnum==-1, ~ repnum==1, name="A_B")
  
  conw <- contrast_weights(con, terms(espec)[[1]])
  glt <- to_glt(conw)
  expect_true(!is.null(glt))
})

test_that("can build a contrast versus the intercept from a convolved term", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  
  con <- unit_contrast(~ repnum==-1, name="A")

  term <- terms(espec)[[1]]
  expect_equal(as.vector(contrast_weights(con, term)$weights), rep(.2,5))
  
})

test_that("can construct a simple pair_contrast", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  
  pc <- pair_contrast(~ repnum == 1, ~ repnum ==2, name="B-A")
  cw <- contrast_weights(pc, terms(espec)[[1]])
  expect_equal(as.vector(cw$weights), c(0,1,-1,0,0))
  
})

test_that("can build a linear contrast from repnum and value_map", {
   facedes$repnum <- factor(facedes$rep_num)
   sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
   espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  
   con <- poly_contrast(~ repnum, degree=1, value_map=list("-1"=0, "1"=1, "2"=2, "3"=3, "4"=4), name="linear_repnum")
   term1 <- terms(espec)[[1]]
   cw <- contrast_weights(con, term1)
   expect_equal(as.vector(contrast_weights(con, term1)$weights), as.vector(poly(c(0,1,2,3,4))))
})

test_that("can build a set of pairwise contrasts", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  levs <- levels(facedes$repnum)
  cset <- pairwise_contrasts(levs)
  expect_equal(length(cset), ncol(combn(length(levs),2)))

})

test_that("can build a one_against_all contrast set", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  levs <- levels(facedes$repnum)
  cset <- one_against_all_contrast(levs, "repnum")
  
  
  expect_equal(length(cset),length(levels(facedes$repnum)))
  
  wtls <- lapply(cset, function(con) {
    contrast_weights(con, terms(espec)[[1]])
  })
  expect_true(!is.null(wtls))
  
})

test_that("can subtract two pairwise contrasts to form an interaction contrast", {
  simple_des <- expand.grid(category=c("face", "scene"), attention=c("attend", "ignored"), replication=c(1,2))
  simple_des$onset <- seq(1,100, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  sframe <- sampling_frame(blocklens=100, TR=2)
  espec <- event_model(onset ~  hrf(category, attention), data=simple_des, block=~run, sampling_frame=sframe)
  con1 <- pair_contrast(~ category=="face", ~ category == "scene", name="face_scene#attend", where=~ attention == "attend")
  con2 <- pair_contrast(~ category=="face", ~ category == "scene", name="face_scene#ignored", where=~ attention == "ignored")
  con3 <- con1 - con2
  expect_true(!is.null(con3))
})


test_that("can contrast two parametric regressors crossed with a factor", {
  simple_des <- expand.grid(category=c("face", "scene"), attention=c("attend", "ignored"), replication=c(1,2))
  simple_des$onset <- seq(1,100, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  simple_des$RT <- rnorm(nrow(simple_des))
  
  sframe <- sampling_frame(blocklens=100, TR=2)
  espec <- event_model(onset ~  hrf(category, RT), data=simple_des, block=~run, sampling_frame=sframe)
  con <- pair_contrast(~ category == "face", ~ category == "scene", name="face_scene_RT")
  cwts <- contrast_weights(con, terms(espec)[[1]])
  expect_true(!is.null(cwts))
})

test_that("can contrast two parametric regressors wrapped in Ident for additive regressors", {
  simple_des <- expand.grid(category=c("face", "scene"), attention=c("attend", "ignored"), replication=c(1,2))
  simple_des$onset <- seq(1,100, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  simple_des$RT1 <- rnorm(nrow(simple_des))
  simple_des$RT2 <- rnorm(nrow(simple_des))
  
  sframe <- sampling_frame(blocklens=100, TR=2)
  espec <- event_model(onset ~  hrf(Ident(RT1,RT2)), data=simple_des, block=~run, sampling_frame=sframe)
  con <- pair_contrast(~ RT1_RT2 == "RT1", ~ RT1_RT2 == "RT2", name="RT1_RT2")
  cwts <- contrast_weights(con, terms(espec)[[1]])
  expect_true(!is.null(cwts))
})

test_that("can contrast two basis functions from a custom multi-phase hrf", {
  simple_des <- expand.grid(trial_type=rep(c("encoding"),15))
  
  simple_des$onset <- seq(1,300, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  
  hrf_encode <- gen_hrf(hrf_spmg1, normalize=TRUE)
  hrf_delay <- gen_hrf(hrf_spmg1, lag=3, width=8, normalize=TRUE)
  hrf_probe <-gen_hrf(hrf_spmg1, lag=11, width=3, normalize=TRUE)  
  hrf_trial <<- gen_hrf_set(hrf_encode, hrf_delay, hrf_probe)
  
  sframe <- sampling_frame(blocklens=250, TR=2)
  espec <- event_model(onset ~  hrf(trial_type, basis=hrf_trial), data=simple_des, block=~run, sampling_frame=sframe)
  con <- pair_contrast(~ basis == "basis1", ~ basis=="basis2", name="test")
  cwts <- contrast_weights(con, terms(espec)[[1]])
  expect_true(!is.null(cwts))
})

test_that("can contrast two parametric regressors wrapped in Ident", {
  simple_des <- expand.grid(category=c("face", "scene"), attention=c("attend", "ignored"), replication=c(1,2))
  simple_des$onset <- seq(1,100, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  simple_des$RT1 <- rnorm(nrow(simple_des))
  simple_des$RT2 <- rnorm(nrow(simple_des))
  
  sframe <- sampling_frame(blocklens=100, TR=2)
  espec <- event_model(onset ~  hrf(Ident(RT1,RT2)), data=simple_des, block=~run, sampling_frame=sframe)
  con <- pair_contrast(~ RT1_RT2 == "RT1", ~ RT1_RT2 == "RT2", name="RT1_RT2")
  cwts <- contrast_weights(con, terms(espec)[[1]])
  expect_true(!is.null(cwts))
})


test_that("can form a simple formula contrast", {
  simple_des <- expand.grid(category=c("face", "scene"), attention=c("attend", "ignored"), replication=c(1,2))
  simple_des$onset <- seq(1,100, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  sframe <- sampling_frame(blocklens=100, TR=2)
  espec <- event_model(onset ~  hrf(category, attention), data=simple_des, block=~run, sampling_frame=sframe)
  con1 <- contrast(~ (face:attend - face:ignored) - (scene:attend - scene:ignored), name="face_scene")
  cwts <- contrast_weights(con1, terms(espec)[[1]])
  expect_equal(as.vector(cwts$weights[,1]), c(1, -1, -1, 1))
  
})

test_that("can form formula contrast with 3 terms", {
  simple_des <- expand.grid(match=c("match", "nonmatch"), condition=c("NOVEL", "REPEAT"), correct=c("correct","incorrect"))
  simple_des$onset <- seq(1,100, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  sframe <- sampling_frame(blocklens=100, TR=2)
  espec <- event_model(onset ~  hrf(match, condition, correct), data=simple_des, block=~run, sampling_frame=sframe)
  #con1 <- contrast(~ (face:ATT:r1 + face:IG:r2) - (scene:ATT:r1 - scene:ignored:r2), name="face_scene")
  con1 <- contrast(
    ~  ((match:NOVEL:correct + match:NOVEL:incorrect) - (nonmatch:NOVEL:correct + nonmatch:NOVEL:incorrect)) -
      ((match:REPEAT:correct + match:REPEAT:incorrect) - (nonmatch:REPEAT:correct + nonmatch:REPEAT:incorrect)), name="cond_by_match")
  cwts <- contrast_weights(con1, terms(espec)[[1]])
  expect_equal(length(as.vector(cwts$weights[,1])), 8)
  
})

test_that("can form formula contrast with two factor terms and one continuous covariate", {
  simple_des <- expand.grid(match=c("match", "nonmatch"), condition=c("NOVEL", "REPEAT"), correct=c(1,2))
  simple_des$onset <- seq(1,100, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  simple_des$correct <- as.factor(simple_des$correct)
  sframe <- sampling_frame(blocklens=100, TR=2)
  espec <- event_model(onset ~  hrf(match, condition, correct), data=simple_des, block=~run, sampling_frame=sframe)
  #con1 <- contrast(~ (face:ATT:r1 + face:IG:r2) - (scene:ATT:r1 - scene:ignored:r2), name="face_scene")
  con1 <- contrast(
    ~  match:NOVEL:1 - nonmatch:NOVEL:1, name="cond_by_match")
  cwts <- contrast_weights(con1, terms(espec)[[1]])
  expect_equal(as.vector(cwts$weights[,1]), c(1, -1, 0,0,0,0,0,0))
  
})





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

