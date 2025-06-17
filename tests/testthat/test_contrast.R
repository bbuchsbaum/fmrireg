options(mc.cores=2)

library(testthat)
library(assertthat)
library(fmrireg)
library(fmrihrf)
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
  sframe <- fmrihrf::sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  con <- pair_contrast(~ repnum==-1, ~ repnum==1, name="A_B")
  
  expect_equal(as.vector(contrast_weights(con, terms(espec)[[1]])$weights), c(1,-1,0,0,0))
})

test_that("can build a simple contrast from a convolved term and convert to glt", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- fmrihrf::sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  con <- pair_contrast(~ repnum==-1, ~ repnum==1, name="A_B")
  
  conw <- contrast_weights(con, terms(espec)[[1]])
  glt <- to_glt(conw)
  expect_true(!is.null(glt))
})

test_that("can build a contrast versus the intercept from a convolved term", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- fmrihrf::sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  
  con <- unit_contrast(~ repnum, name="A")

  term <- terms(espec)[[1]]
  expect_equal(as.vector(contrast_weights(con, term)$weights), rep(.2,5))
  
})

test_that("can construct a simple pair_contrast", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- fmrihrf::sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  
  pc <- pair_contrast(~ repnum == 1, ~ repnum ==2, name="B-A")
  cw <- contrast_weights(pc, terms(espec)[[1]])
  expect_equal(as.vector(cw$weights), c(0,1,-1,0,0))
  
})

test_that("can build a linear contrast from repnum and value_map", {
   facedes$repnum <- factor(facedes$rep_num)
   sframe <- fmrihrf::sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
   espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  
   con <- poly_contrast(~ repnum, degree=1, value_map=list("-1"=0, "1"=1, "2"=2, "3"=3, "4"=4), name="linear_repnum")
   term1 <- terms(espec)[[1]]
   cw <- contrast_weights(con, term1)
   expect_equal(as.vector(contrast_weights(con, term1)$weights), as.vector(poly(c(0,1,2,3,4))))
})

test_that("can build a set of pairwise contrasts", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- fmrihrf::sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
  espec <- event_model(onset ~  hrf(repnum), data=facedes, block=~run, sampling_frame=sframe)
  levs <- levels(facedes$repnum)
  cset <- pairwise_contrasts(levs, facname = "repnum")
  expect_equal(length(cset), ncol(combn(length(levs),2)))

})

test_that("can build a one_against_all contrast set", {
  facedes$repnum <- factor(facedes$rep_num)
  sframe <- fmrihrf::sampling_frame(blocklens=rep(436/2,max(facedes$run)), TR=2)
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
  sframe <- fmrihrf::sampling_frame(blocklens=100, TR=2)
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
  
  sframe <- fmrihrf::sampling_frame(blocklens=100, TR=2)
  espec <- event_model(onset ~  hrf(category), data=simple_des, block=~run, sampling_frame=sframe)
  con <- pair_contrast(~ category == "face", ~ category == "scene", name="face_vs_scene")
  cwts <- contrast_weights(con, terms(espec)[[1]])
  expect_true(!is.null(cwts))
})

test_that("can contrast two parametric regressors wrapped in Ident for additive regressors", {
  simple_des <- expand.grid(category=c("face", "scene"), attention=c("attend", "ignored"), replication=c(1,2))
  simple_des$onset <- seq(1,100, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  simple_des$RT1 <- rnorm(nrow(simple_des))
  simple_des$RT2 <- rnorm(nrow(simple_des))
  
  sframe <- fmrihrf::sampling_frame(blocklens=100, TR=2)
  espec <- event_model(onset ~  hrf(Ident(RT1,RT2)), data=simple_des, block=~run, sampling_frame=sframe)
  
  # --- Determine actual column names for safety ---
  dm_colnames <- colnames(design_matrix(espec))
  term_of_interest <- terms(espec)[[1]] 
  term_condition_tags <- conditions(term_of_interest, expand_basis=FALSE) 
  term_tag_assigned <- attr(term_of_interest, "term_tag")
  
  # For Ident-only terms, term_tag might be NULL according to new naming scheme
  # In that case, column names should be just the variable names (RT1, RT2)
  if (is.null(term_tag_assigned)) {
    # Ident-only case: column names should be just the variable names
    pattern_A_col <- "^RT1$"
    pattern_B_col <- "^RT2$"
  } else {
    # Regular case with term_tag prefix
    if (length(term_condition_tags) < 2) {
      stop("Could not reliably determine column names for Ident(RT1,RT2) term in test.")
    }
    pattern_A_col <- paste0("^", term_tag_assigned, "_", term_condition_tags[1], "$")
    pattern_B_col <- paste0("^", term_tag_assigned, "_", term_condition_tags[2], "$")
  }
  
  # Ensure these patterns actually match columns in the design_matrix
  expect_true(any(grepl(pattern_A_col, dm_colnames)), 
              info = paste("Pattern A:", pattern_A_col, "did not match any of actual colnames:", paste(dm_colnames, collapse=", ")))
  expect_true(any(grepl(pattern_B_col, dm_colnames)),
              info = paste("Pattern B:", pattern_B_col, "did not match any of actual colnames:", paste(dm_colnames, collapse=", ")))

  con <- column_contrast(pattern_A = pattern_A_col, 
                         pattern_B = pattern_B_col, 
                         name="RT1_vs_RT2_cols")
                         
  cwts <- contrast_weights(con, term_of_interest)
  expect_true(!is.null(cwts))
  
  # More specific checks for the weights vector
  weights_vec <- as.vector(cwts$weights)
  col_A_idx <- grep(pattern_A_col, dm_colnames)
  col_B_idx <- grep(pattern_B_col, dm_colnames)
  
  # Check that only one column matches each pattern
  expect_equal(length(col_A_idx), 1, info = "Pattern A should match exactly one column")
  expect_equal(length(col_B_idx), 1, info = "Pattern B should match exactly one column")
  
  if (length(col_A_idx) == 1 && length(col_B_idx) == 1) {
      expect_equal(weights_vec[col_A_idx], 1)
      expect_equal(weights_vec[col_B_idx], -1)
      # Ensure all other weights are zero
      other_indices <- setdiff(seq_along(weights_vec), c(col_A_idx, col_B_idx))
      if (length(other_indices) > 0) {
        expect_true(all(weights_vec[other_indices] == 0))
      }
      expect_equal(sum(weights_vec), 0) # Overall sum should be zero
  }

})

test_that("can contrast two basis functions from a custom multi-phase hrf", {
  # Create a proper factor with multiple levels to avoid the contrasts error
  simple_des <- expand.grid(trial_type=c("encoding", "retrieval"))
  simple_des <- simple_des[rep(1:nrow(simple_des), length.out=15), , drop=FALSE]
  
  simple_des$onset <- seq(1,300, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  
  hrf_encode <- fmrihrf::gen_hrf(fmrihrf::hrf_spmg1, normalize=TRUE)
  hrf_delay <- fmrihrf::gen_hrf(fmrihrf::hrf_spmg1, lag=3, width=8, normalize=TRUE)
  hrf_probe <-fmrihrf::gen_hrf(fmrihrf::hrf_spmg1, lag=11, width=3, normalize=TRUE)  
  hrf_trial <<- fmrihrf::hrf_set(hrf_encode, hrf_delay, hrf_probe)
  
  sframe <- fmrihrf::sampling_frame(blocklens=250, TR=2)
  espec <- event_model(onset ~  hrf(trial_type, basis=hrf_trial), data=simple_des, block=~run, sampling_frame=sframe)
  
  # Updated to use new HRF basis suffix naming scheme: _b01, _b02, etc.
  # Use column_contrast to target specific basis functions
  con <- column_contrast(pattern_A = "_b01$", pattern_B = "_b02$", name="basis1_vs_basis2")
  cwts <- contrast_weights(con, terms(espec)[[1]])
  expect_true(!is.null(cwts))
})

test_that("can form a simple formula contrast", {
  simple_des <- expand.grid(category=c("face", "scene"), attention=c("attend", "ignored"), replication=c(1,2))
  simple_des$onset <- seq(1,100, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  sframe <- fmrihrf::sampling_frame(blocklens=100, TR=2)
  espec <- event_model(onset ~  hrf(category, attention), data=simple_des, block=~run, sampling_frame=sframe)
  
  # Use the new naming scheme: term_tag_condition_tag format
  # The term_tag should be "category_attention" and condition_tags should be "category.face_attention.attend" etc.
  # But for formula contrasts, we still use shortnames which are the old format for backward compatibility
  con1 <- contrast(~ (`face:attend` - `face:ignored`) - (`scene:attend` - `scene:ignored`), name="face_scene")
  cwts <- contrast_weights(con1, terms(espec)[[1]])
  
  # Now the contrast should work correctly with proper weights
  # Order is: face:attend, scene:attend, face:ignored, scene:ignored
  # Formula: (face:attend - face:ignored) - (scene:attend - scene:ignored)
  # Expected: face:attend=1, face:ignored=-1, scene:attend=-1, scene:ignored=1
  expect_equal(as.vector(cwts$weights[,1]), c(1, -1, -1, 1))
  
})

test_that("can form formula contrast with 3 terms", {
  simple_des <- expand.grid(match=c("match", "nonmatch"), condition=c("NOVEL", "REPEAT"), correct=c("correct","incorrect"))
  simple_des$onset <- seq(1,100, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  sframe <- fmrihrf::sampling_frame(blocklens=100, TR=2)
  espec <- event_model(onset ~  hrf(match, condition, correct), data=simple_des, block=~run, sampling_frame=sframe)
  
  # Use the old shortnames format for formula contrasts (backward compatibility)
  con1 <- contrast(
    ~  ((`match:NOVEL:correct` + `match:NOVEL:incorrect`) - (`nonmatch:NOVEL:correct` + `nonmatch:NOVEL:incorrect`)) -
      ((`match:REPEAT:correct` + `match:REPEAT:incorrect`) - (`nonmatch:REPEAT:correct` + `nonmatch:REPEAT:incorrect`)), name="cond_by_match")
  cwts <- contrast_weights(con1, terms(espec)[[1]])
  expect_equal(length(as.vector(cwts$weights[,1])), 8)
  
})

test_that("can form formula contrast with two factor terms and one continuous covariate", {
  simple_des <- expand.grid(match=c("match", "nonmatch"), condition=c("NOVEL", "REPEAT"), correct=c(1,2))
  simple_des$onset <- seq(1,100, length.out=nrow(simple_des))
  simple_des$run <- rep(1,nrow(simple_des))
  simple_des$correct <- as.factor(simple_des$correct)
  sframe <- fmrihrf::sampling_frame(blocklens=100, TR=2)
  espec <- event_model(onset ~  hrf(match, condition, correct), data=simple_des, block=~run, sampling_frame=sframe)
  
  # Use the old shortnames format for formula contrasts
  con1 <- contrast(
    ~  `match:NOVEL:1` - `nonmatch:NOVEL:1`, name="cond_by_match")
  cwts <- contrast_weights(con1, terms(espec)[[1]])
  
  # Check that we get the correct contrast weights
  weights_vec <- as.vector(cwts$weights[,1])
  expect_equal(length(weights_vec), 8)
  
  # The contrast should have +1 for match:NOVEL:1 and -1 for nonmatch:NOVEL:1
  # Find the positions of these conditions in shortnames
  short_names <- shortnames(terms(espec)[[1]])
  match_pos <- which(short_names == "match:NOVEL:1")
  nonmatch_pos <- which(short_names == "nonmatch:NOVEL:1")
  
  expect_equal(weights_vec[match_pos], 1)
  expect_equal(weights_vec[nonmatch_pos], -1)
  # All other weights should be 0
  other_pos <- setdiff(1:length(weights_vec), c(match_pos, nonmatch_pos))
  expect_true(all(weights_vec[other_pos] == 0))
  
})



# 
# test_that("can build a contrast versus the intercept and add to hrfspec", {
#   facedes$repnum <- factor(facedes$rep_num)
#   aux_table <- data.frame(run=rep(1:6, each=218))
#   
#   
#   conf <- contrast_formula(~ `2` - !`1`, id="repnum")
#   sframe <- fmrihrf::sampling_frame(rep(436/2,max(facedes$run)), TR=2)
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


test_that("pair_contrast respects where clause", {
  simple_des <- expand.grid(category=c("face", "scene"), attention=c("attend", "ignored"))
  simple_des$onset <- seq(1, 50, length.out=nrow(simple_des))
  simple_des$run <- 1
  sframe <- fmrihrf::sampling_frame(blocklens=50, TR=2)
  espec <- event_model(onset ~ hrf(category, attention), data=simple_des, block=~run, sampling_frame=sframe)

  con <- pair_contrast(~ category == "face", ~ category == "scene", name="face_vs_scene_attend", where = ~ attention == "attend")
  cw <- contrast_weights(con, terms(espec)[[1]])
  wvec <- as.vector(cw$weights)
  sn <- shortnames(terms(espec)[[1]])
  face_att <- which(sn == "face:attend")
  scene_att <- which(sn == "scene:attend")
  face_ign <- which(sn == "face:ignored")
  scene_ign <- which(sn == "scene:ignored")

  expect_equal(wvec[face_att], 1)
  expect_equal(wvec[scene_att], -1)
  expect_equal(wvec[face_ign], 0)
  expect_equal(wvec[scene_ign], 0)
})

test_that("unit_contrast allows where subsetting", {
  simple_des <- expand.grid(category=c("face", "scene"), attention=c("attend", "ignored"))
  simple_des$onset <- seq(1, 50, length.out=nrow(simple_des))
  simple_des$run <- 1
  sframe <- fmrihrf::sampling_frame(blocklens=50, TR=2)
  espec <- event_model(onset ~ hrf(category, attention), data=simple_des, block=~run, sampling_frame=sframe)

  con <- unit_contrast(~ category, name="face_only_attend", where = ~ attention == "attend" & category == "face")
  cw <- contrast_weights(con, terms(espec)[[1]])
  expect_equal(nrow(cw$weights), 1)
  expect_equal(as.vector(cw$weights)[1], 1)
  expect_true(grepl("face:attend", rownames(cw$weights)))
})

