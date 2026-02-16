# test-convolution-fixes.R
#
# Comprehensive tests for the convolution sampling and durations fixes

test_that("convolution preserves matrix dimensions (no 3x expansion)", {
  # Test case that previously caused 3x dimension expansion
  test_events <- data.frame(
    Onset = rep(c(2.0, 6.0, 10.0), 3),
    Run = rep(1:3, each = 3),
    Condition = rep(c("A", "B", "C"), 3),
    stringsAsFactors = FALSE
  )
  
  TRs_per_run <- 10
  n_runs <- 3
  expected_rows <- TRs_per_run * n_runs
  
  sframe <- fmrihrf::sampling_frame(rep(TRs_per_run, n_runs), TR = 2.0)
  emodel <- event_model(Onset ~ hrf(Condition), 
                        block = ~ Run, 
                        sampling_frame = sframe, 
                        data = test_events)
  
  dmat <- design_matrix(emodel)
  
  # Critical test: dimensions should match total TRs, not 3x
  expect_equal(nrow(dmat), expected_rows,
               info = "Design matrix should have correct dimensions (not 3x expanded)")
  
  # Should have 3 columns (one per condition)
  expect_equal(ncol(dmat), 3,
               info = "Design matrix should have one column per condition")
})

test_that("multi-run events are properly convolved (no zero columns)", {
  # Test case with events distributed across runs
  test_events <- data.frame(
    Onset = c(2.0, 6.0, 10.0,   # Run 1
              2.0, 6.0, 10.0,   # Run 2  
              2.0, 6.0, 10.0),  # Run 3
    Run = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
    CondA = rep(c("old", "new", "lure"), each = 3),
    CondB = rep(c("X", "Y", "Z"), 3),
    stringsAsFactors = FALSE
  )
  
  sframe <- sampling_frame(rep(20, 3), TR = 1.5)
  emodel <- event_model(Onset ~ hrf(CondA, CondB, durations = 1.0), 
                        block = ~ Run, 
                        sampling_frame = sframe, 
                        data = test_events)
  
  dmat <- design_matrix(emodel)
  col_sums <- apply(dmat, 2, sum)
  
  # No column should have zero sum
  expect_true(all(col_sums > 0),
              info = paste("All columns should have non-zero sums. Found zeros in:",
                          paste(names(col_sums)[col_sums == 0], collapse = ", ")))
  
  # Check specific combinations exist
  expect_true("CondA_CondB_CondA.old_CondB.X" %in% colnames(dmat),
              info = "Should have column for old:X combination")
  expect_true("CondA_CondB_CondA.new_CondB.Y" %in% colnames(dmat),
              info = "Should have column for new:Y combination")
})

test_that("durations from hrf() are properly propagated", {
  test_events <- data.frame(
    Onset = c(1.0, 5.0, 9.0),
    Run = rep(1, 3),
    Cond = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )
  
  sframe <- sampling_frame(15, TR = 2.0)
  
  # Test 1: durations in hrf()
  emodel1 <- event_model(Onset ~ hrf(Cond, durations = 3.5), 
                         block = ~ Run, 
                         sampling_frame = sframe, 
                         data = test_events)
  
  event_term1 <- emodel1$terms[[1]]
  expect_equal(unique(event_term1$durations), 3.5,
               info = "Durations specified in hrf() should be propagated")
  expect_equal(length(event_term1$durations), nrow(test_events),
               info = "Should have one duration per event")
  
  # Test 2: durations at model level
  emodel2 <- event_model(Onset ~ hrf(Cond), 
                         block = ~ Run, 
                         sampling_frame = sframe, 
                         data = test_events,
                         durations = 2.5)
  
  event_term2 <- emodel2$terms[[1]]
  expect_equal(unique(event_term2$durations), 2.5,
               info = "Model-level durations should be used when hrf() doesn't specify")
  
  # Test 3: hrf() durations override model-level
  emodel3 <- event_model(Onset ~ hrf(Cond, durations = 1.5), 
                         block = ~ Run, 
                         sampling_frame = sframe, 
                         data = test_events,
                         durations = 2.5)
  
  event_term3 <- emodel3$terms[[1]]
  expect_equal(unique(event_term3$durations), 1.5,
               info = "hrf() durations should override model-level durations")
})

test_that("balanced factorial designs remain balanced after convolution", {
  # Create a perfectly balanced 3x3 factorial design
  conditions <- expand.grid(
    Factor1 = c("A", "B", "C"),
    Factor2 = c("X", "Y", "Z")
  )

  # Create 54 events (18 per run, 2 reps of 9 conditions per run)
  # with strictly increasing onsets within each run
  test_events <- data.frame(
    Onset = c(seq(2, 36, by = 2),   # Run 1: 18 events at 2, 4, ..., 36
              seq(2, 36, by = 2),    # Run 2: 18 events at 2, 4, ..., 36
              seq(2, 36, by = 2)),   # Run 3: 18 events at 2, 4, ..., 36
    Run = rep(1:3, each = 18),
    Factor1 = rep(rep(conditions$Factor1, 2), 3),
    Factor2 = rep(rep(conditions$Factor2, 2), 3),
    stringsAsFactors = FALSE
  )

  sframe <- sampling_frame(rep(100, 3), TR = 1.0)
  emodel <- event_model(Onset ~ hrf(Factor1, Factor2, durations = 1.0),
                        block = ~ Run,
                        sampling_frame = sframe,
                        data = test_events)

  dmat <- design_matrix(emodel)
  col_sums <- apply(dmat, 2, sum)

  # All 9 combinations should be present
  expect_equal(ncol(dmat), 9,
               label = "Should have 9 columns for 3x3 factorial")

  # All columns should have similar sums (balanced design)
  sum_range <- range(col_sums)
  relative_range <- (sum_range[2] - sum_range[1]) / mean(col_sums)
  expect_lt(relative_range, 0.1,
            label = "Column sums should be approximately equal in balanced design")

  # No zero columns
  expect_true(all(col_sums > 0),
              info = "No columns should have zero sum in balanced design")
})

test_that("events outside sampling window are handled correctly", {
  test_events <- data.frame(
    Onset = c(2.0, 10.0, 25.0),  # Last event is outside 20s window
    Run = rep(1, 3),
    Cond = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )
  
  # Create a short sampling frame (10 TRs * 2s = 20s)
  sframe <- sampling_frame(10, TR = 2.0)
  
  expect_warning(
    emodel <- event_model(Onset ~ hrf(Cond), 
                          block = ~ Run, 
                          sampling_frame = sframe, 
                          data = test_events),
    regexp = NA,  # No specific warning expected, but allow for any
    info = "Model creation should handle out-of-window events"
  )
  
  dmat <- design_matrix(emodel)
  col_sums <- apply(dmat, 2, sum)

  # Condition C (onset at 25s) should have minimal contribution
  # since it's outside the 20s window
  cond_C_col <- grep("Cond.C", colnames(dmat), value = TRUE)
  cond_A_col <- grep("Cond.A", colnames(dmat), value = TRUE)
  if (length(cond_C_col) == 1 && length(cond_A_col) == 1) {
    # Out-of-window event should contribute less than in-window events
    expect_lt(col_sums[cond_C_col], col_sums[cond_A_col])
  }
})

test_that("block-specific sampling preserves global event timing", {
  # Events with same relative timing in different runs
  test_events <- data.frame(
    Onset = rep(c(1.0, 3.0, 5.0), 3),
    Run = rep(1:3, each = 3),
    Stimulus = rep(c("face", "house", "face"), 3),
    stringsAsFactors = FALSE
  )
  
  sframe <- sampling_frame(rep(10, 3), TR = 1.0)
  emodel <- event_model(Onset ~ hrf(Stimulus, durations = 0.5), 
                        block = ~ Run, 
                        sampling_frame = sframe, 
                        data = test_events)
  
  dmat <- design_matrix(emodel)
  
  # Extract face and house regressors
  face_col <- grep("face", colnames(dmat), value = TRUE)
  house_col <- grep("house", colnames(dmat), value = TRUE)
  
  expect_length(face_col, 1)
  expect_length(house_col, 1)
  
  # Face events (2 per run) should have higher sum than house (1 per run)
  # More frequent events should have larger column sums
  expect_gt(sum(dmat[[face_col]]), sum(dmat[[house_col]]))

  # Check that responses appear in correct blocks
  block1_rows <- 1:10
  block2_rows <- 11:20
  block3_rows <- 21:30

  # Each block should have signal for both conditions
  expect_gt(sum(dmat[block1_rows, face_col]), 0)  # Block 1 should have face signal
  expect_gt(sum(dmat[block2_rows, face_col]), 0)  # Block 2 should have face signal
  expect_gt(sum(dmat[block3_rows, face_col]), 0)  # Block 3 should have face signal
})

test_that("regression test: original problematic case", {
  # Synthesize data matching the original events_testdat.txt structure:
  # 3 runs x 169 TRs @ TR=1.77, Repetition (3 levels) x Saliency (5 levels)
  set.seed(42)
  runs <- 1:3
  rep_levels <- c("rep1", "rep2", "rep3")
  sal_levels <- c("high", "medhigh", "med", "medlow", "low")
  events_list <- lapply(runs, function(r) {
    onsets <- sort(runif(30, min = 2, max = 169 * 1.77 - 10))
    data.frame(
      Onset = onsets,
      Repetition = factor(rep(rep_levels, length.out = 30), levels = rep_levels),
      Saliency = factor(rep(sal_levels, length.out = 30), levels = sal_levels),
      Run = r,
      stringsAsFactors = FALSE
    )
  })
  events <- do.call(rbind, events_list)
  events$Input <- factor(events$Saliency)
  
  sframe <- sampling_frame(rep(169, 3), TR = 1.77)
  
  # This should now work without zero columns (except for out-of-window events)
  emodel <- event_model(Onset ~ hrf(Repetition, Input, durations = 3.5), 
                        block = ~ Run, 
                        sampling_frame = sframe, 
                        data = events)
  
  dmat <- design_matrix(emodel)
  
  # Check dimensions
  expect_equal(nrow(dmat), 507,
               info = "Should have 507 rows (169 * 3)")
  expect_equal(ncol(dmat), 15,
               info = "Should have 15 columns (3 x 5 factorial)")
  
  # Most columns should be non-zero
  col_sums <- apply(dmat, 2, sum)
  zero_cols <- sum(col_sums == 0)
  expect_lt(zero_cols / ncol(dmat), 0.2,
            label = "fraction of zero columns")
})