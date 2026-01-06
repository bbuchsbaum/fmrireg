# Test contrast generation functions: generate_main_effect_contrast, generate_interaction_contrast

library(testthat)
library(fmrireg)

# ============================================================================
# Test generate_main_effect_contrast
# ============================================================================

test_that("generate_main_effect_contrast works for single factor", {
  # Simple 3-level factor

  des <- data.frame(Time = factor(1:3))

  M <- generate_main_effect_contrast(des, "Time")

  # Should have 3 rows (cells) and 2 columns (L-1 contrasts)
  expect_equal(nrow(M), 3)
  expect_equal(ncol(M), 2)

  # Each row sums should indicate difference coding
  # Difference coding matrix has specific structure
  expect_true(is.matrix(M))
})

test_that("generate_main_effect_contrast works for 2-level factor", {
  des <- data.frame(Cond = factor(c("A", "B")))

  M <- generate_main_effect_contrast(des, "Cond")

  # Should have 2 rows and 1 column (2-1=1 contrast)
  expect_equal(nrow(M), 2)
  expect_equal(ncol(M), 1)
})

test_that("generate_main_effect_contrast works for 4-level factor", {
  des <- data.frame(Region = factor(c("north", "south", "east", "west")))

  M <- generate_main_effect_contrast(des, "Region")

  expect_equal(nrow(M), 4)
  expect_equal(ncol(M), 3)  # 4-1 = 3 contrasts
})

test_that("generate_main_effect_contrast errors with multiple factors", {
  des <- data.frame(
    A = factor(1:2),
    B = factor(1:3)
  )

  expect_error(
    generate_main_effect_contrast(des, c("A", "B")),
    regexp = "exactly one factor"
  )
})

test_that("generate_main_effect_contrast errors with invalid factor name", {
  des <- data.frame(A = factor(1:3))

  expect_error(
    generate_main_effect_contrast(des, "nonexistent"),
    regexp = "all.*factors.*names"
  )
})

# ============================================================================
# Test generate_interaction_contrast
# ============================================================================

test_that("generate_interaction_contrast works for 2x2 factorial", {
  des <- expand.grid(
    Time = factor(1:2),
    Cond = factor(c("A", "B"))
  )

  I <- generate_interaction_contrast(des, c("Time", "Cond"))

  # For 2x2: (2-1)*(2-1) = 1 contrast column
  expect_equal(nrow(I), 4)  # 4 cells
  expect_equal(ncol(I), 1)  # 1 interaction contrast
})

test_that("generate_interaction_contrast works for 3x2 factorial", {
  des <- expand.grid(
    Time = factor(1:3),
    Cond = factor(c("face", "scene"))
  )

  I <- generate_interaction_contrast(des, c("Time", "Cond"))

  # For 3x2: (3-1)*(2-1) = 2 contrasts
  expect_equal(nrow(I), 6)  # 6 cells
  expect_equal(ncol(I), 2)  # 2 interaction contrasts
})

test_that("generate_interaction_contrast works for 4x2 factorial (example from docs)", {
  des <- expand.grid(
    Time = factor(1:4),
    Cond = factor(c("face", "scene"))
  )

  I <- generate_interaction_contrast(des, c("Time", "Cond"))

  # For 4x2: (4-1)*(2-1) = 3 contrasts
  expect_equal(nrow(I), 8)  # 8 cells
  expect_equal(ncol(I), 3)  # 3 interaction contrasts
})

test_that("generate_interaction_contrast works for 3x3 factorial", {
  des <- expand.grid(
    FactorA = factor(1:3),
    FactorB = factor(c("X", "Y", "Z"))
  )

  I <- generate_interaction_contrast(des, c("FactorA", "FactorB"))

  # For 3x3: (3-1)*(3-1) = 4 contrasts
  expect_equal(nrow(I), 9)   # 9 cells
  expect_equal(ncol(I), 4)   # 4 interaction contrasts
})

test_that("generate_interaction_contrast with single factor gives main effect", {
  des <- expand.grid(Time = factor(1:4), Cond = factor(c("A", "B")))

  # Single factor should reproduce main effect
  I_single <- generate_interaction_contrast(des, "Time")
  M_main <- generate_main_effect_contrast(
    data.frame(Time = factor(1:4)),
    "Time"
  )

  # The contrast structure for Time should have same number of contrast columns
  expect_equal(ncol(I_single), ncol(M_main))  # Both = 3
})

test_that("generate_interaction_contrast errors with invalid factor names", {
  des <- expand.grid(A = factor(1:2), B = factor(1:2))

  expect_error(
    generate_interaction_contrast(des, c("A", "C")),
    regexp = "all.*factors.*names"
  )
})

test_that("generate_interaction_contrast works with 3-way interaction", {
  des <- expand.grid(
    A = factor(1:2),
    B = factor(1:2),
    C = factor(1:2)
  )

  I <- generate_interaction_contrast(des, c("A", "B", "C"))

  # For 2x2x2: (2-1)*(2-1)*(2-1) = 1 three-way interaction contrast
  expect_equal(nrow(I), 8)   # 8 cells
  expect_equal(ncol(I), 1)   # 1 three-way interaction contrast
})

test_that("generate_interaction_contrast works for 3-way with different levels", {
  des <- expand.grid(
    A = factor(1:3),
    B = factor(1:2),
    C = factor(1:4)
  )

  I <- generate_interaction_contrast(des, c("A", "B", "C"))

  # For 3x2x4: (3-1)*(2-1)*(4-1) = 6 contrasts
  expect_equal(nrow(I), 24)  # 24 cells
  expect_equal(ncol(I), 6)   # 6 contrasts
})

# ============================================================================
# Test contrast orthogonality and properties
# ============================================================================

test_that("main effect contrasts are orthogonal to grand mean", {
  des <- data.frame(Cond = factor(1:4))
  M <- generate_main_effect_contrast(des, "Cond")

  # Each contrast should sum to zero (orthogonal to grand mean)
  col_sums <- colSums(M)
  expect_true(all(abs(col_sums) < 1e-10))
})

test_that("interaction contrasts preserve cell structure", {
  des <- expand.grid(A = factor(1:2), B = factor(1:3))

  I <- generate_interaction_contrast(des, c("A", "B"))

  # Rows should correspond to cells in lexicographic order
  expect_equal(nrow(I), nrow(des))
})

test_that("contrast generation is consistent with repeated calls", {
  des <- expand.grid(Time = factor(1:3), Cond = factor(c("A", "B")))

  I1 <- generate_interaction_contrast(des, c("Time", "Cond"))
  I2 <- generate_interaction_contrast(des, c("Time", "Cond"))

  expect_identical(I1, I2)
})

# ============================================================================
# Edge cases
# ============================================================================

test_that("contrast generation handles single-cell factor in interaction", {
  # One factor has only 1 level - unusual but should work
  des <- expand.grid(
    A = factor(1:3),
    B = factor("only_one")
  )

  # For 3x1: (3-1)*(1-1) = 0 interaction contrasts
  I <- generate_interaction_contrast(des, c("A", "B"))
  expect_equal(ncol(I), 0)  # No interaction possible
})

test_that("generate_main_effect_contrast handles unbalanced factor", {
  # Unbalanced design (more entries for some levels)
  des <- data.frame(
    Cond = factor(c("A", "A", "B", "A", "B", "C"))
  )

  # Should still work - uses unique levels
  M <- generate_main_effect_contrast(des, "Cond")

  # 3 levels -> 2 contrasts
  expect_equal(ncol(M), 2)
})

test_that("contrast matrices can be used for hypothesis testing", {
  # Verify the contrast matrix can be applied to coefficient vector
  des <- expand.grid(Time = factor(1:3), Cond = factor(c("A", "B")))

  M <- generate_main_effect_contrast(des, "Time")
  I <- generate_interaction_contrast(des, c("Time", "Cond"))

  # Simulate coefficients for each cell
  set.seed(123)
  betas <- rnorm(nrow(des))

  # Apply contrasts
  main_effects <- crossprod(M, betas)
  interaction_effects <- crossprod(I, betas)

  expect_equal(length(main_effects), ncol(M))
  expect_equal(length(interaction_effects), ncol(I))
})

# ============================================================================
# Test with real fMRI design patterns
# ============================================================================

test_that("contrast generation works with typical fMRI factorial design", {
  # Typical 2 (stimulus type) x 3 (repetition) design
  des <- expand.grid(
    stimulus = factor(c("face", "scene")),
    repetition = factor(1:3)
  )

  # Main effects
  M_stim <- generate_main_effect_contrast(
    data.frame(stimulus = factor(c("face", "scene"))),
    "stimulus"
  )
  M_rep <- generate_main_effect_contrast(
    data.frame(repetition = factor(1:3)),
    "repetition"
  )

  # Interaction
  I <- generate_interaction_contrast(des, c("stimulus", "repetition"))

  # Verify dimensions

  expect_equal(ncol(M_stim), 1)  # 2-1 = 1

  expect_equal(ncol(M_rep), 2)  # 3-1 = 2
  expect_equal(ncol(I), 2)      # (2-1)*(3-1) = 2
  expect_equal(nrow(I), 6)      # 2*3 = 6 cells
})
