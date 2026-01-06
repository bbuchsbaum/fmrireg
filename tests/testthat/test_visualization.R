# Test visualization functions: autoplot.Reg, design_plot, correlation_map

library(testthat)
library(fmrireg)
library(ggplot2)

# ============================================================================
# Test autoplot.Reg
# ============================================================================

test_that("autoplot.Reg works with single-basis HRF regressor", {
  # Create a simple regressor with canonical HRF
  hrf_obj <- fmrihrf::HRF_SPMG1
  reg <- fmrihrf::regressor(onsets = c(5, 15, 25), hrf = hrf_obj)

  # Should produce a ggplot object without error
  p <- autoplot(reg)
  expect_s3_class(p, "ggplot")

  # Check that the plot has the expected layers
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomLine"), logical(1))))
})

test_that("autoplot.Reg works with multi-basis HRF regressor", {
  # Create a regressor with multi-basis HRF (e.g., FIR)
  hrf_fir <- fmrihrf::hrf_fir(nbasis = 5, span = 20)
  reg <- fmrihrf::regressor(onsets = c(5, 15, 25), hrf = hrf_fir)

  # Should produce a ggplot object without error
  p <- autoplot(reg)
  expect_s3_class(p, "ggplot")

  # Multi-basis should have faceting
  expect_true(!is.null(p$facet))
})

test_that("autoplot.Reg works with custom grid", {
  hrf_obj <- fmrihrf::HRF_SPMG1
  reg <- fmrihrf::regressor(onsets = c(5, 15), hrf = hrf_obj)

  # Provide a custom evaluation grid
  custom_grid <- seq(0, 40, by = 0.5)
  p <- autoplot(reg, grid = custom_grid)

  expect_s3_class(p, "ggplot")
})

test_that("autoplot.Reg works with different precision", {
  hrf_obj <- fmrihrf::HRF_SPMG1
  reg <- fmrihrf::regressor(onsets = c(10), hrf = hrf_obj)

  # Test different precision levels
  p_fine <- autoplot(reg, precision = 0.05)
  p_coarse <- autoplot(reg, precision = 0.5)

  expect_s3_class(p_fine, "ggplot")
  expect_s3_class(p_coarse, "ggplot")
})

test_that("autoplot.Reg handles regressor with no onsets gracefully", {
  # Edge case: empty onsets
  hrf_obj <- fmrihrf::HRF_SPMG1
  reg <- fmrihrf::regressor(onsets = numeric(0), hrf = hrf_obj)

  # Should still work, using default grid
  p <- autoplot(reg)
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.Reg works with various HRF types", {
  onsets <- c(5, 20)

  # Test with different HRF implementations
  hrf_types <- list(
    spmg1 = fmrihrf::HRF_SPMG1,
    gamma = fmrihrf::hrf_gamma(shape = 6, rate = 1),
    gaussian = fmrihrf::hrf_gaussian(mean = 5, sd = 2)
  )

  for (name in names(hrf_types)) {
    reg <- fmrihrf::regressor(onsets = onsets, hrf = hrf_types[[name]])
    p <- autoplot(reg)
    expect_s3_class(p, "ggplot", info = paste("Failed for HRF type:", name))
  }
})

# ============================================================================
# Test .correlation_map_common (internal helper)
# ============================================================================

test_that(".correlation_map_common works with basic matrix input", {
  # Create a simple design matrix
  DM <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  colnames(DM) <- paste0("Reg_", 1:5)

  # Test Pearson correlation
  p <- fmrireg:::.correlation_map_common(DM, method = "pearson")
  expect_s3_class(p, "ggplot")

  # Test Spearman correlation
  p_spearman <- fmrireg:::.correlation_map_common(DM, method = "spearman")
  expect_s3_class(p_spearman, "ggplot")
})

test_that(".correlation_map_common works with half_matrix option", {
  DM <- matrix(rnorm(50 * 4), nrow = 50, ncol = 4)
  colnames(DM) <- paste0("Col_", 1:4)

  # Full matrix
  p_full <- fmrireg:::.correlation_map_common(DM, half_matrix = FALSE)
  expect_s3_class(p_full, "ggplot")

  # Half matrix (lower triangle)
  p_half <- fmrireg:::.correlation_map_common(DM, half_matrix = TRUE)
  expect_s3_class(p_half, "ggplot")
})

test_that(".correlation_map_common handles absolute_limits option", {
  DM <- matrix(rnorm(30 * 3), nrow = 30, ncol = 3)
  colnames(DM) <- c("A", "B", "C")

  # Force scale to -1..+1
  p_abs <- fmrireg:::.correlation_map_common(DM, absolute_limits = TRUE)
  expect_s3_class(p_abs, "ggplot")

  # Data-driven scale
  p_data <- fmrireg:::.correlation_map_common(DM, absolute_limits = FALSE)
  expect_s3_class(p_data, "ggplot")
})

test_that(".correlation_map_common requires at least 2 columns", {
  DM_single <- matrix(rnorm(20), nrow = 20, ncol = 1)
  expect_error(
    fmrireg:::.correlation_map_common(DM_single),
    regexp = "must have at least 2 columns"
  )
})

test_that(".correlation_map_common handles constant columns gracefully", {
  # One column is constant (will produce NA correlations)
  DM <- cbind(rnorm(50), rep(1, 50), rnorm(50))
  colnames(DM) <- c("Var", "Const", "Var2")

  # Should not error (warning is suppressed internally)
  expect_no_error(
    p <- fmrireg:::.correlation_map_common(DM)
  )
  expect_s3_class(p, "ggplot")
})

# ============================================================================
# Test design_plot (Shiny app - limited testing)
# ============================================================================

test_that("design_plot validates input types", {
  # Should error with non-fmri_model input
  expect_error(
    design_plot("not_a_model"),
    regexp = "fmri_model"
  )
})

test_that("design_plot validates term_name parameter", {
  skip_on_cran()
  skip_if_not_installed("shiny")

  # Create minimal fmri_model for testing
  sframe <- fmrihrf::sampling_frame(blocklens = c(50, 50), TR = 2)
  event_table <- data.frame(
    onset = c(10, 30, 60, 80),
    condition = factor(c("A", "B", "A", "B")),
    run = c(1, 1, 2, 2)
  )

  base_mod <- baseline_model(basis = "bs", degree = 3, sframe = sframe, intercept = "runwise")
  ev_mod <- event_model(
    onset ~ hrf(condition),
    data = event_table,
    block = ~ run,
    sampling_frame = sframe
  )
  fmri_mod <- fmri_model(ev_mod, base_mod)

  # Invalid term name should error
  expect_error(
    design_plot(fmri_mod, term_name = "nonexistent_term"),
    regexp = "term_name not found"
  )
})

# ============================================================================
# Integration test: visualization with real model
# ============================================================================

test_that("visualization integrates with model construction pipeline", {
  skip_on_cran()

  # Create a simple model
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  events <- data.frame(
    onset = c(10, 25, 50, 75),
    cond = factor(c("face", "scene", "face", "scene")),
    run = rep(1, 4)
  )

  ev_mod <- event_model(
    onset ~ hrf(cond),
    data = events,
    block = ~ run,
    sampling_frame = sframe
  )

  # Get design matrix and create correlation map
  dm <- as.matrix(design_matrix(ev_mod))
  if (ncol(dm) >= 2) {
    p <- fmrireg:::.correlation_map_common(dm)
    expect_s3_class(p, "ggplot")
  }

  # Test that regressors from the model can be plotted
  terms_list <- terms(ev_mod)
  if (length(terms_list) > 0) {
    first_term <- terms_list[[1]]
    reg <- first_term$regressors[[1]]
    if (inherits(reg, "Reg")) {
      p <- autoplot(reg)
      expect_s3_class(p, "ggplot")
    }
  }
})
