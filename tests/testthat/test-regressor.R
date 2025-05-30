# Test suite for Regressor functionality

library(testthat)
library(fmrireg)
library(ggplot2)

# Helper to get default HRF span safely
get_default_hrf_span <- function(hrf = HRF_SPMG1) {
  attr(hrf, "span") %||% 24 # Use 24 as a fallback
}

# === Test Core Reg Constructor (Tickets A-1, B-1, B-2, B-3) ===

test_that("Reg constructor creates object with correct class and structure", {
  ons <- c(10, 20, 30)
  reg <- fmrireg:::Reg(onsets = ons)
  expect_s3_class(reg, "Reg")
  expect_s3_class(reg, "list")
  expect_named(reg, c("onsets", "duration", "amplitude", "hrf", "span", "summate"))
  expect_equal(reg$onsets, ons)
  expect_true(inherits(reg$hrf, "HRF"))
  expect_equal(reg$duration, rep(0, 3))
  expect_equal(reg$amplitude, rep(1, 3))
  expect_true(is.numeric(reg$span))
  expect_true(is.logical(reg$summate))
})

test_that("Reg constructor recycles scalar duration and amplitude", {
  ons <- c(10, 20, 30)
  reg <- fmrireg:::Reg(onsets = ons, duration = 2, amplitude = 5)
  expect_equal(reg$duration, rep(2, 3))
  expect_equal(reg$amplitude, rep(5, 3))
})

test_that("Reg constructor handles vector duration and amplitude", {
  ons <- c(10, 20, 30)
  dur <- c(1, 2, 3)
  amp <- c(4, 5, 6)
  reg <- fmrireg:::Reg(onsets = ons, duration = dur, amplitude = amp)
  expect_equal(reg$duration, dur)
  expect_equal(reg$amplitude, amp)
})

test_that("Reg constructor errors on mismatched vector lengths", {
  ons <- c(10, 20, 30)
  # Use regexp to match the error structure flexibly
  expect_error(fmrireg:::Reg(onsets = ons, duration = c(1, 2)), 
               regexp = "Length mismatch for duration: got 2, expected 3")
  expect_error(fmrireg:::Reg(onsets = ons, amplitude = c(1, 2)), 
               regexp = "Length mismatch for amplitude: got 2, expected 3")
})

test_that("Reg constructor filters zero and NA amplitudes (B-3)", {
  ons <- c(10, 20, 30, 40, 50)
  amp <- c(1, 0, 2, NA, 3)
  dur <- c(1, 2, 3, 4, 5) # Provide vector to check filtering alignment
  
  reg <- fmrireg:::Reg(onsets = ons, amplitude = amp, duration = dur)
  
  expect_equal(reg$onsets, c(10, 30, 50))
  expect_equal(reg$amplitude, c(1, 2, 3))
  expect_equal(reg$duration, c(1, 3, 5))
})

test_that("Reg constructor handles all-zero/NA amplitudes", {
  ons <- c(10, 20, 30)
  amp <- c(0, NA, 0)
  reg <- regressor(onsets = ons, amplitude = amp)

  expect_equal(length(reg$onsets), 0)
  expect_equal(length(reg$duration), 0)
  expect_equal(length(reg$amplitude), 0)
})

test_that("Reg constructor handles empty onsets", {
  reg <- fmrireg:::Reg(numeric(0))
  
  expect_true(inherits(reg, "Reg"))
  expect_equal(length(reg$onsets), 0)
  expect_equal(length(reg$duration), 0) 
  expect_equal(length(reg$amplitude), 0)
})

test_that("Reg constructor handles NA onset (for null_regressor case)", {
   # Should not warn for the explicit NA case, should produce empty fields
   expect_warning(reg <- fmrireg:::Reg(onsets = NA, amplitude=0), regexp = NA) 
   expect_equal(length(reg$onsets), 0)
   expect_equal(length(reg$duration), 0)
   expect_equal(length(reg$amplitude), 0)
})

# === Test shift.Reg (Ticket B-4) ===

test_that("shift.Reg shifts onsets correctly", {
  ons <- c(10, 20, 30)
  reg_orig <- fmrireg:::Reg(onsets = ons)
  
  reg_plus5 <- shift(reg_orig, 5)
  expect_s3_class(reg_plus5, "Reg")
  expect_equal(reg_plus5$onsets, c(15, 25, 35))
  expect_equal(reg_plus5$hrf, reg_orig$hrf) # Other fields remain the same
  
  reg_minus3 <- shift(reg_orig, -3)
  expect_s3_class(reg_minus3, "Reg")
  expect_equal(reg_minus3$onsets, c(7, 17, 27))
  
  reg_zero <- shift(reg_orig, 0)
  expect_equal(reg_zero$onsets, reg_orig$onsets)
})

test_that("shift.Reg preserves original class structure", {
  ons <- c(10, 20, 30)
  # Create using the public wrapper to get 'regressor' class
  reg_orig_pub <- regressor(onsets = ons)
  expect_true(inherits(reg_orig_pub, "Reg"))
  
  reg_shifted <- shift(reg_orig_pub, 5)
  expect_true(inherits(reg_shifted, "Reg"))
  expect_equal(class(reg_shifted), class(reg_orig_pub))
})

test_that("shift.Reg handles empty regressors", {
  reg_empty1 <- fmrireg:::Reg(numeric(0))
  reg_empty2 <- fmrireg:::Reg(c(10, 20), amplitude = 0)
  
  expect_warning(shift(reg_empty1, 5), regexp = NA) # No warning expected
  expect_equal(length(shift(reg_empty1, 5)$onsets), 0)
  
  expect_warning(shift(reg_empty2, 5), regexp = NA)
  expect_equal(length(shift(reg_empty2, 5)$onsets), 0)
})

# === Test Legacy Wrappers (Ticket A-2') ===

test_that("regressor() facade works", {
  ons <- c(10, 20, 30)
  reg_pub <- regressor(onsets = ons, amplitude = 2)
  expect_true(inherits(reg_pub, "Reg"))
  expect_equal(reg_pub$amplitude, c(2, 2, 2)) # Check args passed through
})


# === Test evaluate.Reg (Tickets C-1, C-4) ===

# More evaluate tests needed, especially parity tests (F-3)

test_that("evaluate.Reg basic loop execution", {
    ons <- c(10, 30, 50)
    reg <- regressor(ons, HRF_SPMG1)
    grid <- seq(0, 70, by=2)
    
    # Currently only testing loop method works without error
    expect_no_error(val_loop <- evaluate(reg, grid, method="loop"))
    expect_true(is.numeric(val_loop))
    expect_equal(length(val_loop), length(grid))
    
    # Test multi-basis HRF
    reg_multi <- regressor(ons, HRF_SPMG3)
    expect_no_error(val_multi_loop <- evaluate(reg_multi, grid, method="loop"))
    expect_true(is.matrix(val_multi_loop))
    expect_equal(nrow(val_multi_loop), length(grid))
    expect_equal(ncol(val_multi_loop), 3)
})

test_that("evaluate.Reg handles empty regressor", {
    reg_empty <- fmrireg:::Reg(numeric(0))
    grid <- seq(0, 50, by=1)
    res <- evaluate(reg_empty, grid)
    expect_true(is.matrix(res))
    expect_equal(dim(res), c(length(grid), 1)) # Should default to nbasis=1
    expect_true(all(res == 0))
    
    # Empty with multi-basis hrf
    reg_empty_multi <- fmrireg:::Reg(numeric(0), hrf=HRF_SPMG3)
    res_multi <- evaluate(reg_empty_multi, grid)
    expect_true(is.matrix(res_multi))
    expect_equal(dim(res_multi), c(length(grid), 3))
    expect_true(all(res_multi == 0))
})

# === Test Evaluation Parity (Ticket F-3) ===

test_that("evaluate.Reg methods (loop, fft, conv) give consistent results", {
  ons <- c(10, 30, 50)
  grid <- seq(0, 70, by=2)
  
  # Scenario 1: Simple SPMG1
  reg_spmg1 <- regressor(ons, HRF_SPMG1)
  
  val_loop <- evaluate(reg_spmg1, grid, method = "loop")
  val_fft  <- evaluate(reg_spmg1, grid, method = "fft")
  val_conv <- evaluate(reg_spmg1, grid, method = "conv")
  
  expect_equal(val_fft, val_loop, tolerance = .1, info = "SPMG1: FFT vs Loop")
  expect_equal(val_conv, val_loop, tolerance = .1, info = "SPMG1: Conv vs Loop")

  # Scenario 2: Multi-basis SPMG3
  reg_spmg3 <- regressor(ons, HRF_SPMG3)
  
  val_multi_loop <- evaluate(reg_spmg3, grid, method = "loop")
  val_multi_fft  <- evaluate(reg_spmg3, grid, method = "fft")
  val_multi_conv <- evaluate(reg_spmg3, grid, method = "conv")
  
  expect_equal(val_multi_fft, val_multi_loop, tolerance = .1, info = "SPMG3: FFT vs Loop")
  expect_equal(val_multi_conv, val_multi_loop, tolerance = .1, info = "SPMG3: Conv vs Loop")

  # Scenario 3: With duration 
  reg_dur <- regressor(ons, HRF_SPMG1, duration = 5)
  
  val_dur_loop <- evaluate(reg_dur, grid, method = "loop")
  val_dur_fft  <- evaluate(reg_dur, grid, method = "fft")
  val_dur_conv <- evaluate(reg_dur, grid, method = "conv")
  
  # Tolerance might need adjustment depending on precision differences
  expect_equal(val_dur_fft, val_dur_loop, tolerance = .1, info = "Duration: FFT vs Loop")
  expect_equal(val_dur_conv, val_dur_loop, tolerance = .1, info = "Duration: Conv vs Loop")
  
  # Scenario 4: With amplitude modulation
  reg_amp <- regressor(ons, HRF_SPMG1, amplitude = c(1, -0.5, 2))
  
  val_amp_loop <- evaluate(reg_amp, grid, method = "loop")
  val_amp_fft  <- evaluate(reg_amp, grid, method = "fft")
  val_amp_conv <- evaluate(reg_amp, grid, method = "conv")
  
  expect_equal(val_amp_fft, val_amp_loop, tolerance = .1, info = "Amplitude: FFT vs Loop")
  expect_equal(val_amp_conv, val_amp_loop, tolerance = .1, info = "Amplitude: Conv vs Loop")
  
  # Note: More scenarios (Rconv, edge cases, sparse) still needed for full F-3 completion
})

# === Test Plotting & Printing (Tickets E-1, E-2) ===

test_that("print.Reg runs without error", {
  reg <- regressor(c(10, 20))
  expect_no_error(print(reg))
  
  reg_empty <- fmrireg:::Reg(numeric(0))
  expect_no_error(print(reg_empty))
})

test_that("autoplot.Reg runs without error", {
  skip_if_not_installed("ggplot2")
  reg <- regressor(c(10, 20))
  expect_s3_class(autoplot(reg), "ggplot")
  
  reg_multi <- regressor(c(10, 20), hrf=HRF_SPMG3)
  expect_s3_class(autoplot(reg_multi), "ggplot")
  
  reg_empty <- fmrireg:::Reg(numeric(0))
  expect_s3_class(autoplot(reg_empty), "ggplot")
})

test_that("generate an event model with one observation per level", {
  # Create mock data since lopdes is not available
  lopdes <- data.frame(
    WordPresentationOnset = rep(seq(1000, 10000, by = 1000), 5),
    Target = factor(rep(c("A", "B", "C", "D", "E"), each = 10)),
    Run = rep(1:5, each = 10)
  )
  
  sframe <- sampling_frame(blocklens=rep(401,5), TR=1.5)
  lopdes$onset <- lopdes$WordPresentationOnset/1000
  lopdes$Target <- factor(lopdes$Target)
  # Suppress warnings specifically for this call, as empty cells might trigger them
  expect_warning(ev <- event_model(onset ~ hrf(Target), data=lopdes, block= ~ Run, sampling_frame=sframe),
                 regexp = NA) # Expect warnings, but don't fail test if they occur
  expect_true(!is.null(ev))
})

test_that("evaluate.Reg errors when FFT size would be huge", {
  reg <- regressor(onsets = 0, hrf = HRF_SPMG1)
  grid <- seq(0, 1, by = 1)
  expect_error(
    evaluate(reg, grid, precision = 1e-6, method = "fft"),
    regexp = "FFT size"
  )
})
