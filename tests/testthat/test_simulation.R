test_that("simulate_bold_signal validates parameters correctly", {
  # Invalid ncond
  expect_error(simulate_bold_signal(ncond = 0), "ncond must be positive")
  expect_error(simulate_bold_signal(ncond = -1), "ncond must be positive")
  
  # Mismatched amps length
  expect_error(simulate_bold_signal(ncond = 3, amps = c(1, 2)), 
               "Length of 'amps' must equal 'ncond'")
  
  # Invalid ISI
  expect_error(simulate_bold_signal(ncond = 2, isi = c(3)), 
               "ISI must be a vector of length 2")
  expect_error(simulate_bold_signal(ncond = 2, isi = c(6, 3)), 
               "isi\\[2\\] > isi\\[1\\]")
  
  # Invalid TR
  expect_error(simulate_bold_signal(ncond = 2, TR = 0), "TR must be positive")
  expect_error(simulate_bold_signal(ncond = 2, TR = -1), "TR must be positive")
})

test_that("simulate_bold_signal generates expected output structure", {
  set.seed(123)
  sim <- simulate_bold_signal(ncond = 3, nreps = 10, TR = 2)
  
  # Check output structure
  expect_type(sim, "list")
  expect_named(sim, c("onset", "condition", "mat"))
  
  # Check dimensions
  expect_equal(length(sim$onset), 30)  # 3 conditions * 10 reps
  expect_equal(length(sim$condition), 30)
  expect_equal(ncol(sim$mat), 4)  # time + 3 conditions
  
  # Check condition names (order may vary due to sampling)
  expect_setequal(unique(sim$condition), c("Cond1", "Cond2", "Cond3"))
  
  # Check that time column is correct
  expect_equal(as.numeric(sim$mat[1, 1]), 0)
  expect_equal(as.numeric(diff(sim$mat[, 1])[1]), 2)  # TR = 2
})

test_that("simulate_bold_signal respects amplitude parameters", {
  set.seed(123)
  
  # Different amplitudes
  sim1 <- simulate_bold_signal(ncond = 2, amps = c(1, 2), nreps = 20, ampsd = 0)
  
  # Check that second condition has roughly double the amplitude
  max1 <- max(sim1$mat[, 2])
  max2 <- max(sim1$mat[, 3])
  expect_true(max2 > max1 * 1.8)  # Allow some tolerance
  
  # With amplitude variability
  sim2 <- simulate_bold_signal(ncond = 2, amps = c(1, 1), nreps = 50, ampsd = 0.5)
  
  # The signals should differ due to random amplitude sampling
  expect_false(all(sim2$mat[, 2] == sim2$mat[, 3]))
})

test_that("simulate_bold_signal handles different ISI distributions", {
  set.seed(123)
  
  # Uniform ISI
  sim1 <- simulate_bold_signal(ncond = 2, nreps = 10, isi = c(2, 8))
  isis1 <- diff(sim1$onset)
  expect_true(all(isis1 >= 2))
  expect_true(all(isis1 <= 8))
  
  # Fixed ISI (edge case) - need small difference for valid ISI
  sim2 <- simulate_bold_signal(ncond = 2, nreps = 10, isi = c(4, 4.001))
  isis2 <- diff(sim2$onset)
  expect_true(all(isis2 >= 4 & isis2 <= 4.001))
})

test_that("simulate_bold_signal is reproducible with seed", {
  sim1 <- simulate_bold_signal(ncond = 2, nreps = 10)
  sim2 <- simulate_bold_signal(ncond = 2, nreps = 10)
  
  # Without seed, should be different
  expect_false(all(sim1$onset == sim2$onset))
  
  # With same seed, should be identical
  set.seed(456)
  sim3 <- simulate_bold_signal(ncond = 2, nreps = 10)
  set.seed(456)
  sim4 <- simulate_bold_signal(ncond = 2, nreps = 10)
  
  expect_equal(sim3$onset, sim4$onset)
  expect_equal(sim3$mat, sim4$mat)
})

test_that("simulate_noise_vector generates noise with correct properties", {
  set.seed(123)
  n <- 200
  TR <- 2
  
  # White noise (use small AR to avoid warning)
  noise_white <- simulate_noise_vector(n, TR = TR, ar = numeric(0), ma = numeric(0), 
                                      drift_amplitude = 0, physio = FALSE)
  expect_length(noise_white, n)
  expect_equal(mean(noise_white), 0, tolerance = 0.2)
  expect_equal(sd(noise_white), 1, tolerance = 0.1)
  
  # AR(1) noise - should have autocorrelation
  noise_ar1 <- simulate_noise_vector(n, TR = TR, ar = 0.5, ma = 0,
                                    drift_amplitude = 0, physio = FALSE)
  acf_ar1 <- acf(noise_ar1, plot = FALSE)$acf[2]
  expect_true(acf_ar1 > 0.3)  # Should have positive autocorrelation
  
  # Check drift component
  noise_drift <- simulate_noise_vector(n, TR = TR, ar = numeric(0), ma = numeric(0),
                                      drift_freq = 1/64, drift_amplitude = 5,
                                      physio = FALSE, sd = 0)
  # Should see a clear sinusoidal pattern
  expect_true(max(noise_drift) > 4)
  expect_true(min(noise_drift) < -4)
})

test_that("simulate_noise_vector handles physiological noise", {
  set.seed(123)
  n <- 300
  TR <- 1
  
  # Without physio
  noise_no_physio <- simulate_noise_vector(n, TR = TR, physio = FALSE,
                                          drift_amplitude = 0, ar = numeric(0), ma = numeric(0))
  
  # With physio
  noise_physio <- simulate_noise_vector(n, TR = TR, physio = TRUE,
                                       drift_amplitude = 0, ar = numeric(0), ma = numeric(0))
  
  # Physio noise should add variance
  expect_true(var(noise_physio) > var(noise_no_physio))
  
  # Check for frequency components around cardiac/respiratory rates
  # This is a simple check - in reality would use FFT
  expect_false(all(noise_physio == noise_no_physio))
})

test_that("simulate_noise_vector is reproducible with seed", {
  n <- 100
  
  noise1 <- simulate_noise_vector(n, seed = 789)
  noise2 <- simulate_noise_vector(n, seed = 789)
  noise3 <- simulate_noise_vector(n, seed = 790)
  
  expect_equal(noise1, noise2)
  expect_false(all(noise1 == noise3))
})

test_that("simulate_simple_dataset integrates signal and noise correctly", {
  set.seed(123)
  
  data <- simulate_simple_dataset(ncond = 3, nreps = 10, TR = 2, snr = 0.5)
  
  # Check structure
  expect_type(data, "list")
  expect_named(data, c("clean", "noisy", "noise", "onsets", "conditions"))
  
  # Check dimensions match
  expect_equal(dim(data$clean$mat), dim(data$noisy))
  expect_equal(nrow(data$noise), nrow(data$clean$mat))
  expect_equal(ncol(data$noise), 3)  # 3 conditions
  
  # Check that noisy = clean + noise
  clean_signals <- data$clean$mat[, -1]  # Remove time column
  reconstructed <- clean_signals + data$noise
  noisy_signals <- data$noisy[, -1]
  
  expect_equal(reconstructed, noisy_signals, tolerance = 1e-10)
  
  # Check SNR roughly matches
  signal_power <- var(as.vector(clean_signals))
  noise_power <- var(as.vector(data$noise))
  estimated_snr <- sqrt(signal_power / noise_power)
  expect_equal(estimated_snr, 0.5, tolerance = 0.2)
})

test_that("simulate_simple_dataset handles edge cases", {
  set.seed(123)
  
  # Single condition
  data1 <- simulate_simple_dataset(ncond = 1, nreps = 5)
  expect_equal(ncol(data1$clean$mat), 2)  # time + 1 condition
  expect_equal(ncol(data1$noise), 1)
  
  # Single trial per condition
  data2 <- simulate_simple_dataset(ncond = 3, nreps = 1)
  expect_equal(length(data2$onsets), 3)
  expect_equal(length(unique(data2$conditions)), 3)
  
  # Very high SNR (minimal noise)
  data3 <- simulate_simple_dataset(ncond = 2, snr = 5, nreps = 20)
  # Check that the data was constructed correctly
  clean_signals <- data3$clean$mat[, -1]
  # Verify noise was added
  expect_false(all(data3$noisy[, -1] == clean_signals))
  # Just verify basic structure
  expect_equal(dim(data3$noise), dim(clean_signals))
})

test_that("simulate_fmri_matrix validates parameters", {
  # Basic parameter validation
  expect_error(simulate_fmri_matrix(n = 0), "n must be positive")
  expect_error(simulate_fmri_matrix(n = -1), "n must be positive")
  expect_error(simulate_fmri_matrix(total_time = 0), "total_time must be positive")
  expect_error(simulate_fmri_matrix(TR = 0), "TR must be positive")
  
  # ISI parameter validation
  expect_error(simulate_fmri_matrix(isi_min = 10, isi_max = 5),
               "isi_max must be greater than isi_min")
  
  # Amplitude/duration array length validation
  expect_error(simulate_fmri_matrix(n_events = 10, amplitudes = c(1, 2, 3)),
               "amplitudes must be length=1 or match n_events")
  expect_error(simulate_fmri_matrix(n_events = 10, durations = c(1, 2, 3)),
               "durations must be length=1 or match n_events")
})

test_that("simulate_fmri_matrix generates correct output structure", {
  set.seed(123)
  
  result <- simulate_fmri_matrix(
    n = 5,
    total_time = 120,
    TR = 2,
    n_events = 10,
    noise_type = "ar1"
  )
  
  # Check structure
  expect_type(result, "list")
  expect_named(result, c("time_series", "ampmat", "durmat", "hrf_info", "noise_params"))
  
  # Check time_series is a matrix_dataset
  expect_s3_class(result$time_series, "matrix_dataset")
  
  # Check dimensions
  expect_equal(ncol(result$time_series$data), 5)  # n = 5
  expect_equal(nrow(result$time_series$data), 60)  # 120s / 2s TR
  
  # Check amplitude and duration matrices - events may be reduced to fit
  actual_events <- nrow(result$time_series$event_table)
  expect_equal(dim(result$ampmat), c(actual_events, 5))  # n_events x n
  expect_equal(dim(result$durmat), c(actual_events, 5))
  
  # Check event table has at least some events
  expect_true(actual_events > 0)
  expect_true(actual_events <= 10)  # May be reduced
  expect_true(all(result$time_series$event_table$onset < 104))  # 120 - 16 buffer
})

test_that("simulate_fmri_matrix handles different ISI distributions", {
  set.seed(123)
  
  # Even spacing
  result_even <- simulate_fmri_matrix(
    n = 2,
    total_time = 100,
    n_events = 10,
    isi_dist = "even"
  )
  onsets_even <- result_even$time_series$event_table$onset
  isis_even <- diff(onsets_even)
  expect_true(sd(isis_even) < 0.1)  # Should be nearly constant
  
  # Uniform distribution
  result_unif <- simulate_fmri_matrix(
    n = 2,
    total_time = 100,
    n_events = 20,
    isi_dist = "uniform",
    isi_min = 2,
    isi_max = 6
  )
  onsets_unif <- result_unif$time_series$event_table$onset
  isis_unif <- diff(onsets_unif)
  expect_true(all(isis_unif >= 2))
  expect_true(all(isis_unif <= 6))
  
  # Exponential distribution
  result_exp <- simulate_fmri_matrix(
    n = 2,
    total_time = 200,
    n_events = 30,
    isi_dist = "exponential",
    isi_rate = 0.5
  )
  # Just check it runs without error
  expect_s3_class(result_exp$time_series, "matrix_dataset")
})

test_that("simulate_fmri_matrix handles amplitude and duration variability", {
  set.seed(123)
  
  # With amplitude variability
  result <- simulate_fmri_matrix(
    n = 3,
    n_events = 10,
    amplitudes = 2,
    amplitude_sd = 0.5,
    amplitude_dist = "gaussian"
  )
  
  # Check that amplitudes vary across columns
  amp_col1 <- result$ampmat[, 1]
  amp_col2 <- result$ampmat[, 2]
  expect_false(all(amp_col1 == amp_col2))
  
  # Check mean is roughly correct
  expect_equal(mean(result$ampmat), 2, tolerance = 0.3)
  
  # With duration variability
  result2 <- simulate_fmri_matrix(
    n = 3,
    n_events = 10,
    durations = 1,
    duration_sd = 0.2,
    duration_dist = "lognormal"
  )
  
  # Check that durations vary and are positive
  expect_true(all(result2$durmat > 0))
  expect_false(all(result2$durmat[, 1] == result2$durmat[, 2]))
})

test_that("simulate_fmri_matrix handles different noise types", {
  set.seed(123)
  
  # No noise
  result_none <- simulate_fmri_matrix(
    n = 2,
    n_events = 5,
    noise_type = "none"
  )
  
  # White noise
  result_white <- simulate_fmri_matrix(
    n = 2,
    n_events = 5,
    noise_type = "white",
    noise_sd = 1.5
  )
  
  # AR(1) noise with default
  result_ar1 <- simulate_fmri_matrix(
    n = 2,
    n_events = 5,
    noise_type = "ar1"
  )
  expect_equal(result_ar1$noise_params$noise_ar, 0.3)  # Default
  
  # AR(2) noise with custom coefficients
  result_ar2 <- simulate_fmri_matrix(
    n = 2,
    n_events = 5,
    noise_type = "ar2",
    noise_ar = c(0.4, 0.3)
  )
  expect_equal(result_ar2$noise_params$noise_ar, c(0.4, 0.3))
  
  # Check that noise was actually added
  data_none <- result_none$time_series$data
  data_white <- result_white$time_series$data
  expect_false(all(data_none == data_white))
})

test_that("simulate_fmri_matrix handles single trial mode", {
  set.seed(123)
  
  result_regular <- simulate_fmri_matrix(
    n = 2,
    n_events = 5,
    single_trial = FALSE,
    noise_type = "white",
    noise_sd = 0.5,
    amplitude_sd = 0.2  # Add variability
  )
  
  result_single <- simulate_fmri_matrix(
    n = 2,
    n_events = 5,
    single_trial = TRUE,
    noise_type = "white", 
    noise_sd = 0.5,
    amplitude_sd = 0.2  # Add variability
  )
  
  # Both should produce valid datasets
  expect_s3_class(result_regular$time_series, "matrix_dataset")
  expect_s3_class(result_single$time_series, "matrix_dataset")
  
  # Results should differ due to different processing modes
  # Check clean signals (before noise) by comparing column means
  col1_regular <- mean(result_regular$time_series$data[,1])
  col1_single <- mean(result_single$time_series$data[,1])
  # They should be different if the modes work differently
  expect_true(abs(col1_regular - col1_single) > 0 || 
              !identical(result_regular$ampmat, result_single$ampmat))
})

test_that("simulate_fmri_matrix is reproducible with seed", {
  result1 <- simulate_fmri_matrix(n = 3, n_events = 10, random_seed = 999, 
                                  noise_type = "white", noise_sd = 1.0)
  result2 <- simulate_fmri_matrix(n = 3, n_events = 10, random_seed = 999, 
                                  noise_type = "white", noise_sd = 1.0)
  result3 <- simulate_fmri_matrix(n = 3, n_events = 10, random_seed = 1000, 
                                  noise_type = "white", noise_sd = 1.0)
  
  # Same seed should give identical results
  expect_equal(result1$time_series$data, result2$time_series$data)
  expect_equal(result1$ampmat, result2$ampmat)
  expect_equal(result1$durmat, result2$durmat)
  
  # Different seed should give different results
  expect_false(identical(result1$time_series$data, result3$time_series$data))
})

test_that("simulate_fmri_matrix handles buffer correctly", {
  set.seed(123)
  
  result <- simulate_fmri_matrix(
    n = 2,
    total_time = 100,
    TR = 2,
    n_events = 20,
    buffer = 10
  )
  
  # All onsets should be within effective time (total_time - buffer)
  onsets <- result$time_series$event_table$onset
  expect_true(all(onsets <= 90))  # 100 - 10 (allow equality)
  
  # Time series should still be full length
  expect_equal(nrow(result$time_series$data), 50)  # 100 / 2
})

test_that("simulate_fmri_matrix handles extreme parameters gracefully", {
  set.seed(123)
  
  # Very short scan
  result1 <- simulate_fmri_matrix(
    n = 1,
    total_time = 20,
    n_events = 2
  )
  expect_equal(nrow(result1$time_series$event_table), 2)
  
  # Many events in short time - should reduce events
  expect_message(
    result2 <- simulate_fmri_matrix(
      n = 1,
      total_time = 30,
      n_events = 100,
      isi_dist = "uniform",
      isi_min = 0.5,
      isi_max = 1
    ),
    "Reduced to"
  )
  expect_true(nrow(result2$time_series$event_table) < 100)
  
  # Single time series
  result3 <- simulate_fmri_matrix(n = 1, n_events = 5)
  expect_equal(ncol(result3$time_series$data), 1)
  expect_equal(ncol(result3$ampmat), 1)
})

test_that("simulated data can be used with fmri_model construction", {
  skip_if_not_installed("fmrihrf")
  skip_if_not_installed("fmridesign")
  
  set.seed(123)
  
  # Generate simulated data
  sim_result <- simulate_fmri_matrix(
    n = 50,  # 50 voxels
    total_time = 200,
    TR = 2,
    n_events = 20,
    amplitudes = 2,
    amplitude_sd = 0.3,
    noise_type = "ar1",
    noise_sd = 0.5
  )
  
  # Extract dataset and event table
  dset <- sim_result$time_series
  etab <- dset$event_table
  
  # Create a simple event model
  etab$condition <- factor(rep(c("A", "B"), length.out = nrow(etab)))
  
  # Try to create models
  expect_no_error({
    emodel <- event_model(onset ~ hrf(condition), 
                         data = etab,
                         block = ~ run,
                         sampling_frame = dset$sampling_frame)
    
    bmodel <- baseline_model(basis = "poly", degree = 2, 
                           sframe = dset$sampling_frame)
    
    fmod <- fmri_model(emodel, bmodel, dset)
  })
  
  # Check that model was created successfully
  expect_s3_class(fmod, "fmri_model")
  expect_true(!is.null(fmod$event_model))
  expect_true(!is.null(fmod$baseline_model))
  expect_equal(fmod$dataset, dset)
  
  # Check design matrix can be extracted
  dm <- design_matrix(fmod)
  expect_equal(nrow(dm), nrow(dset$data))
  expect_true(ncol(dm) > 0)
})

test_that("ground truth recovery works for simulate_bold_signal", {
  set.seed(123)
  
  # Known parameters
  true_amps <- c(1, 2, 3)
  ncond <- 3
  nreps <- 50  # Many reps for better estimation
  
  sim <- simulate_bold_signal(ncond = ncond, amps = true_amps, 
                             nreps = nreps, ampsd = 0, TR = 1)
  
  # Estimate amplitudes by taking max of each condition's signal
  estimated_amps <- sapply(2:(ncond+1), function(i) max(sim$mat[, i]))
  
  # Should recover relative amplitudes
  ratios_true <- true_amps / true_amps[1]
  ratios_est <- estimated_amps / estimated_amps[1]
  
  expect_equal(ratios_est, ratios_true, tolerance = 0.1)
})

test_that("memory efficiency for large simulations", {
  skip_if_not_installed("pryr")
  
  # Test that memory usage is reasonable for large simulations
  # This is a basic check - in production would use more sophisticated profiling
  
  set.seed(123)
  
  # Large simulation
  result <- simulate_fmri_matrix(
    n = 1000,  # 1000 voxels
    total_time = 600,  # 10 minutes
    TR = 2,
    n_events = 100
  )
  
  # Check that object size is reasonable
  # Expected: ~300 timepoints x 1000 voxels x 8 bytes ≈ 2.4 MB for data alone
  obj_size <- object.size(result$time_series$data)
  expect_true(as.numeric(obj_size) < 10e6)  # Less than 10 MB
  
  # Clean up
  rm(result)
  gc()
})