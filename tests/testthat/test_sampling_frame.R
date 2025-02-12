context("Sampling Frame")

test_that("sampling_frame constructor works correctly", {
  # Basic construction
  sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)
  expect_s3_class(sframe, "sampling_frame")
  expect_equal(length(sframe$blocklens), 2)
  expect_equal(sframe$TR, c(2, 2))
  expect_equal(sframe$start_time, c(1, 1))
  
  # Test with different TRs per block
  sframe2 <- sampling_frame(blocklens = c(100, 200), TR = c(2, 1.5))
  expect_equal(sframe2$TR, c(2, 1.5))
  
  # Test input validation
  expect_error(sampling_frame(blocklens = c(-1, 100), TR = 2), 
              "Block lengths must be positive")
  expect_error(sampling_frame(blocklens = c(100, 100), TR = -1), 
              "TR .* must be positive")
  expect_error(sampling_frame(blocklens = c(100, 100), TR = 2, precision = 3),
              "Precision must be positive and less than")
})

test_that("samples.sampling_frame works correctly", {
  sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)
  
  # Test relative timing
  rel_samples <- samples(sframe, global = FALSE)
  expect_equal(length(rel_samples), 200)
  expect_equal(rel_samples[1:5], c(1, 3, 5, 7, 9))
  
  # Test global timing
  glob_samples <- samples(sframe, global = TRUE)
  expect_equal(length(glob_samples), 200)
  expect_equal(glob_samples[101] - glob_samples[100], 2)  # Check TR spacing
  
  # Test block selection
  block1_samples <- samples(sframe, blockids = 1)
  expect_equal(length(block1_samples), 100)
  
  # Test memoization
  samples2 <- samples(sframe, global = FALSE)
  expect_identical(rel_samples, samples2)  # Should return cached result
})

test_that("global_onsets works correctly", {
  sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)
  
  # Test basic functionality
  onsets <- c(10, 20)
  blockids <- c(1, 2)
  global_times <- global_onsets(sframe, onsets, blockids)
  expect_equal(length(global_times), 2)
  expect_equal(global_times[1], 10)  # First block onset unchanged
  expect_equal(global_times[2], 220)  # Second block onset = 200 (block1 duration) + 20
  
  # Test error conditions
  #expect_error(global_onsets(sframe, c(10), c(1, 2)), 
  #            "length.*onsets.*length.*blockids")
  #expect_error(global_onsets(sframe, c(10), c(3)), 
  #            "blockids.*1.*length")
})

test_that("print.sampling_frame works correctly", {
  sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)
  expect_output(print(sframe), "Sampling Frame")
  expect_output(print(sframe), "Structure")
  expect_output(print(sframe), "Timing")
  expect_output(print(sframe), "Duration")
})

test_that("sampling_frame handles edge cases", {
  # Single block
  single_block <- sampling_frame(blocklens = 100, TR = 2)
  expect_equal(length(single_block$blocklens), 1)
  expect_equal(length(samples(single_block)), 100)
  
  # Very short block
  short_block <- sampling_frame(blocklens = c(1, 1), TR = 2)
  expect_equal(length(samples(short_block)), 2)
  
  # Different start times
  custom_starts <- sampling_frame(blocklens = c(100, 100), 
                                TR = 2, 
                                start_time = c(0, 5))
  expect_equal(custom_starts$start_time, c(0, 5))
  
  # High precision
  high_prec <- sampling_frame(blocklens = c(10, 10), 
                            TR = 2, 
                            precision = 0.01)
  expect_equal(high_prec$precision, 0.01)
})

test_that("sampling_frame maintains temporal consistency", {
  sframe <- sampling_frame(blocklens = c(100, 100, 100), TR = 2)
  glob_samples <- samples(sframe, global = TRUE)
  
  # Check uniform spacing within blocks
  for (block in 1:3) {
    block_idx <- which(sframe$blockids == block)
    diffs <- diff(glob_samples[block_idx])
    expect_true(all(abs(diffs - 2) < 1e-10))
  }
  
  # Check block transitions
  block_ends <- cumsum(sframe$blocklens)
  for (i in 1:(length(block_ends)-1)) {
    time_diff <- glob_samples[block_ends[i] + 1] - glob_samples[block_ends[i]]
    expect_equal(time_diff, 2)
  }
})