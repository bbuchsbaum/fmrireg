library(testthat)
library(fmrireg)
library(rlang)

test_that("create_event_model handles various variable specifications", {
  # Setup test data
  set.seed(123)
  
  # Define block structure properly
  blocklens <- c(5, 5)  # 5 scans per block
  TR <- 2
  
  # Create proper sampling frame
  sf <- sampling_frame(blocklens = blocklens, TR = TR)
  
  # Create events with proper block-wise onsets
  events <- data.frame(
    x = rnorm(10),
    y = rnorm(10),
    condition = factor(rep(c("A", "B"), each=5)),
    modulator = seq(0, 1, length.out=10),
    onset = rep(seq(0, 8, by=2), 2),  # Onsets reset for each block
    block = rep(1:2, each=5)
  )
  
  # Extract onsets and blocks
  onsets <- events$onset
  blocks <- events$block
  
  # Test 1: Character string variables
  event_terms <- list(
    list(
      variables = c("x", "y"),
      hrf = "spmg1"
    )
  )
  
  model1 <- create_event_model(
    event_terms = event_terms,
    events = events,
    onsets = onsets,
    block = blocks,
    sampling_frame = sf
  )
  
  expect_s3_class(model1, "event_model")
  expect_equal(length(terms(model1)), 1)
  
  # Verify block structure
  expect_equal(length(unique(model1$model_spec$blockids)), 2)
  expect_equal(model1$model_spec$blockvals, blocks)
  
  # Test that sampling frame matches block structure
  expect_equal(length(sf$blocklens), 2)
  expect_equal(sum(sf$blocklens), nrow(events))
  
  # Test 2: Mix of character strings and expressions
  event_terms <- list(
    list(
      variables = c("condition"),
      hrf = "spmg1"
    ),
    list(
      variables = exprs(Poly(modulator, 2)),
      hrf = "gamma"
    )
  )
  
  model2 <- create_event_model(
    event_terms = event_terms,
    events = events,
    onsets = onsets,
    block = blocks,
    sampling_frame = sf
  )
  
  expect_s3_class(model2, "event_model")
  expect_equal(length(terms(model2)), 2)
  
  # Test 3: Multiple expressions in one term
  event_terms <- list(
    list(
      variables = exprs(condition, Poly(modulator, 2)),
      hrf = "spmg1"
    )
  )
  
  model3 <- create_event_model(
    event_terms = event_terms,
    events = events,
    onsets = onsets,
    block = blocks,
    sampling_frame = sf
  )
  
  expect_s3_class(model3, "event_model")
  expect_equal(length(terms(model3)), 1)
  
  # Test 4: Complex model with multiple terms and transformations
  event_terms <- list(
    list(
      variables = c("condition"),
      hrf = "spmg1"
    ),
    list(
      variables = exprs(Poly(x, 2)),
      hrf = "gamma"
    ),
    list(
      variables = exprs(condition, modulator),
      hrf = "spmg1"
    ),
    list(
      variables = exprs(Standardized(y)),
      hrf = "spmg1"
    )
  )
  
  model4 <- create_event_model(
    event_terms = event_terms,
    events = events,
    onsets = onsets,
    block = blocks,
    sampling_frame = sf
  )
  
  expect_s3_class(model4, "event_model")
  expect_equal(length(terms(model4)), 4)
  
  # Test 5: Error handling for invalid variable names
  event_terms <- list(
    list(
      variables = c("nonexistent"),
      hrf = "spmg1"
    )
  )
  
  expect_error(
    create_event_model(
      event_terms = event_terms,
      events = events,
      onsets = onsets,
      block = blocks,
      sampling_frame = sf
    )
  )
  
  # Test 6: Mix of basis functions
  event_terms <- list(
    list(
      variables = exprs(BSpline(x, 3)),
      hrf = "spmg1"
    ),
    list(
      variables = exprs(Poly(y, 2), Standardized(modulator)),
      hrf = "gamma"
    )
  )
  
  model6 <- create_event_model(
    event_terms = event_terms,
    events = events,
    onsets = onsets,
    block = blocks,
    sampling_frame = sf
  )
  
  expect_s3_class(model6, "event_model")
  expect_equal(length(terms(model6)), 2)
  
  # Test 7: Check design matrix dimensions
  dm <- design_matrix(model6)
  expect_true(nrow(dm) == length(onsets))
  
  # Test 8: Check formula reconstruction
  expect_true(inherits(model6$model_spec$formula, "formula"))
})

test_that("create_event_model handles HRF specifications correctly", {
  # Setup test data
  set.seed(123)
  events <- data.frame(
    x = rnorm(10),
    y = rnorm(10)
  )
  onsets <- seq(0, 90, length.out=10)
  blocks <- rep(1:2, each=5)
  sf <- sampling_frame(c(5,5), TR=2)
  
  # Test different HRF specifications
  event_terms <- list(
    list(
      variables = c("x"),
      hrf = "spmg1"
    ),
    list(
      variables = c("y"),
      hrf = list(hrf = "gamma", parameters = list(shape = 6, rate = 1))
    )
  )
  
  model <- create_event_model(
    event_terms = event_terms,
    events = events,
    onsets = onsets,
    block = blocks,
    sampling_frame = sf
  )
  
  expect_s3_class(model, "event_model")
  expect_equal(length(terms(model)), 2)
})

test_that("create_event_model handles subset specifications", {
  # Setup test data
  set.seed(123)
  events <- data.frame(
    x = rnorm(10),
    y = rnorm(10),
    condition = factor(rep(c("A", "B"), each=5))
  )
  onsets <- seq(0, 90, length.out=10)
  blocks <- rep(1:2, each=5)
  sf <- sampling_frame(c(5,5), TR=2)
  
  # Test with subset
  event_terms <- list(
    list(
      variables = c("x", "y"),
      hrf = "spmg1",
      subset = quote(condition == "A")
    )
  )
  
  model <- create_event_model(
    event_terms = event_terms,
    events = events,
    onsets = onsets,
    block = blocks,
    sampling_frame = sf
  )
  
  expect_s3_class(model, "event_model")
})

test_that("create_event_model handles block-wise onsets correctly", {
  # Setup test data with explicit block structure
  blocklens <- c(5, 5)
  TR <- 2
  sf <- sampling_frame(blocklens = blocklens, TR = TR)
  
  events <- data.frame(
    x = rnorm(10),
    onset = rep(seq(0, 8, by=2), 2),  # Onsets reset for each block
    block = rep(1:2, each=5)
  )
  
  # Test that global onsets are computed correctly
  global_ons <- global_onsets(sf, events$onset, events$block)
  expect_equal(length(global_ons), nrow(events))
  expect_true(all(diff(global_ons[1:5]) > 0))  # First block onsets increase
  expect_true(all(diff(global_ons[6:10]) > 0))  # Second block onsets increase
  
  # Create model
  event_terms <- list(
    list(
      variables = c("x"),
      hrf = "spmg1"
    )
  )
  
  model <- create_event_model(
    event_terms = event_terms,
    events = events,
    onsets = events$onset,
    block = events$block,
    sampling_frame = sf
  )
  
  expect_s3_class(model, "event_model")
  
  # Check that design matrix has correct dimensions based on sampling frame
  dm <- design_matrix(model)
  expect_equal(nrow(dm), sum(blocklens))
})

test_that("event_model handles ~1 block formula correctly", {
  # Setup test data
  set.seed(123)
  
  # Create simple dataset
  events <- data.frame(
    x = rnorm(10),
    onset = seq(0, 18, by=2)
  )
  
  # Create sampling frame for a single block
  sf <- sampling_frame(blocklens = 10, TR = 2)
  
  # Create event terms
  event_terms <- list(
    list(
      variables = c("x"),
      hrf = "spmg1"
    )
  )
  
  # Model with block = ~1
  model1 <- create_event_model(
    event_terms = event_terms,
    events = events,
    onsets = events$onset,
    block = ~1,
    sampling_frame = sf
  )
  
  # Model with explicit block vector
  model2 <- create_event_model(
    event_terms = event_terms,
    events = events,
    onsets = events$onset,
    block = rep(1, nrow(events)),
    sampling_frame = sf
  )
  
  # Tests
  expect_s3_class(model1, "event_model")
  
  # Check that all events are in a single block
  expect_equal(length(unique(model1$model_spec$blockids)), 1)
  expect_equal(model1$model_spec$blockids, rep(1, nrow(events)))
  
  # Compare with explicit vector model
  expect_equal(model1$model_spec$blockids, model2$model_spec$blockids)
  expect_equal(model1$model_spec$blockvals, model2$model_spec$blockvals)
  
  # Design matrices should be identical
  dm1 <- design_matrix(model1)
  dm2 <- design_matrix(model2)
  expect_equal(dm1, dm2)
  
  # Also test with event_model.formula directly
  form <- onset ~ hrf(x)
  model3 <- event_model(form, events, block = ~1, sampling_frame = sf)
  expect_s3_class(model3, "event_model")
  expect_equal(length(unique(model3$model_spec$blockids)), 1)
})
