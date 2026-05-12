test_that("fmridesign integration works correctly", {
  skip_if_not_installed("fmridesign")
  
  # Create test data
  event_data <- data.frame(
    condition = factor(c("A", "B", "A", "B", "A", "B")),
    onsets = c(1, 15, 30, 45, 60, 75),
    run = c(1, 1, 1, 2, 2, 2)
  )
  
  # Create sampling frame
  sframe <- sampling_frame(blocklens = c(50, 50), TR = 2)
  
  # Test event model creation
  ev_model <- event_model(
    onsets ~ hrf(condition),
    data = event_data,
    block = ~run,
    sampling_frame = sframe
  )
  
  expect_s3_class(ev_model, "event_model")
  expect_true(!is.null(ev_model$design_matrix))
  
  # Test baseline model creation
  base_model <- baseline_model(
    basis = "bs",
    degree = 3,
    sframe = sframe
  )
  
  expect_s3_class(base_model, "baseline_model")
  
  # Test fmri dataset creation
  X <- matrix(rnorm(100 * 100), 100, 100)  # 100 time points, 100 voxels
  dset <- matrix_dataset(X, TR = 2, run_length = c(50, 50), 
                         event_table = event_data)
  
  # Test fmri model creation
  fmri_mod <- fmri_model(ev_model, base_model, dset)
  expect_s3_class(fmri_mod, "fmri_model")
  
  # Test design matrix extraction
  dm <- design_matrix(fmri_mod)
  expect_equal(nrow(dm), 100)
  expect_true(ncol(dm) > 0)
  
  # Test model fitting
  fit <- fmri_lm(fmri_mod, dset, strategy = "chunkwise", nchunks = 2)
  expect_s3_class(fit, "fmri_lm")
  
  # Test coefficient extraction
  beta <- coef(fit, include_baseline = TRUE)
  expect_equal(ncol(beta), ncol(dm))
  expect_equal(nrow(beta), 100)  # number of voxels
})

test_that("event model formula interface works", {
  skip_if_not_installed("fmridesign")
  
  event_data <- data.frame(
    cond = factor(c("A", "B", "C", "A", "B", "C")),
    rt = c(0.5, 0.7, 0.6, 0.8, 0.9, 0.7),
    onsets = c(1, 10, 20, 30, 40, 50),
    run = rep(1, 6)
  )
  
  sframe <- sampling_frame(blocklens = 60, TR = 2)
  
  # Test multiple terms
  ev_model <- event_model(
    onsets ~ hrf(cond) + hrf(rt),
    data = event_data,
    block = ~run,
    sampling_frame = sframe
  )
  
  expect_true(length(ev_model$terms) == 2)
  expect_true("cond" %in% names(ev_model$terms))
  expect_true("rt" %in% names(ev_model$terms))
})

test_that("contrast weights are properly exported", {
  skip_if_not_installed("fmridesign")
  
  event_data <- data.frame(
    condition = factor(c("A", "B", "A", "B")),
    onsets = c(1, 10, 20, 30),
    run = c(1, 1, 1, 1)
  )
  
  sframe <- sampling_frame(blocklens = 40, TR = 2)
  
  ev_model <- event_model(
    onsets ~ hrf(condition),
    data = event_data,
    block = ~run,
    sampling_frame = sframe
  )
  
  # Test that contrast_weights method is available
  expect_true(exists("contrast_weights"))
  cw <- contrast_weights(ev_model)
  # contrast_weights may return NULL for simple event models without explicit contrasts
  expect_true(is.null(cw) || is.list(cw))
})

test_that("list-built hrfspec subsets match formula interface", {
  skip_if_not_installed("fmridesign")

  keep_a <- function(x) x == "A"

  event_data <- data.frame(
    onset = c(1, 10, 20, 30),
    condition = factor(c("A", "B", "A", "B")),
    duration = 1,
    run = rep(1, 4)
  )

  sframe <- sampling_frame(blocklens = 40, TR = 2)

  formula_eq <- event_model(
    onset ~ hrf(condition, id = "stim", subset = condition == "A"),
    data = event_data,
    block = ~run,
    sampling_frame = sframe,
    durations = event_data$duration
  )

  list_eq <- event_model(
    list(stim = hrf(condition, id = "stim", subset = condition == "A")),
    data = event_data,
    block = ~run,
    sampling_frame = sframe,
    durations = event_data$duration
  )

  expect_equal(names(formula_eq$terms), names(list_eq$terms))
  expect_equal(colnames(design_matrix(formula_eq)), colnames(design_matrix(list_eq)))
  expect_equal(
    unname(as.matrix(design_matrix(formula_eq))),
    unname(as.matrix(design_matrix(list_eq)))
  )

  formula_helper <- event_model(
    onset ~ hrf(condition, id = "stim", subset = keep_a(condition)),
    data = event_data,
    block = ~run,
    sampling_frame = sframe,
    durations = event_data$duration
  )

  list_helper <- event_model(
    list(stim = hrf(condition, id = "stim", subset = keep_a(condition))),
    data = event_data,
    block = ~run,
    sampling_frame = sframe,
    durations = event_data$duration
  )

  expect_equal(names(formula_helper$terms), names(list_helper$terms))
  expect_equal(colnames(design_matrix(formula_helper)), colnames(design_matrix(list_helper)))
  expect_equal(
    unname(as.matrix(design_matrix(formula_helper))),
    unname(as.matrix(design_matrix(list_helper)))
  )
})
