test_that("benchmark dataset loading works", {
  skip_if_not_installed("fmrireg")
  
  # Test that we can list datasets
  datasets_info <- list_benchmark_datasets()
  expect_true(is.data.frame(datasets_info))
  expect_true(nrow(datasets_info) >= 5)
  expect_true(all(c("Dataset", "Description") %in% colnames(datasets_info)))
  
  # Test loading a specific dataset
  data <- load_benchmark_dataset("BM_Canonical_HighSNR")
  expect_true(is.list(data))
  expect_true("Y_noisy" %in% names(data))
  expect_true("event_onsets" %in% names(data))
  expect_true("condition_labels" %in% names(data))
  expect_true("true_betas_condition" %in% names(data))
  
  # Check dimensions are reasonable
  expect_true(is.matrix(data$Y_noisy))
  expect_true(ncol(data$Y_noisy) == 100)  # 100 voxels
  expect_true(nrow(data$Y_noisy) > 100)   # Should have many time points
  
  # Test metadata loading
  metadata <- load_benchmark_dataset("metadata")
  expect_true(is.list(metadata))
  expect_true("creation_date" %in% names(metadata))
  expect_true("description" %in% names(metadata))
})

test_that("benchmark dataset summary works", {
  skip_if_not_installed("fmrireg")
  
  summary_info <- get_benchmark_summary("BM_Canonical_HighSNR")
  expect_true(is.list(summary_info))
  expect_true("dimensions" %in% names(summary_info))
  expect_true("experimental_design" %in% names(summary_info))
  expect_true("hrf_information" %in% names(summary_info))
  
  # Check dimensions
  dims <- summary_info$dimensions
  expect_equal(dims$n_voxels, 100)
  expect_equal(dims$n_conditions, 3)
  expect_true(dims$n_events > 0)
})

test_that("design matrix creation works", {
  skip_if_not_installed("fmrireg")
  
  # Test design matrix creation
  X <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", fmrihrf::HRF_SPMG1)
  expect_true(is.matrix(X))
  expect_true(ncol(X) == 4)  # 3 conditions + intercept
  expect_true(nrow(X) > 100)  # Should have many time points
  expect_true("Intercept" %in% colnames(X))
  
  # Test without intercept
  X_no_int <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", fmrihrf::HRF_SPMG1, 
                                                  include_intercept = FALSE)
  expect_true(ncol(X_no_int) == 3)  # 3 conditions only
  expect_false("Intercept" %in% colnames(X_no_int))
})

test_that("benchmark dimensions and scan metadata are internally consistent", {
  datasets <- load_benchmark_dataset("all")

  for (dataset in datasets) {
    expect_equal(dataset$run_length, nrow(dataset$Y_noisy))
    expect_equal(dataset$core_data_args$run_length, nrow(dataset$Y_noisy))
    expect_length(dataset$event_durations, length(dataset$event_onsets))
    expect_false(anyNA(dataset$event_durations))
    expect_equal(ncol(dataset$true_betas_condition), ncol(dataset$Y_noisy))
    expect_equal(ncol(dataset$true_amplitudes_trial), ncol(dataset$Y_noisy))

    if (!is.null(dataset$Y_clean)) {
      expect_equal(nrow(dataset$Y_clean), nrow(dataset$Y_noisy))
    }
    if (!is.null(dataset$true_durations_trial)) {
      expect_equal(ncol(dataset$true_durations_trial), ncol(dataset$Y_noisy))
    }
    if (!is.null(dataset$X_list_true_hrf)) {
      expect_true(all(vapply(
        dataset$X_list_true_hrf,
        NROW,
        integer(1)
      ) == nrow(dataset$Y_noisy)))
    }
  }
})

test_that("all complete canonical oracles recover their declared betas", {
  complete_oracles <- c(
    "BM_Canonical_HighSNR",
    "BM_Canonical_LowSNR"
  )

  for (dataset_name in complete_oracles) {
    dataset <- load_benchmark_dataset(dataset_name)
    X <- create_design_matrix_from_benchmark(
      dataset_name,
      dataset$true_hrf_parameters$hrf_object
    )
    estimated <- qr.solve(X, dataset$Y_clean)[-1, , drop = FALSE]

    expect_identical(dataset$oracle_contract$level, "complete_condition_beta")
    expect_lt(
      max(abs(estimated - dataset$true_betas_condition)),
      dataset$oracle_contract$tolerance
    )
  }
})

test_that("every declared HRF recipe reconstructs exactly on a dense grid", {
  metadata <- load_benchmark_dataset("metadata")
  recipes <- metadata$hrf_recipes
  expect_named(
    recipes,
    c("HRF_SPMG1", "HRF_SPMG2", "HRF_SPMG3", "variant1", "variant2")
  )

  clean_base <- function(hrf) {
    fmrihrf::as_hrf(
      function(time) hrf(time),
      name = attr(hrf, "name"),
      nbasis = fmrihrf::nbasis(hrf),
      span = attr(hrf, "span")
    )
  }
  expected <- list(
    HRF_SPMG1 = fmrihrf::HRF_SPMG1,
    HRF_SPMG2 = fmrihrf::HRF_SPMG2,
    HRF_SPMG3 = fmrihrf::HRF_SPMG3,
    variant1 = fmrihrf::gen_hrf(
      clean_base(fmrihrf::HRF_SPMG1),
      lag = 1, width = 1.2, normalize = TRUE
    ),
    variant2 = fmrihrf::gen_hrf(
      clean_base(fmrihrf::HRF_SPMG1),
      lag = -0.5, width = 0.8, normalize = TRUE
    )
  )
  dense_grid <- seq(0, 40, by = 0.01)
  required_fields <- c(
    "schema_version", "recipe_name", "generator", "base_hrf_name",
    "lag", "width", "precision", "half_life", "summate", "normalize",
    "name_override", "span_override", "fmrihrf_version"
  )

  for (recipe_name in names(recipes)) {
    recipe <- recipes[[recipe_name]]
    expect_true(all(required_fields %in% names(recipe)), info = recipe_name)
    expect_identical(recipe$recipe_name, recipe_name)
    expect_identical(
      recipe$fmrihrf_version,
      as.character(utils::packageVersion("fmrihrf"))
    )
    reconstructed <- fmrireg:::.reconstruct_hrf_object(recipe)
    expect_equal(
      fmrihrf::evaluate(reconstructed, dense_grid),
      fmrihrf::evaluate(expected[[recipe_name]], dense_grid),
      tolerance = 1e-12,
      info = recipe_name
    )
  }

  collect_groups <- function(parameters) {
    if (is.list(parameters) && "hrf_object_name" %in% names(parameters)) {
      return(list(parameters))
    }
    unlist(lapply(parameters, collect_groups), recursive = FALSE)
  }
  for (dataset in load_benchmark_dataset("all")) {
    for (group in collect_groups(dataset$true_hrf_parameters)) {
      expect_identical(group$hrf_recipe, recipes[[group$hrf_object_name]])
      expect_equal(
        fmrihrf::evaluate(group$hrf_object, dense_grid),
        fmrihrf::evaluate(expected[[group$hrf_object_name]], dense_grid),
        tolerance = 1e-12,
        info = group$hrf_object_name
      )
    }
  }
})

test_that("unknown or inconsistent HRF recipes fail closed", {
  recipe <- load_benchmark_dataset("metadata")$hrf_recipes$variant1

  unknown_name <- recipe
  unknown_name$recipe_name <- "unknown"
  expect_error(
    fmrireg:::.reconstruct_hrf_object(unknown_name),
    "Unknown benchmark HRF recipe name"
  )

  unknown_base <- recipe
  unknown_base$base_hrf_name <- "HRF_NOT_REAL"
  expect_error(
    fmrireg:::.reconstruct_hrf_object(unknown_base),
    "Unknown base HRF"
  )

  incompatible_version <- recipe
  incompatible_version$fmrihrf_version <- "0.0.0"
  expect_error(
    fmrireg:::.reconstruct_hrf_object(incompatible_version),
    "requires fmrihrf 0.0.0"
  )

  expect_error(
    fmrireg:::.reconstruct_hrf_object("variant1"),
    "recipe must be a named list"
  )
})

test_that("wrong HRF is worse by a scale-sensitive metric", {
  dataset <- load_benchmark_dataset("BM_Canonical_HighSNR")
  X_correct <- create_design_matrix_from_benchmark(
    "BM_Canonical_HighSNR",
    fmrihrf::HRF_SPMG1
  )
  X_wrong <- create_design_matrix_from_benchmark(
    "BM_Canonical_HighSNR",
    fmrihrf::HRF_GAUSSIAN
  )

  beta_correct <- qr.solve(X_correct, dataset$Y_noisy)[-1, , drop = FALSE]
  beta_wrong <- qr.solve(X_wrong, dataset$Y_noisy)[-1, , drop = FALSE]
  rmse_correct <- sqrt(mean((beta_correct - dataset$true_betas_condition)^2))
  rmse_wrong <- sqrt(mean((beta_wrong - dataset$true_betas_condition)^2))

  expect_lt(rmse_correct, rmse_wrong)
})

test_that("performance evaluation works", {
  skip_if_not_installed("fmrireg")
  
  # Load dataset and create simple test
  data <- load_benchmark_dataset("BM_Canonical_HighSNR")
  true_betas <- data$true_betas_condition
  
  # Create some fake estimated betas with more realistic variation
  # Add noise to true betas AND add some spatial variation to avoid zero variance
  set.seed(123)
  estimated_betas <- true_betas + matrix(rnorm(prod(dim(true_betas)), 0, 0.1), 
                                        nrow = nrow(true_betas))
  
  # Add some spatial variation to make correlations meaningful
  for (i in 1:nrow(estimated_betas)) {
    spatial_trend <- seq(-0.1, 0.1, length.out = ncol(estimated_betas))
    estimated_betas[i, ] <- estimated_betas[i, ] + spatial_trend
  }
  
  # Evaluate performance
  performance <- evaluate_method_performance("BM_Canonical_HighSNR", 
                                            estimated_betas, 
                                            "Test_Method")
  
  expect_true(is.list(performance))
  expect_true("overall_metrics" %in% names(performance))
  expect_true("condition_metrics" %in% names(performance))
  expect_true("voxel_metrics" %in% names(performance))
  
  # Check that correlation is reasonable (should be positive but not perfect due to added noise)
  expect_true(performance$overall_metrics$correlation > 0.5)
  expect_true(performance$overall_metrics$correlation < 1.0)
  expect_true(performance$overall_metrics$mse > 0)
})

test_that("all benchmark datasets can be loaded", {
  skip_if_not_installed("fmrireg")
  
  # Get list of all datasets
  datasets_info <- list_benchmark_datasets()
  dataset_names <- datasets_info$Dataset
  
  # Try to load each dataset
  for (dataset_name in dataset_names) {
    data <- load_benchmark_dataset(dataset_name)
    expect_true(is.list(data))
    expect_true("Y_noisy" %in% names(data))
    expect_true("event_onsets" %in% names(data))
    expect_true("condition_labels" %in% names(data))
    
    # Check that dimensions are consistent with simulation
    # The simulation uses ceiling(total_time / TR) for the number of timepoints
    expected_timepoints <- ceiling(data$total_time / data$TR)
    expect_equal(nrow(data$Y_noisy), expected_timepoints)
    expect_equal(length(data$event_onsets), length(data$condition_labels))
  }
})

test_that("error handling works correctly", {
  skip_if_not_installed("fmrireg")
  
  # Test invalid dataset name
  expect_error(load_benchmark_dataset("NonExistent_Dataset"))
  expect_error(get_benchmark_summary("NonExistent_Dataset"))
  
  # Test performance evaluation with wrong dimensions
  data <- load_benchmark_dataset("BM_Canonical_HighSNR")
  wrong_betas <- matrix(1, nrow = 2, ncol = 50)  # Wrong dimensions
  expect_error(evaluate_method_performance("BM_Canonical_HighSNR", wrong_betas, "Test"))
}) 
