# tests/testthat/test_fmridsl.R
library(bidser)
library(fmrireg)
library(testthat)

context("fmriDSL Loading and Model Building")

test_that("Minimal fmriDSL config loads and builds fmri_model", {

  # 1. Setup Mock BIDS Data using bidser
  # -------------------------------------------------
  # Define participants
  participants_df <- tibble::tibble(participant_id = c("01"), age = c(28))

  # Define the file structure (raw data only for simplicity)
  # Suffixes are BIDS entities like 'bold', 'events', 'T1w'. Extensions added by create_mock_bids.
  file_structure_df <- tibble::tribble(
    ~subid, ~session, ~datatype, ~task,     ~run, ~suffix,  ~fmriprep,
    "01",   NA,       "anat",    NA,        NA,   "T1w",    FALSE,
    "01",   NA,       "func",    "taskA",   "01", "bold",   FALSE,
    "01",   NA,       "func",    "taskA",   "01", "events", FALSE
  )

  # Define event data content for the events.tsv file
  # Predict the relative path create_mock_bids will use
  event_filename_1 <- bidser:::generate_bids_filename(subid = "01", task = "taskA", run = "01", suffix = "events.tsv")
  event_path_1 <- file.path("sub-01", "func", event_filename_1)

  event_data_list <- list()
  event_data_list[[event_path_1]] <- tibble::tibble(
    onset = c(1.0, 5.0, 10.0, 15.0),
    duration = c(0.5, 0.5, 0.5, 0.5),
    trial_type = c("condA", "condB", "condA", "condB"),
    response_time = c(0.8, 0.9, 0.75, 0.85),
    run = c(1, 1, 1, 1) # Add run column needed for block mapping
  )

  # Create a temporary directory for the BIDS stub
  # Using tempdir() for broader compatibility than local_tempdir()
  temp_bids_dir <- tempfile(pattern = "fmridsl_bids_test_")
  dir.create(temp_bids_dir, recursive = TRUE)
  # Ensure cleanup happens when this test_that block finishes
  on.exit(unlink(temp_bids_dir, recursive = TRUE, force = TRUE), add = TRUE)

  # Create the mock BIDS structure on disk
  mock_proj_stub <- bidser::create_mock_bids(
    project_name = "MockTest",
    participants = participants_df,
    file_structure = file_structure_df,
    event_data = event_data_list,
    create_stub = TRUE,
    stub_path = temp_bids_dir
  )

  # Check if stub creation seemed successful (optional intermediate check)
  expect_true(fs::dir_exists(temp_bids_dir))
  expect_true(fs::file_exists(file.path(temp_bids_dir, "participants.tsv")))
  expect_true(fs::file_exists(file.path(temp_bids_dir, event_path_1)))

  # 2. Define Minimal YAML Configuration
  # -------------------------------------------------
  # Use sprintf to insert the temporary BIDS path
  # Ensure path separators are forward slashes for compatibility
  bids_path_for_yaml <- gsub("\\\\", "/", temp_bids_dir)

  yaml_content <- sprintf('
dataset:
  path: "%s"        # Path to the mock BIDS dataset
  subjects:
    include: [sub-01] # Specify the subject
  tasks: [taskA]      # Specify the task

events:
  onset_column: onset       # Column name for onsets in events.tsv
  duration_column: duration # Column name for durations
  block_column: run         # Column name for run identifier

# Define HRFs (optional if using defaults, but good practice)
hrfs:
  canonical:
    type: HRF_SPMG1      # Use SPM canonical HRF

# Define Terms (link model variables to event columns and HRF)
terms:
  stim:                   # Name of the model term
    type: hrf             # Standard HRF convolution
    variables: [condition]# Model variable name
    hrf: canonical        # Link to the HRF defined above

# Define the Model
model:
  name: "SimpleTestModel"
  variable_mapping:
    condition: trial_type # Map model variable 'condition' to BIDS column 'trial_type'
    run: run              # Map model variable 'run' to BIDS column 'run'
  terms:
    - stim                # Include the 'stim' term in this model
  variable_roles:         # Explicitly declare roles (good practice)
    factors: [condition]
  baseline:               # Minimal baseline
    basis: constant
    intercept: runwise
', bids_path_for_yaml)

  # Write YAML to a temporary file
  temp_yaml_file <- tempfile(fileext = ".yaml")
  writeLines(yaml_content, temp_yaml_file)
  # Ensure cleanup
  on.exit(unlink(temp_yaml_file), add = TRUE)

  # 3. Load and Validate Configuration
  # -------------------------------------------------
  config <- NULL
  expect_no_error({
    config <- load_fmri_config(temp_yaml_file)
  })

  # Basic checks on the loaded config object
  expect_s3_class(config, "fmri_config")
  expect_true(config$validated)
  expect_equal(config$subjects, "sub-01") # Check subject selection
  expect_equal(config$tasks, "taskA")     # Check task selection
  expect_equal(config$spec$model$name, "SimpleTestModel")
  expect_true("stim" %in% config$spec$model$terms)
  expect_true("condition" %in% names(config$spec$model$variable_mapping))
  expect_equal(config$spec$model$variable_mapping$condition, "trial_type")
  expect_true("condition" %in% config$variable_roles$factors) # Check inferred/declared role

  # 4. Build fmri_model from Config
  # -------------------------------------------------
  fmri_model_obj <- NULL
  expect_no_error({
    # Need to select the specific subject ID as stored in config$subjects
    # Note: build_fmri_model_from_config expects the subject ID *without* "sub-" prefix
    subject_id_for_build <- gsub("^sub-", "", config$subjects[1])
    fmri_model_obj <- build_fmri_model_from_config(config, subject_id_for_build)
  })

  # Basic checks on the built model object
  expect_s3_class(fmri_model_obj, "fmri_model")
  expect_s3_class(fmri_model_obj$event_model, "event_model")
  expect_s3_class(fmri_model_obj$baseline_model, "baseline_model")

  # Check terms in the built event model
  event_terms_built <- terms(fmri_model_obj$event_model)
  expect_true("stim" %in% names(event_terms_built))
  # Check conditions generated by the term
  expect_true(all(c("condition#condA", "condition#condB") %in% conditions(event_terms_built$stim)))

  # Check baseline model components
  baseline_terms_built <- terms(fmri_model_obj$baseline_model)
  # Should have 'constant' (run-wise intercept) and possibly drift if not constant
  expect_true("constant" %in% names(baseline_terms_built))
  # Check number of baseline columns (1 intercept per run + drift terms)
  # For 'constant' basis and 'runwise' intercept, expect 1 column per run.
  # Here, 1 run -> 1 column.
  expect_equal(ncol(design_matrix(fmri_model_obj$baseline_model)), 1)

  # Check dimensions of the full design matrix (optional but good)
  # Needs run_length and TR from the mock data setup
  # Mock data: n_volumes = 50, run_length = c(50), TR = 2.0 (from dataset_description)
  # This part requires adapting the setup_temp_bids or knowing the NIfTI details
  # Let's calculate expected rows based on the mock `event_data_list` and `TR`
  # Max onset = 15, duration = 0.5. Let's assume a reasonable HRF span (e.g., 32s)
  # Max time needed approx 15 + 32 = 47s. Number of TRs = ceil(47 / 2) = 24.
  # However, fmri_model uses the sampling_frame. Let's build it based on the mock file_structure.
  # The `create_mock_bids` doesn't easily provide run_length. We need to enhance mock setup
  # or make assumptions. Let's *assume* run length from the builder's data loading step.
  # The builder uses load_and_prepare_subject_data, which *should* get run_lengths.
  # Let's inspect the built model's sampling frame.
  built_sframe <- fmri_model_obj$event_model$sampling_frame
  expect_equal(sum(built_sframe$blocklens), 25) # Based on default 50 volumes and TR 2.0 in setup

  dmat <- design_matrix(fmri_model_obj)
  # Expected columns: 2 event conditions + 1 baseline intercept
  expected_cols <- 2 + 1
  expect_equal(ncol(dmat), expected_cols)
  expect_equal(nrow(dmat), sum(built_sframe$blocklens))

})