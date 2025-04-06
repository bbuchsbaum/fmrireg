# tests/testthat/helper-bids.R (or place directly in test file)

# Helper to create a temporary BIDS structure for testing
# Returns the root path of the temporary directory
setup_temp_bids <- function(
    subjects = c("sub-01"),
    tasks = c("task-test"),
    runs_per_task = 1,
    n_volumes = 50,
    event_cols = c("onset", "duration", "trial_type", "response_time"),
    confound_cols = c("trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z", "csf", "white_matter", "motion_outlier"),
    session = NULL, # Optional session name (e.g., "ses-01")
    tr = 2.0
) {

  # Create a unique temporary directory for this test run
  # Using testthat's local temp dir is safer for cleanup
  # Requires testthat >= 3.0.0 for local_tempdir
  # If older, use tempdir() and manage cleanup manually with on.exit
  # root_dir <- local_tempdir(pattern = "bids_test_") # Needs testthat >= 3.0.0
  # Manual temp dir for broader compatibility:
  root_dir <- tempfile(pattern = "bids_test_")
  dir.create(root_dir, recursive = TRUE)
  # Ensure cleanup happens when the calling test context exits
  # This requires the helper to be called *within* a test_that block or context
  # Alternatively, return the path and have the test call unlink.
  # Let's return the path and rely on the test for cleanup.

  # 1. Create dataset_description.json
  desc_content <- sprintf('{
    "Name": "Test BIDS Dataset",
    "BIDSVersion": "1.6.0",
    "Authors": ["Test Generator"],
    "RepetitionTime": %f
  }', tr)
  writeLines(desc_content, file.path(root_dir, "dataset_description.json"))

  # 2. Create participants.tsv
  part_df <- data.frame(participant_id = subjects)
  # Add dummy age column if needed
  if (nrow(part_df) > 0) part_df$age <- sample(20:40, nrow(part_df), replace = TRUE)
  write.table(part_df, file.path(root_dir, "participants.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # 3. Create subject/session/func directories and files
  for (sub in subjects) {
    sub_path <- file.path(root_dir, sub)
    session_path <- if (!is.null(session)) file.path(sub_path, session) else sub_path
    func_path <- file.path(session_path, "func")
    dir.create(func_path, recursive = TRUE)

    for (task in tasks) {
      for (run in 1:runs_per_task) {
        run_str <- sprintf("run-%02d", run)
        base_name <- paste(sub,
                           if (!is.null(session)) session else NULL,
                           paste0("task-", task),
                           run_str,
                           "bold", sep = "_")
        base_name_conf <- paste(sub,
                                if (!is.null(session)) session else NULL,
                                paste0("task-", task),
                                run_str,
                                "desc-confounds_timeseries", sep = "_") # FMRIPREP style

        # Create dummy events.tsv
        n_events <- sample(5:15, 1) # Random number of events
        events_df <- data.frame(matrix(runif(n_events * (length(event_cols)-2)), n_events, length(event_cols)-2))
        names(events_df) <- setdiff(event_cols, c("onset", "duration"))
        # Ensure specific types if needed for tests
        if("trial_type" %in% names(events_df)) events_df$trial_type <- sample(c("A","B"), n_events, replace=TRUE)
        if("response_time" %in% names(events_df)) events_df$response_time <- rnorm(n_events, 1, 0.2)

        events_df$onset <- sort(runif(n_events, 0, n_volumes * tr - 10)) # Ensure onsets are sorted
        events_df$duration <- runif(n_events, 0.5, 3)
        events_df <- events_df[, event_cols] # Ensure correct order
        write.table(events_df, file.path(func_path, paste0(base_name, "_events.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

        # Create dummy confounds.tsv
        confounds_df <- data.frame(matrix(rnorm(n_volumes * length(confound_cols)), n_volumes, length(confound_cols)))
        names(confounds_df) <- confound_cols
        write.table(confounds_df, file.path(func_path, paste0(base_name_conf, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

        # Create dummy bold.nii.gz (just an empty file for existence check)
        # For real tests needing data, neuroim2 would be used here.
        file.create(file.path(func_path, paste0(base_name, ".nii.gz")))
      }
    }
  }
  return(root_dir)
}


test_that("load_fmri_config handles minimal valid configuration", {
  # Setup temporary BIDS structure
  bids_root <- setup_temp_bids(
    subjects = c("sub-01"),
    tasks = c("task-test"),
    runs_per_task = 1,
    event_cols = c("onset", "duration", "trial_type"),
    confound_cols = c("csf", "white_matter") # Only provide needed confounds
  )
  # Ensure cleanup of the temporary directory
  on.exit(unlink(bids_root, recursive = TRUE, force = TRUE), add = TRUE)

  # Define minimal YAML content
  yaml_content <- sprintf('
dataset:
  path: "%s" # Use the temporary path
  subjects:
    include: [sub-01]
  tasks: [task-test] # Use the task name from setup_temp_bids

events:
  onset_column: onset
  duration_column: duration
  block_column: run # Need a block column in events data or mapping

terms:
  stim_effect:
    type: hrf
    variables: [condition] # This is the *model* variable name

model:
  name: "MinimalModel"
  variable_mapping:
    condition: trial_type # Map model condition
    run: .run # Map model run
              # OR map to an actual column if it exists
  terms:
    - stim_effect
  variable_roles: # Optional, but good practice
    factors: [condition]
  baseline: # Explicitly minimal baseline
    basis: constant
    intercept: global
    confound_variables: [] # No confounds used
', gsub("\\\\", "/", bids_root)) # Ensure forward slashes for path

  # Write YAML to a temporary file
  yaml_file <- tempfile(fileext = ".yaml")
  writeLines(yaml_content, yaml_file)
  on.exit(unlink(yaml_file), add = TRUE)

  # Load the configuration
  # Mocking read_events/read_confounds minimally if needed by build_config_from_ior
  # Depending on bidser implementation details, mocking might still be required
  # if it doesn't handle the simple structure correctly.
  # Assuming build_config_from_ior uses bidser functions that work with the temp structure:
  config <- NULL
  expect_no_error(config <- load_fmri_config(yaml_file))

  # Assertions on the loaded config object
  expect_s3_class(config, "fmri_config")
  expect_true(config$validated)
  expect_equal(config$subjects, c("sub-01"))
  expect_equal(config$tasks, c("task-test"))
  expect_equal(config$spec$model$name, "MinimalModel")
  expect_true("stim_effect" %in% config$spec$model$terms)
  expect_true("condition" %in% config$variable_roles$factors)
  expect_equal(config$spec$model$baseline$basis, "constant")
  expect_equal(length(config$spec$model$baseline$confound_variables), 0)
  expect_equal(config$events_info$mapping$onset_column, "onset")
})

# source("helper-bids.R") # If helper is separate

test_that("load_fmri_config handles features like parametric, contrasts, confounds", {
  # Setup BIDS structure with relevant columns
  bids_root <- setup_temp_bids(
    subjects = c("sub-01", "sub-02"),
    tasks = c("task-mixed"),
    runs_per_task = 2,
    event_cols = c("onset", "duration", "trial_type", "rt", "accuracy", "run_id"), # Added rt, accuracy, run_id
    confound_cols = c("csf", "white_matter", "trans_x", "rot_y", "motion_outlier01", "cosine00")
  )
  on.exit(unlink(bids_root, recursive = TRUE, force = TRUE), add = TRUE)

  # Define YAML using various features
  yaml_content <- sprintf('
dataset:
  path: "%s"
  subjects:
    include: [sub-01] # Select only sub-01
  tasks: [task-mixed]
  runs: [run-01, run-02] # Select both runs

events:
  onset_column: onset
  duration_column: duration
  block_column: run_id # Map to the column in the generated events file

hrfs:
  canonical:
    type: HRF_SPMG1
  shifted:
    type: HRF_SPMG1
    lag: 1.5 # Apply a lag to this HRF

confounds:
  include: ["csf", "white_matter", "^trans.*", "^rot.*"] # Select some confounds
  exclude: ["motion_outlier.*"] # Exclude outliers

terms:
  stim_effect:
    type: hrf
    variables: [condition] # Model variable name
    hrf: canonical
    subset: "accuracy == 1" # Use only correct trials for this term
  rt_mod:
    type: parametric
    variables: [condition, rt] # Modulate condition by rt
    hrf: shifted # Use the lagged HRF
    transform: [center] # Center the modulator
    basis:
      type: Poly
      parameters: { degree: 1 } # Linear basis for RT effect
    subset: "accuracy == 1"
  motion:
    type: nuisance # Include motion without convolution
    variables: [mx, my, mz, rx, ry, rz] # Model variable names for motion

contrasts:
  condA_vs_condB:
    type: formula
    expression: "condition[A] - condition[B]" # Contrast model variables

model:
  name: "FeatureTestModel"
  baseline:
    basis: bspline
    degree: 3
    intercept: runwise
    confound_variables: [wm, csf] # Use model variable names for baseline confounds
  variable_mapping:
    condition: trial_type # BIDS column
    rt: response_time    # BIDS column
    accuracy: accuracy     # BIDS column (will be inferred as factor, but we override)
    run_id: run_id         # BIDS column for block
    # Confound mappings
    wm: white_matter
    csf: csf
    mx: trans_x
    my: trans_y
    mz: trans_z
    rx: rot_x
    ry: rot_y
    rz: rot_z
  variable_roles:
    factors: [condition] # trial_type is already factor-like
    parametric: [rt, accuracy] # Treat numeric accuracy as parametric here
  terms: # Select terms defined globally
    - stim_effect
    - rt_mod
    - motion # Include motion nuisance term
  contrasts: # Select contrasts defined globally
    - condA_vs_condB
', gsub("\\\\", "/", bids_root))

  # Write YAML file
  yaml_file <- tempfile(fileext = ".yaml")
  writeLines(yaml_content, yaml_file)
  on.exit(unlink(yaml_file), add = TRUE)

  # Load config (mocking might still be needed if bidser has issues with temp files)
  # Assuming build_config_from_ior works with the temp structure:
  config <- NULL
  expect_no_error(config <- load_fmri_config(yaml_file))

  # --- Assertions ---
  expect_s3_class(config, "fmri_config")
  expect_true(config$validated)
  expect_equal(config$subjects, "sub-01") # Check subject filtering
  expect_equal(config$tasks, "task-mixed")

  # Check confound selection
  expect_true(!is.null(config$confounds_info))
  # Note: The exact selected columns depend on the generated mock data column names
  # We expect csf, white_matter, trans_x, rot_y to be selected
  # We expect motion_outlier01 to be excluded
  expect_true(all(c("csf", "white_matter", "trans_x", "rot_y") %in% config$confounds_info$columns))
  expect_false("motion_outlier01" %in% config$confounds_info$columns)
  expect_false("cosine00" %in% config$confounds_info$columns) # Not included

  # Check baseline confounds used in the model
  expect_equal(config$spec$model$baseline$confound_variables, c("wm", "csf"))

  # Check variable roles
  expect_true("condition" %in% config$variable_roles$factors)
  expect_true("rt" %in% config$variable_roles$parametric)
  expect_true("accuracy" %in% config$variable_roles$parametric) # Check override

  # Check selected terms and contrasts
  expect_true(all(c("stim_effect", "rt_mod", "motion") %in% config$spec$model$terms))
  expect_true("condA_vs_condB" %in% config$spec$model$contrasts)

  # Check term details were parsed correctly (e.g., parametric basis)
  expect_equal(config$spec$terms$rt_mod$basis$type, "Poly")
  expect_equal(config$spec$terms$rt_mod$transform, list("center"))

  # Check HRF details
  expect_equal(config$spec$hrfs$shifted$lag, 1.5)
})