library(testthat)
library(bidser)
library(mockery)


# Helper function to create mock event data
create_mock_events <- function(
    n_events = 20,
    onset_interval = 10,
    duration = 2,
    conditions = list(
      attend = list(values = c("on", "off"), each = n_events / 2),
      stimulus = list(values = c("face", "house"), rep = n_events / 2)
    ),
    continuous = list(
      RT = list(min = 0.5, max = 1.5),
      accuracy = list(prob = 0.8, n = n_events)
    )
) {
  onsets <- seq(0, by = onset_interval, length.out = n_events)
  events <- data.frame(
    onset = onsets,
    duration = rep(duration, n_events),
    block = seq_len(n_events)
  )
  # Add conditions
  for (cond in names(conditions)) {
    cond_info <- conditions[[cond]]
    if (!is.null(cond_info$values)) {
      if (!is.null(cond_info$each)) {
        events[[cond]] <- rep(cond_info$values, each = cond_info$each)
      } else if (!is.null(cond_info$rep)) {
        events[[cond]] <- rep(cond_info$values, cond_info$rep)
      }
    }
  }
  # Add continuous variables
  for (var in names(continuous)) {
    var_info <- continuous[[var]]
    if (!is.null(var_info$min) && !is.null(var_info$max)) {
      events[[var]] <- runif(n_events, min = var_info$min, max = var_info$max)
    } else if (!is.null(var_info$prob) && !is.null(var_info$n)) {
      events[[var]] <- rbinom(n_events, size = var_info$n, prob = var_info$prob)
    }
  }
  events
}

# Helper function to create mock confound data
create_mock_confounds <- function(
    n_volumes = 200,
    variables = list(
      motion1 = list(mean = 0, sd = 1),
      motion2 = list(mean = 0, sd = 1)
    )
) {
  confounds <- lapply(variables, function(var) {
    rnorm(n_volumes, mean = var$mean, sd = var$sd)
  })
  as.data.frame(confounds)
}

# Create S3 methods for bids_project class
#' @export
participants.bids_project <- function(x, ...) {
  # Return a tibble with participant_id column
  tibble::tibble(participant_id = x$part_df$participant_id)
}

#' @export
tasks.bids_project <- function(x, ...) {
  x$config$tasks
}

#' @export
read_events.bids_project <- function(x, subid = ".*", task = ".*") {
  # Validate inputs
  pt <- participants(x)
  idx <- grep(subid, pt$participant_id)
  
  if (length(idx) == 0) {
    stop(paste("no matching participants for 'subid' regex: ", subid))
  }
  
  sids <- pt$participant_id[idx]
  
  taskset <- tasks(x)
  task.idx <- grep(task, taskset)
  
  if (length(task.idx) == 0) {
    stop(paste("no matching tasks for 'task' regex: ", task))
  }
  
  tasks <- taskset[task.idx]
  
  # Generate events for each task/subject/run combination
  lapply(tasks, function(t) {
    lapply(sids, function(sid) {
      # Generate events for each run
      run_events <- lapply(seq_len(x$config$runs), function(run) {
        events <- do.call(create_mock_events, x$config$event_spec)
        # Add required columns in correct order
        events$.task <- t
        events$.run <- run
        events$.subid <- sid
        events
      }) %>% dplyr::bind_rows()
      
      if (nrow(run_events) > 0) {
        run_events
      } else {
        NULL
      }
    }) %>% 
      dplyr::bind_rows() %>%
      dplyr::group_by(.task, .run, .subid) %>%
      tidyr::nest()
  }) %>% dplyr::bind_rows()
}

#' @export
print.bids_events <- function(x, ...) {
  cat("BIDS Events Data\n")
  cat("Participant:", x$participant_id, "\n")
  cat("Task:", x$task, "\n")
  cat("Runs:", length(x$data), "\n")
  cat("Events per run:", sapply(x$data, nrow), "\n")
}

#' @export
read_confounds.bids_project <- function(x, subid, task=NULL, cvars=NULL, nest=TRUE, ...) {
  # Validate subject exists
  pt <- participants(x)
  idx <- grep(subid, pt$participant_id)
  
  if (length(idx) == 0) {
    stop(paste("no matching participants found for regex: ", subid))
  }
  
  # Generate confounds for each run
  run_confounds <- lapply(seq_len(x$config$runs), function(run) {
    do.call(create_mock_confounds, 
            c(list(n_volumes = x$config$n_volumes), 
              x$config$confound_spec))
  })
  
  if (nest) {
    data.frame(
      participant_id = rep(subid, x$config$runs),
      run = seq_len(x$config$runs),
      session = rep(1, x$config$runs),
      data = I(run_confounds),
      stringsAsFactors = FALSE
    )
  } else {
    run_confounds
  }
}

#' @export
print.bids_project <- function(x, ...) {
  cat("BIDS Project:", x$name, "\n")
  cat("Path:", x$path, "\n")
  cat("Subjects:", length(x$config$subjects), "\n")
  cat("Tasks:", paste(x$config$tasks, collapse=", "), "\n")
  cat("Runs per task:", x$config$runs, "\n")
}

read_events.bids_project <- function(x, subid = ".*", task = ".*", ...) {
  # Get matching participant IDs.
  pt <- participants(x)
  sids <- grep(subid, pt$participant_id, value = TRUE)
  
  # Get matching tasks.
  tset <- tasks(x)
  tasks <- grep(task, tset, value = TRUE)
  
  # For each combination of task, subject, and run, generate events.
  all_events <- lapply(tasks, function(t) {
    lapply(sids, function(sid) {
      lapply(seq_len(x$config$runs), function(run) {
        ev <- do.call(create_mock_events, x$config$event_spec)
        ev$.task <- t
        ev$.run <- run
        ev$.subid <- sid
        ev
      }) %>% dplyr::bind_rows()
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
  
  # Now group and nest the events so that the result has four columns:
  # .task, .run, .subid, and data.
  nested <- all_events %>%
    dplyr::group_by(.task, .run, .subid) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.task, .run, .subid)
  
  nested
}

# Main mock BIDS project creator
create_mock_bids_project <- function(
    subjects = c("sub-01", "sub-02"),
    tasks = c("task-1"),
    runs = 2,
    n_volumes = 200,
    event_spec = list(
      n_events = 20,
      onset_interval = 10,
      duration = 2,
      conditions = list(
        attend = list(values = c("on", "off"), each = 10),
        stimulus = list(values = c("face", "house"), rep = 10)
      ),
      continuous = list(
        RT = list(min = 0.5, max = 1.5),
        accuracy = list(prob = 0.8, n = 1)
      )
    ),
    confound_spec = list(
      variables = list(
        motion1 = list(mean = 0, sd = 1),
        motion2 = list(mean = 0, sd = 1)
      )
    )
) {
  # Create mock project structure
  mock_proj <- structure(
    list(
      name = "mock_project",
      path = "/mock/path",
      has_fmriprep = TRUE,
      has_sessions = FALSE,
      part_df = data.frame(
        participant_id = subjects,
        stringsAsFactors = FALSE
      ),
      config = list(
        subjects = subjects,
        tasks = tasks,
        runs = runs,
        n_volumes = n_volumes,
        event_spec = event_spec,
        confound_spec = confound_spec
      )
    ),
    class = "bids_project"
  )
  
  # Add mock methods directly to the object
  mock_proj$tasks <- function(...) tasks
  mock_proj$participants <- function(...) data.frame(participant_id = subjects)
  # mock_proj$read_events <- function(subid, task, ...) {
  #   events <- do.call(create_mock_events, event_spec)
  #   events$.task <- task
  #   events$.run <- 1:runs
  #   events$.subid <- subid
  #   events
  # }
  mock_proj$read_events <- function(subid, task, ...) {
    # Get matching participant IDs
    pt <- participants(mock_proj)
    sid_idx <- grep(subid, pt$participant_id)
    sids <- pt$participant_id[sid_idx]
    
    # Get matching tasks
    tset <- tasks(mock_proj)
    task_idx <- grep(task, tset)
    tasks <- tset[task_idx]
    
    # For each combination of task, subject, and run, generate events
    all_events <- lapply(tasks, function(t) {
      lapply(sids, function(sid) {
        lapply(seq_len(mock_proj$config$runs), function(run) {
          ev <- do.call(create_mock_events, mock_proj$config$event_spec)
          ev$.task <- t
          ev$.run <- run
          ev$.subid <- sid
          ev
        }) %>% dplyr::bind_rows()
      }) %>% dplyr::bind_rows()
    }) %>% dplyr::bind_rows()
    
    # Group by .task, .run, and .subid and nest the events into a "data" column
    nested <- all_events %>%
      dplyr::group_by(.task, .run, .subid) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::group_by(.task, .run, .subid)
    
    nested
  }
  mock_proj$read_confounds <- function(subid, ...) {
    data.frame(
      participant_id = rep(subid, runs),
      run = 1:runs,
      session = rep(1, runs),
      data = I(lapply(1:runs, function(r) {
        do.call(create_mock_confounds, 
                c(list(n_volumes = n_volumes), 
                  confound_spec))
      })),
      stringsAsFactors = FALSE
    )
  }
  
  mock_proj
}

# Example usage:
test_that("mock project creation is flexible", {
  # Default mock project
  mock_proj <- create_mock_bids_project()
  events <- bidser::read_events(mock_proj, "sub-01", "task-1")
  # Check the number of rows in the resulting grouped tibble, should be 2 (for 2 runs)
  expect_equal(nrow(events), 2) 
  expect_equal(nrow(events$data[[1]]), 20)  # 20 events per run
  
  # Custom mock project
  custom_proj <- create_mock_bids_project(
    subjects = c("sub-01", "sub-02", "sub-03"),
    tasks = c("task-1", "task-2"),
    runs = 3,
    event_spec = list(
      n_events = 30,
      onset_interval = 8,
      duration = 1.5,
      conditions = list(
        condition = list(values = c("A", "B", "C"), each = 10)
      ),
      continuous = list(
        score = list(min = 0, max = 100)
      )
    ),
    confound_spec = list(
      variables = list(
        fd = list(mean = 0.2, sd = 0.1),
        dvars = list(mean = 1, sd = 0.5)
      )
    )
  )
  
  events <- bidser::read_events(custom_proj, "sub-01", "task-1")
  expect_true(inherits(events, "data.frame")) # read_events returns grouped tibble now
  expect_equal(nrow(events), 3) # 3 runs
  expect_equal(nrow(events$data[[1]]), 30)  # 30 events per run
  expect_true("condition" %in% names(events$data[[1]]))
  expect_true("score" %in% names(events$data[[1]]))
  
  confounds <- bidser::read_confounds(custom_proj, "sub-01")
  expect_true(all(c("fd", "dvars") %in% names(confounds$data[[1]])))
})

# Test full config loading
test_that("load_fmri_config works correctly with minimal spec", {
  # Use a real temp directory
  temp_bids_dir <- tempfile("mockbids_")
  dir.create(temp_bids_dir, recursive = TRUE)
  on.exit(unlink(temp_bids_dir, recursive = TRUE), add = TRUE)

  mock_proj <- create_mock_bids_project(subjects = c("sub-01", "sub-02"), tasks = "task-1")
  # Add a dummy participants.tsv required by bidser
  writeLines(c("participant_id", "sub-01", "sub-02"), file.path(temp_bids_dir, "participants.tsv"))
  
  yaml_content <- '
dataset:
  path: "PLACEHOLDER_PATH"
  subjects:
    include: ["sub-01", "sub-02"]
  tasks: ["task-1"]

events: # Top-level events section REQUIRED
  onset: onset
  duration: duration
  block: block
  stimulus: { column: stimulus } # Variable used in regressor

confounds:
  include: ["motion1", "motion2"]

regressors: # Top-level regressors section REQUIRED
  baseline:
    type: hrf
    variables: [stimulus]
    hrf: HRF_SPMG1

model:
  name: "test_model"
  events: # Model-level event mapping
    onset: onset # Using column name directly
    duration: duration
    block: block
    stimulus: stimulus # Map model variable stimulus to event column stimulus
  regressors: # Reference top-level regressor
    baseline:

# HRFs section is optional
'

  # Inject real path
  yaml_content <- sub("PLACEHOLDER_PATH", temp_bids_dir, yaml_content, fixed = TRUE)
  
  # Print the final YAML content for debugging if needed
  # message("Minimal YAML Content:\n", yaml_content)
  
  yaml_file <- tempfile(fileext = ".yaml")
  writeLines(yaml_content, yaml_file)
  on.exit(unlink(yaml_file), add = TRUE) # Ensure temp file cleanup
  
  # Mock bidser functions that interact with the filesystem/project
  local_mocked_bindings(
    # Mock bids_project to prevent actual filesystem scanning beyond basic checks
    # It still needs the path to exist for the initial check
    bids_project = function(path, ...) {
        # Basic check if path exists, return the mock_proj structure
        if (!fs::dir_exists(path)) stop("Mock BIDS path doesn't exist: ", path)
        # Add the real path to the mock object for print methods etc.
        mock_proj$path <- path
        mock_proj
    },
    # participants and tasks can use the functions within mock_proj
    participants = function(proj, ...) proj$participants(),
    tasks = function(proj, ...) proj$tasks(),
    # Mock read_events to return correctly structured data
    read_events = function(proj, subid, task, ...) {
        # Simulate read_events structure expected by build_config_from_ior
        # Returns a list containing a dataframe for the first subject/task
        list(data=list(create_mock_events()))
    },
    # Mock read_confounds similarly
    read_confounds = function(proj, subid, ...) { 
        list(data=list(data.frame(motion1 = rnorm(10), motion2 = rnorm(10), motion_outlier=rnorm(10))))
    },
    .package = "bidser"
  )
  
  config <- load_fmri_config(yaml_file)
  
  expect_s3_class(config, "fmri_config")
  expect_equal(config$subjects, c("sub-01", "sub-02"))
  expect_equal(config$tasks, "task-1")
  
  # Additional checks for confounds
  expect_true(!is.null(config$confounds_info))
  expect_equal(config$confounds_info$columns, c("motion1", "motion2"))
})

# Test full config loading with comprehensive YAML
test_that("load_fmri_config handles complete specification", {
  # Use a real temp directory
  temp_bids_dir <- tempfile("mockbids_")
  dir.create(temp_bids_dir, recursive = TRUE)
  on.exit(unlink(temp_bids_dir, recursive = TRUE), add = TRUE)

  mock_proj <- create_mock_bids_project(subjects = c("sub-01", "sub-02"), tasks = "task-1")
  # Add dummy participants.tsv
  writeLines(c("participant_id", "sub-01", "sub-02"), file.path(temp_bids_dir, "participants.tsv"))
  
  yaml_content <- '
dataset:
  path: "PLACEHOLDER_PATH"
  subjects:
    include: ["sub-01", "sub-02"]
  tasks: ["task-1"]
  scan_params:
    TR:
      default: 2.0
      overrides:
        - value: 1.5
          pattern: "task-1_run-02"
    run_length:
      default:
        "task-1": 200
      overrides:
        - value: 180
          pattern: "sub-02_task-1_run-02"

events: # Top level events section
  onset: onset
  duration: duration
  block: block
  attend: { column: attend }
  stimulus: { column: stimulus }
  RT: { column: RT }
  accuracy: { column: accuracy }

confounds:
  include: ["^motion[12]$"]
  exclude: ["motion_outlier"]

regressors: # Top level regressors section
  main_effect:
    type: hrf
    variables: [attend, stimulus] # Uses variables mapped in top-level events
    hrf: HRF_SPMG1 # Assuming default built-in HRF
  rt_mod:
    type: hrf_parametric
    variables: [RT] # Uses variable mapped in top-level events
    hrf: HRF_SPMG1
    basis:
      type: Poly
      parameters:
        degree: 1

model:
  name: "full_model"
  factors: [attend, stimulus] # Explicitly declare factors
  parametric: [RT, accuracy]  # Explicitly declare parametrics
  events: # Model-specific event mapping (can reference top-level)
    onset: onset
    duration: duration
    block: block
    attend: attend       
    stimulus: stimulus   
    RT: RT               
    accuracy: accuracy   
  regressors: # Reference regressors defined above
    main_effect: 
    rt_mod:
'

  # Inject real path
  yaml_content <- sub("PLACEHOLDER_PATH", temp_bids_dir, yaml_content, fixed = TRUE)
  
  # message("Complete YAML Content:\n", yaml_content)
  
  yaml_file <- tempfile(fileext = ".yaml")
  writeLines(yaml_content, yaml_file)
  on.exit(unlink(yaml_file), add = TRUE)
  
  # Set up mocked bindings
  local_mocked_bindings(
    bids_project = function(path, ...) {
        if (!fs::dir_exists(path)) stop("Mock BIDS path doesn't exist: ", path)
        mock_proj$path <- path
        mock_proj
    },
    participants = function(proj, ...) proj$participants(),
    tasks = function(proj, ...) proj$tasks(),
    read_events = function(proj, subid, task, ...) {
        list(data=list(create_mock_events()))
    },
    read_confounds = function(proj, subid, ...) {
        # Create mock confounds data matching potential include patterns
        confounds_df <- data.frame(
            motion1 = rnorm(10),
            motion2 = rnorm(10),
            motion_outlier = sample(0:1, 10, replace=TRUE),
            other_var = rnorm(10)
        )
        list(data = list(confounds_df))
    },
    .package = "bidser"
  )
  
  config <- load_fmri_config(yaml_file)
  
  expect_s3_class(config, "fmri_config")
  expect_equal(config$subjects, c("sub-01", "sub-02"))
  expect_equal(config$tasks, "task-1")
  
  # Additional checks for confounds
  expect_true(!is.null(config$confounds_info))
  expect_equal(config$confounds_info$columns, c("motion1", "motion2"))
})

test_that("mock read_events matches bidser structure exactly", {
  mock_proj <- create_mock_bids_project(
    subjects = c("sub-01", "sub-02"),
    tasks = c("task-1", "task-2"),
    runs = 2
  )
  
  events <- read_events(mock_proj, "sub-.*", "task-.*")
  
  # Check structure
  expect_true(inherits(events, "data.frame"))
  expect_equal(names(events), c(".task", ".run", ".subid", "data"))
  
  # Check grouping
  expect_true(dplyr::is_grouped_df(events))
  expect_equal(dplyr::group_vars(events), c(".task", ".run", ".subid"))
  
  # Check contents
  expect_equal(unique(events$.task), c("task-1", "task-2"))
  expect_equal(unique(events$.run), c(1, 2))
  expect_equal(unique(events$.subid), c("sub-01", "sub-02"))
  
  # Check nested data
  first_data <- events$data[[1]]
  expect_true(is.data.frame(first_data))
  expect_true("onset" %in% names(first_data))
  expect_true("duration" %in% names(first_data))
})

test_that("participants method works correctly", {
  mock_proj <- create_mock_bids_project()
  
  # Set up mocked binding for participants
  local_mocked_bindings(
    participants = function(proj) {
      tibble::tibble(participant_id = proj$config$subjects)
    },
    .package = "bidser"
  )
  
  # Test direct participants call
  parts <- bidser::participants(mock_proj)
  expect_true(is.data.frame(parts))
  expect_equal(names(parts), "participant_id")
  expect_equal(parts$participant_id, c("sub-01", "sub-02"))
})
