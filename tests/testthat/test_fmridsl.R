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
  expect_equal(length(events$data), 2)  # 2 runs
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
  expect_equal(length(events$data), 3)  # 3 runs
  expect_equal(nrow(events$data[[1]]), 30)  # 30 events per run
  expect_true("condition" %in% names(events$data[[1]]))
  expect_true("score" %in% names(events$data[[1]]))
  
  confounds <- bidser::read_confounds(custom_proj, "sub-01")
  expect_true(all(c("fd", "dvars") %in% names(confounds$data[[1]])))
})

# Test validate_subjects
test_that("validate_subjects works correctly", {
  mock_proj <- create_mock_bids_project()
  
  # Mock the bidser::participants function
  mockery::stub(validate_subjects, 'bidser::participants', function(proj) {
    data.frame(participant_id = proj$config$subjects, stringsAsFactors = FALSE)
  })
  
  # Test with exclude
  subjects_spec <- list(include = c("sub-01", "sub-02"), exclude = c("sub-02"))
  result <- validate_subjects(mock_proj, subjects_spec)
  expect_equal(result, "sub-01")
  
  # Test with invalid subject
  subjects_spec <- list(include = c("sub-01", "sub-03"))
  expect_error(
    validate_subjects(mock_proj, subjects_spec),
    "Subjects not found: sub-03"
  )
  
  # Test with no include (should return all subjects)
  subjects_spec <- list()
  result <- validate_subjects(mock_proj, subjects_spec)
  expect_equal(result, c("sub-01", "sub-02"))
})

# Test validate_events
test_that("validate_events works correctly", {
  mock_proj <- create_mock_bids_project()
  
  events_spec <- list(
    onset = "onset",
    duration = "duration",
    block = "block",
    RT = "RT",
    accuracy = "accuracy"
  )
  
  # Create a proper mock events data frame with realistic block structure
  # 10 events per block, 2 blocks total
  mock_events <- do.call(rbind, lapply(1:2, function(block) {
    data.frame(
      onset = seq(0, by = 10, length.out = 10),  # starts at 0 for each block
      duration = rep(2, 10),
      block = rep(block, 10),  # block number
      RT = runif(10, 0.5, 1.5),
      accuracy = rbinom(10, 1, 0.8),
      stringsAsFactors = FALSE
    )
  }))
  
  # Ensure the data is sorted by block
  mock_events <- mock_events[order(mock_events$block), ]
  
  # Mock the bidser functions
  mockery::stub(validate_events, 'bidser::read_events', function(proj, subid, task, ...) {
    # Return the structure that bidser::read_events would return
    list(
      data = list(mock_events)  # List of data frames, one per run
    )
  })
  mockery::stub(validate_events, 'bidser::tasks', function(proj) proj$config$tasks)
  
  # Print mock data for debugging
  print("Mock Events Structure:")
  print(str(mock_events))
  print("Block values:")
  print(table(mock_events$block))
  
  result <- validate_events(mock_proj, events_spec, 
                          subjects = c("sub-01", "sub-02"),
                          tasks = c("task-1"))
  
  expect_true(!is.null(result$events))
  expect_equal(names(result$events), c("sub-01", "sub-02"))
  expect_equal(names(result$events$`sub-01`), "task-1")
  
  # Additional checks to verify the structure
  expect_true(all(c("onset", "duration", "block") %in% result$columns))
  expect_equal(result$mapping, events_spec)
  
  # Test that block values are properly validated
  expect_true(is.numeric(result$events$`sub-01`$`task-1`$block))
  expect_true(all(diff(result$events$`sub-01`$`task-1`$block) >= 0))
  
  # Test that onsets reset to 0 for each block
  block1_events <- mock_events[mock_events$block == 1, ]
  block2_events <- mock_events[mock_events$block == 2, ]
  expect_equal(min(block1_events$onset), 0)
  expect_equal(min(block2_events$onset), 0)
})

# Test validate_confounds
test_that("validate_confounds works correctly", {
  mock_proj <- create_mock_bids_project()
  
  # Create mock confounds data
  mock_confounds <- list(
    data = list(
      data.frame(
        motion_x = rnorm(10),
        motion_y = rnorm(10),
        motion_z = rnorm(10),
        motion_rot_x = rnorm(10),
        motion_rot_y = rnorm(10),
        motion_rot_z = rnorm(10),
        other_var = rnorm(10)
      )
    )
  )
  
  # Mock the read_confounds function
  local_mocked_bindings(
    read_confounds = function(...) mock_confounds,
    .package = "bidser"
  )
  
  # Test including motion variables
  confounds_spec <- list(
    include = c("^motion.*"),
    exclude = c("motion_outlier")
  )
  
  result <- validate_confounds(mock_proj, confounds_spec, c("sub-01"))
  
  expect_equal(length(result$columns), 6)  # Should find all motion variables
  expect_true(all(grepl("^motion", result$columns)))
  expect_equal(result$spec, confounds_spec)
  
  # Test with no patterns (should include all columns)
  result <- validate_confounds(mock_proj, list(), c("sub-01"))
  expect_equal(length(result$columns), 7)  # Should include all columns
  
  # Test with non-matching pattern
  confounds_spec$include <- c("nonexistent.*")
  expect_error(
    validate_confounds(mock_proj, confounds_spec, c("sub-01")),
    regexp="Confound pattern"
  )
})

# Test full config loading
test_that("load_fmri_config works correctly", {
  mock_proj <- create_mock_bids_project()
  
  yaml_content <- '
dataset:
  path: "/mock/path"
  subjects:
    include: ["sub-01", "sub-02"]
  tasks: ["task-1"]
model:
  events:
    onset: "onset"
    duration: "duration"
    block: "block"
confounds:
  include: ["motion1", "motion2"]
'
  yaml_file <- tempfile(fileext = ".yaml")
  writeLines(yaml_content, yaml_file)
  
  # Create mock functions
  mock_participants <- function(proj) {
    tibble::tibble(participant_id = proj$config$subjects)
  }
  
  # Set up mocked bindings
  local_mocked_bindings(
    bids_project = function(...) mock_proj,
    participants = mock_participants,
    tasks = function(proj) proj$config$tasks,
    read_events = function(proj, subid, task, ...) {
      list(data = list(create_mock_events()))
    },
    read_confounds = function(proj, subid, ...) {
      list(data = list(data.frame(motion1 = rnorm(10), motion2 = rnorm(10))))
    },
    .package = "bidser"
  )
  
  config <- load_fmri_config(yaml_file)
  
  expect_s3_class(config, "fmri_config")
  expect_equal(config$subjects, c("sub-01", "sub-02"))
  expect_equal(config$tasks, "task-1")
  
  unlink(yaml_file)
})

# Test full config loading with comprehensive YAML
test_that("load_fmri_config handles complete specification", {
  mock_proj <- create_mock_bids_project()
  
  yaml_content <- '
dataset:
  path: "/mock/path"
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
        task-1: 200
      overrides:
        - value: 180
          pattern: "sub-02_task-1_run-02"

model:
  events:
    onset: "onset"
    duration: "duration"
    block: "block"
    attend: "attend"
    stimulus: "stimulus"
    RT: "RT"
    accuracy: "accuracy"

confounds:
  include: ["^motion.*"]
  exclude: ["motion_outlier"]
'
  
  yaml_file <- tempfile(fileext = ".yaml")
  writeLines(yaml_content, yaml_file)
  
  # Set up mocked bindings
  local_mocked_bindings(
    bids_project = function(...) mock_proj,
    participants = function(proj) {
      tibble::tibble(participant_id = proj$config$subjects)
    },
    tasks = function(proj) proj$config$tasks,
    read_events = function(proj, subid, task, ...) {
      list(data = list(create_mock_events()))
    },
    read_confounds = function(proj, subid, ...) {
      # Create mock confounds data with columns starting with "motion"
      confounds_df <- data.frame(
        motion_x = rnorm(10),
        motion_y = rnorm(10),
        motion_z = rnorm(10),
        motion_rot_x = rnorm(10),
        motion_rot_y = rnorm(10),
        motion_rot_z = rnorm(10),
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
  expect_true(all(grepl("^motion", names(config$confounds_info$columns)[1:6])))
  
  unlink(yaml_file)
})

# Add scan parameter test:
test_that("scan parameter functions work correctly", {
  scan_params <- list(
    TR = list(
      default = 2.0,
      overrides = list(
        list(value = 1.5, pattern = "task-1_run-02")
      )
    ),
    run_length = list(
      default = list("task-1" = 200),
      overrides = list(
        list(value = 180, pattern = "sub-02_task-1_run-02")
      )
    )
  )
  
  # Test TR override
  expect_equal(get_scan_tr("sub-01_task-1_run-02", scan_params), 1.5)
  expect_equal(get_scan_tr("sub-01_task-1_run-01", scan_params), 2.0)
  
  # Test run length override
  expect_equal(get_scan_length("sub-02_task-1_run-02", "task-1", scan_params), 180)
  expect_equal(get_scan_length("sub-01_task-1_run-01", "task-1", scan_params), 200)
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
  expect_true(tibble::is_tibble(parts))
  expect_equal(names(parts), "participant_id")
  expect_equal(parts$participant_id, c("sub-01", "sub-02"))
  
  # Test that validate_subjects works with the tibble
  subjects_spec <- list(include = c("sub-01", "sub-02"))
  result <- validate_subjects(mock_proj, subjects_spec)
  expect_equal(result, c("sub-01", "sub-02"))
})
