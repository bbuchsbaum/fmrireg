#' Load and validate fMRI analysis configuration
#' 
#' @param yaml_file Path to YAML configuration file
#' @return An fmri_config object containing validated configuration and components
#' @export
load_fmri_config <- function(yaml_file) {
  # Read YAML spec
  spec <- yaml::read_yaml(yaml_file)
  
  # Load and validate BIDS project
  proj <- bidser::bids_project(spec$dataset$path, fmriprep = TRUE)
  
  # Validate subjects
  subjects <- validate_subjects(proj, spec$dataset$subjects)
  
  # Validate tasks
  tasks <- validate_tasks(proj, spec$dataset$tasks)
  
  # Validate events
  events_info <- validate_events(proj, spec$model$events, subjects, tasks)
  
  # Validate confounds if specified
  confounds_info <- NULL
  if (!is.null(spec$confounds)) {
    confounds_info <- validate_confounds(proj, spec$confounds, subjects)
  }
  
  # Create validated config object
  config <- structure(
    list(
      spec = spec,
      project = proj,
      subjects = subjects,
      tasks = tasks,
      events_info = events_info,
      confounds_info = confounds_info,
      validated = TRUE
    ),
    class = "fmri_config"
  )
  
  config
}

# Validation functions

#' @keywords internal
validate_subjects <- function(proj, subjects_spec) {
  available <- bidser::participants(proj)
  # Extract participant_id column from tibble
  available_ids <- available$participant_id
  
  if (!is.null(subjects_spec$include)) {
    missing <- setdiff(subjects_spec$include, available_ids)
    if (length(missing) > 0) {
      stop("Subjects not found: ", paste(missing, collapse=", "))
    }
    included <- subjects_spec$include
  } else {
    included <- available_ids
  }
  
  if (!is.null(subjects_spec$exclude)) {
    included <- setdiff(included, subjects_spec$exclude)
  }
  
  included
}

#' @keywords internal
validate_tasks <- function(proj, tasks_spec) {
  if (is.null(tasks_spec)) return(NULL)
  
  available <- bidser::tasks(proj)
  missing <- setdiff(tasks_spec, available)
  if (length(missing) > 0) {
    stop("Tasks not found: ", paste(missing, collapse=", "))
  }
  
  tasks_spec
}

#' @keywords internal
validate_events <- function(proj, events_spec, subjects, tasks) {
  # Get example events file for each subject to validate
  events_list <- lapply(subjects, function(subj) {
    task_events <- lapply(tasks, function(task) {
      events <- bidser::read_events(proj, 
                                  subid = subj,
                                  task = task)
      if (length(events) == 0) {
        stop(sprintf("No event files found for subject %s, task %s", subj, task))
      }
      events$data[[1]]
    })
    names(task_events) <- tasks
    task_events
  })
  names(events_list) <- subjects
  
  # Validate required columns exist in all event files
  for (subj in subjects) {
    for (task in tasks) {
      events_data <- events_list[[subj]][[task]]
      available_cols <- names(events_data)
      
      required_cols <- c(events_spec$onset, 
                        events_spec$duration,
                        events_spec$block)
      missing_cols <- setdiff(required_cols, available_cols)
      if (length(missing_cols) > 0) {
        stop(sprintf("Required columns not found in events file for subject %s, task %s: %s", 
                    subj, task, paste(missing_cols, collapse=", ")))
      }
      
      # Validate block values are numeric and increasing within each run
      block_vals <- events_data[[events_spec$block]]
      if (!is.numeric(block_vals) && !is.integer(block_vals)) {
        if (is.factor(block_vals)) {
          block_vals <- as.numeric(as.character(block_vals))
        } else {
          stop(sprintf("Block values must be numeric for subject %s, task %s", subj, task))
        }
      }
      
      # Check blocks are non-decreasing
      if (is.unsorted(block_vals)) {
        browser()
        stop(sprintf("Block values must be non-decreasing for subject %s, task %s", subj, task))
      }
    }
  }
  
  # Check other variables (needed for regressors)
  other_vars <- setdiff(names(events_spec), 
                       c("onset", "duration", "block"))
  for (var in other_vars) {
    if (!events_spec[[var]] %in% available_cols) {
      stop("Column '", events_spec[[var]], "' not found in events file")
    }
  }
  
  # Return validated events info and the loaded event data
  list(
    columns = available_cols,
    mapping = events_spec,
    events = events_list  # Store validated event data for later use
  )
}

#' @keywords internal
validate_confounds <- function(proj, confounds_spec, subjects) {
  # Get example confounds file to validate patterns
  example_confounds <- bidser::read_confounds(proj, subjects[1])
  if (length(example_confounds) == 0) {
    stop("No confound files found for subject ", subjects[1])
  }
  
  available_cols <- names(example_confounds$data[[1]])
  
  # Process include patterns
  included_cols <- character()
  if (!is.null(confounds_spec$include)) {
    for (pattern in confounds_spec$include) {
      matches <- grep(pattern, available_cols, value = TRUE)
      if (length(matches) == 0) {
        stop("Confound patterns match no variables: ", pattern)
      }
      included_cols <- c(included_cols, matches)
    }
  } else {
    included_cols <- available_cols
  }
  
  # Process exclude patterns
  if (!is.null(confounds_spec$exclude)) {
    for (pattern in confounds_spec$exclude) {
      matches <- grep(pattern, included_cols, value = TRUE)
      included_cols <- setdiff(included_cols, matches)
    }
  }
  
  # Return validated info
  list(
    spec = confounds_spec,
    columns = included_cols
  )
}

#' Get TR for a specific scan
#' @param scan_name Full name of scan (e.g., "sub-001_task-memory_run-01")
#' @param scan_params Scan parameters from config
#' @return TR value for this scan
#' @keywords internal
get_scan_tr <- function(scan_name, scan_params) {
  # Start with default
  tr <- scan_params$TR$default
  
  # Check overrides
  if (!is.null(scan_params$TR$overrides)) {
    for (override in scan_params$TR$overrides) {
      if (grepl(override$pattern, scan_name)) {
        tr <- override$value
        break
      }
    }
  }
  
  tr
}

#' Get run length for a specific scan
#' @param scan_name Full name of scan
#' @param task_name Task name
#' @param scan_params Scan parameters from config
#' @return Run length for this scan
#' @keywords internal
get_scan_length <- function(scan_name, task_name, scan_params) {
  # Start with default for this task
  length <- scan_params$run_length$default[[task_name]]
  
  # Check overrides
  if (!is.null(scan_params$run_length$overrides)) {
    for (override in scan_params$run_length$overrides) {
      if (grepl(override$pattern, scan_name)) {
        length <- override$value
        break
      }
    }
  }
  
  length
}

#' @export
print.fmri_config <- function(x, ...) {
  # Helper function to format lists
  format_list <- function(lst, indent = "  ") {
    if (length(lst) == 0) return("none")
    paste(paste0(indent, "- ", lst), collapse = "\n")
  }
  
  # Helper for section headers
  section_header <- function(title) {
    paste0("\n", crayon::bold(crayon::blue("=== ", title, " ===")))
  }
  
  # Helper for subsection headers
  subsection_header <- function(title) {
    paste0(crayon::bold(crayon::cyan("--- ", title, " ---")))
  }
  
  cat(crayon::bold(crayon::green("fMRI Analysis Configuration\n")))
  
  # Dataset Information
  cat(section_header("Dataset"))
  cat("\nProject Path:", x$project$path)
  cat("\nSubjects:", length(x$subjects), "total")
  cat("\n", format_list(x$subjects))
  
  # Tasks
  cat(section_header("Tasks"))
  if (length(x$tasks) > 0) {
    cat("\n", format_list(x$tasks))
  } else {
    cat("\nNo tasks specified (using all available)")
  }
  
  # Events Information
  cat(section_header("Events"))
  cat("\n", subsection_header("Mappings"))
  for (mapping in names(x$events_info$mapping)) {
    cat(sprintf("\n  %s -> %s", mapping, x$events_info$mapping[[mapping]]))
  }
  
  # Confounds (if present)
  if (!is.null(x$confounds_info)) {
    cat(section_header("Confounds"))
    if (!is.null(x$confounds_info$spec$include)) {
      cat("\n", subsection_header("Include Patterns"))
      cat("\n", format_list(x$confounds_info$spec$include))
    }
    if (!is.null(x$confounds_info$spec$exclude)) {
      cat("\n", subsection_header("Exclude Patterns"))
      cat("\n", format_list(x$confounds_info$spec$exclude))
    }
  }
  
  # Model Specification
  if (!is.null(x$spec$model)) {
    cat(section_header("Model"))
    
    # Factors
    if (!is.null(x$spec$model$factors)) {
      cat("\n", subsection_header("Factors"))
      cat("\n", format_list(x$spec$model$factors))
    }
    
    # Parametric Modulators
    if (!is.null(x$spec$model$parametric)) {
      cat("\n", subsection_header("Parametric Modulators"))
      cat("\n", format_list(x$spec$model$parametric))
    }
    
    # HRFs
    if (!is.null(x$spec$hrfs)) {
      cat("\n", subsection_header("HRF Specifications"))
      for (hrf_name in names(x$spec$hrfs)) {
        hrf <- x$spec$hrfs[[hrf_name]]
        cat(sprintf("\n  %s: %s", hrf_name, hrf$type))
        if (!is.null(hrf$parameters)) {
          params <- paste(names(hrf$parameters), hrf$parameters, sep = "=", collapse = ", ")
          cat(sprintf(" (%s)", params))
        }
      }
    }
    
    # Contrasts
    if (!is.null(x$spec$model$contrasts)) {
      cat("\n", subsection_header("Contrasts"))
      for (con_name in names(x$spec$model$contrasts)) {
        con <- x$spec$model$contrasts[[con_name]]
        cat(sprintf("\n  %s: %s", con_name, con$type))
        if (!is.null(con$expression)) {
          cat(sprintf(" [%s]", con$expression))
        }
      }
    }
  }
  
  # Validation Status
  cat(section_header("Validation"))
  cat("\nStatus:", if(x$validated) crayon::green("✓ Validated") else crayon::red("✗ Not Validated"))
  cat("\n")
}
