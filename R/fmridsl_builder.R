
# --- Helper Functions ---

#' @keywords internal
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

# get_hrf_from_config (Assume defined as before)
# apply_transformations (Assume defined as before)
# create_parametric_basis_object (Assume defined as before)

# --- Function to Load and Prepare Subject Data ---

#' Load Events and Confounds for a Subject
#'
#' @param config Validated `fmri_config` object.
#' @param subject_id Subject ID string.
#' @return A list containing:
#'   - `events`: Combined data frame of events across runs for the subject.
#'   - `confounds`: Combined data frame of selected confounds across runs.
#'   - `run_lengths`: Vector of scan counts per run.
#'   - `TR`: Effective TR for the subject.
#' @keywords internal
load_and_prepare_subject_data <- function(config, subject_id) {

  tasks_to_load <- config$tasks
  runs_to_load <- config$runs # May be "all" or specific run IDs

  # --- Load Events ---
  all_subject_events <- list()
  loaded_run_lengths <- list()
  effective_TR <- config$spec$dataset$scan_params$TR$default %||% NA_real_

  message(sprintf("Loading data for subject %s, tasks: %s", subject_id, paste(tasks_to_load, collapse=", ")))

  tryCatch({
      events_read <- bidser::read_events(config$project, subid = subject_id, task = paste(tasks_to_load, collapse="|"))

      # Filter by run if necessary
      if (!is.null(runs_to_load) && !identical(runs_to_load, "all")) {
          run_nums_to_keep <- as.numeric(gsub("run-", "", runs_to_load))
          events_read <- events_read %>% dplyr::filter(as.numeric(.run) %in% run_nums_to_keep)
      }

      if (nrow(events_read) == 0) stop("No matching event files found/read.")

      # Unnest and combine
      all_subject_events <- events_read %>%
          tidyr::unnest(data) %>%
          dplyr::arrange(.run) # Ensure ordered by run

      # --- Determine Run Lengths and TR ---
      # This part needs robust logic based on BIDS metadata or config overrides
      # Placeholder: Infer from scan_params config or assume constant
      default_run_lengths <- config$spec$dataset$scan_params$run_length$default %||% list()
      run_length_overrides <- config$spec$dataset$scan_params$run_length$overrides %||% list()
      tr_overrides <- config$spec$dataset$scan_params$TR$overrides %||% list()

      run_info <- all_subject_events %>% distinct(.run, .task) %>% arrange(as.numeric(.run))
      final_run_lengths <- rep(NA_integer_, nrow(run_info))

      for(i in 1:nrow(run_info)){
           run_id_str <- paste0("run-", sprintf("%02d", as.numeric(run_info$.run[i]))) # Create BIDS-like run string
           task_id_str <- run_info$.task[i]
           
           # Check overrides first (more specific)
           matched_override <- FALSE
           
           # Run Length Overrides
           for(ovr in run_length_overrides) {
                pattern <- ovr$pattern
                # Construct a string to match against, e.g., "sub-01_task-nback_run-01"
                test_str <- paste(subject_id, task_id_str, run_id_str, sep="_") 
                if (grepl(pattern, test_str)) {
                    final_run_lengths[i] <- ovr$value
                    matched_override <- TRUE
                    break
                }
           }
           
           # Use task default if no override matched
           if (!matched_override) {
               final_run_lengths[i] <- default_run_lengths[[task_id_str]] %||% NA_integer_
           }
           
           # TR Override (apply same logic, update effective_TR if needed - assuming 1 TR for now)
            # ... similar loop for tr_overrides ...
            # if (tr_override_found) effective_TR <- ovr$value
           
      }
      
      if(anyNA(final_run_lengths)) stop("Could not determine run lengths for all runs.")
      if(anyNA(effective_TR)) stop("Could not determine TR.")
      
      loaded_run_lengths <- final_run_lengths


  }, error = function(e) {
      stop(sprintf("Error loading/processing event data for subject %s: %s", subject_id, e$message))
  })


  # --- Load Confounds ---
  all_subject_confounds <- NULL
  if (length(config$confounds_info$columns) > 0) {
      tryCatch({
          confounds_read <- bidser::read_confounds(config$project, subid = subject_id) # Add task/run?

          # Filter by run if necessary (matching events_read filtering)
          if (!is.null(runs_to_load) && !identical(runs_to_load, "all")) {
              run_nums_to_keep <- as.numeric(gsub("run-", "", runs_to_load))
              confounds_read <- confounds_read %>% dplyr::filter(as.numeric(run) %in% run_nums_to_keep)
          }

          if (nrow(confounds_read) == 0) stop("No matching confound files found/read.")

          # Unnest, select columns, and combine
          all_subject_confounds <- confounds_read %>%
               tidyr::unnest(data) %>%
               dplyr::select(all_of(config$confounds_info$columns)) %>%
               dplyr::arrange(as.numeric(run)) # Ensure order matches events/run_lengths

           # Check if number of rows matches total run lengths
           if (nrow(all_subject_confounds) != sum(loaded_run_lengths)) {
               stop(sprintf("Total number of confound rows (%d) does not match total expected scans (%d).",
                     nrow(all_subject_confounds), sum(loaded_run_lengths)))
           }

      }, error = function(e) {
          stop(sprintf("Error loading/processing confound data for subject %s: %s", subject_id, e$message))
      })
  }


  return(list(
    events = all_subject_events,
    confounds = all_subject_confounds, # Might be NULL
    run_lengths = loaded_run_lengths,
    TR = effective_TR
  ))
}

#' Build fmri_model Object from fmri_config
#'
#' Constructs an `fmri_model` object for a specific subject based on the
#' validated `fmri_config`.
#'
#' @param config A validated `fmri_config` object.
#' @param subject_id The subject ID string.
#' @return An `fmri_model` object.
#' @export
build_fmri_model_from_config <- function(config, subject_id) {

  if (!inherits(config, "fmri_config") || !isTRUE(config$validated)) {
    stop("'config' must be a validated object.")
  }
  if (!subject_id %in% config$subjects) {
    stop(sprintf("Subject '%s' not in config.", subject_id))
  }

  message(sprintf("Building fmri_model for subject: %s, model: %s", subject_id, config$spec$model$name))

  # --- 1. Load Subject Data ---
  subject_data <- load_and_prepare_subject_data(config, subject_id)
  subject_events_df <- subject_data$events
  subject_confounds_df <- subject_data$confounds
  subject_run_lengths <- subject_data$run_lengths
  subject_TR <- subject_data$TR

  # --- 2. Prepare Mapped Data Frame ---
  model_data <- data.frame(row.names = seq_len(nrow(subject_events_df)))
  model_variable_map <- config$spec$model$variable_mapping
  model_factors <- config$variable_roles$factors
  model_parametric <- config$variable_roles$parametric

  for (model_var in names(model_variable_map)) {
      bids_col <- model_variable_map[[model_var]]
      if (bids_col %in% names(subject_events_df)) {
         col_data_raw <- subject_events_df[[bids_col]]
         # Apply role conversion
         if (model_var %in% model_factors && !is.factor(col_data_raw)) {
            model_data[[model_var]] <- factor(col_data_raw)
         } else if (model_var %in% model_parametric && !is.numeric(col_data_raw)) {
            # Attempt conversion for parametric, warn if issues
            col_data_num <- suppressWarnings(as.numeric(as.character(col_data_raw)))
            if (anyNA(col_data_num) && !all(is.na(col_data_raw))) { # Check if conversion introduced NAs where original wasn't NA
               warning(sprintf("Conversion to numeric for parametric variable '%s' (column '%s') produced NAs.", model_var, bids_col))
            }
            model_data[[model_var]] <- col_data_num
         } else {
            # Use raw data type
            model_data[[model_var]] <- col_data_raw
         }
      } else {
         # This variable must map to confounds (or error)
         # Confounds are handled in baseline, just ensure mapping exists
         if (!(bids_col %in% names(subject_confounds_df %||% data.frame()))) {
              stop(sprintf("Mapped BIDS column '%s' for model variable '%s' not found in subject's events or loaded confounds.", bids_col, model_var))
         }
         # We don't add confound columns to model_data here
      }
  }
  # Ensure essential columns are present in model_data for event_model call
  essential_event_cols <- config$spec$events[c("onset_column", "duration_column", "block_column")]
  for(col_name in essential_event_cols) {
      if (!col_name %in% names(model_data)) {
          if(col_name %in% names(subject_events_df)) {
              model_data[[col_name]] <- subject_events_df[[col_name]]
          } else {
              stop(sprintf("Essential event column '%s' is missing from event data.", col_name))
          }
      }
  }
  # Make sure block column is numeric
  block_col_mapped <- config$spec$events$block_column
  if (!is.numeric(model_data[[block_col_mapped]])) {
      model_data[[block_col_mapped]] <- as.numeric(as.character(model_data[[block_col_mapped]]))
      if(anyNA(model_data[[block_col_mapped]])) stop("Block column could not be converted to numeric.")
  }


  # --- 3. Build sampling_frame ---
  sampling_frame <- sampling_frame(blocklens = subject_run_lengths, TR = subject_TR)

  # --- 4. Build baseline_model ---
  baseline_spec <- config$spec$model$baseline
  baseline_confound_model_vars <- baseline_spec$confound_variables %||% list()
  baseline_nuisance_list <- NULL

  if (length(baseline_confound_model_vars) > 0) {
     if (is.null(subject_confounds_df)) {
         stop("Baseline model requires confounds, but none were loaded/selected.")
     }
     bids_conf_cols_for_baseline <- unlist(model_variable_map[baseline_confound_model_vars])
     confounds_for_baseline <- subject_confounds_df[, bids_conf_cols_for_baseline, drop = FALSE]

     # Split by run
     run_ids_for_confounds <- rep(1:length(subject_run_lengths), times=subject_run_lengths)
     baseline_nuisance_list <- split.data.frame(confounds_for_baseline, run_ids_for_confounds)
     baseline_nuisance_list <- lapply(baseline_nuisance_list, as.matrix)
  }

  baseline_model <- baseline_model(
    basis = baseline_spec$basis,
    degree = baseline_spec$degree,
    sframe = sampling_frame,
    intercept = baseline_spec$intercept,
    nuisance_list = baseline_nuisance_list
  )

  # --- 5. Build event_model ---
  event_model_terms_specs <- list() # To hold hrfspec objects
  event_term_names <- config$spec$model$terms %||% list()

  # Retrieve globally defined contrasts selected for this model
  selected_contrast_names <- config$spec$model$contrasts %||% list()
  model_contrasts <- config$spec$contrasts[selected_contrast_names] %||% list()

  for (term_name in event_term_names) {
      term_def <- config$spec$terms[[term_name]]
      if (is.null(term_def)) {
          warning(sprintf("Term '%s' skipped (not defined globally).", term_name)); next
      }

      hrf_name <- term_def$hrf %||% "canonical"
      hrf_object_base <- get_hrf_from_config(hrf_name, config)
      term_lag <- term_def$lag %||% (attr(hrf_object_base, "lag") %||% 0) # Use term lag > HRF lag > 0
      
      # Apply lag if needed
      hrf_object <- if (term_lag != 0) {
          gen_hrf_lagged(hrf_object_base, lag = term_lag)
      } else {
          hrf_object_base
      }

      subset_expr <- if (!is.null(term_def$subset)) parse(text = term_def$subset)[[1]] else NULL

      vars_for_spec <- list()
      term_data_env <- list2env(as.list(model_data), parent = baseenv()) # Env for this term

      if (term_def$type == "parametric") {
          if (length(term_def$variables) < 1) stop("Parametric term needs variables.")
          modulator_var_name <- tail(term_def$variables, 1)
          event_selector_vars <- head(term_def$variables, -1)

          if (!modulator_var_name %in% names(model_data)) stop(sprintf("Modulator '%s' not in model_data.", modulator_var_name))

          modulator_data_raw <- model_data[[modulator_var_name]]
          modulator_data_transformed <- apply_transformations(
              modulator_data_raw,
              term_def$transform,
              subject_ids = model_data[[config$spec$events$block_column]]
          )

          # Handle basis expansion
          if (!is.null(term_def$basis)) {
              basis_object <- create_parametric_basis_object(
                  term_def$basis,
                  modulator_data_transformed,
                  modulator_var_name
              )
              # Add selectors first, then the basis object
              for(sel_var in event_selector_vars) vars_for_spec[[length(vars_for_spec)+1]] <- rlang::sym(sel_var)
              vars_for_spec[[length(vars_for_spec)+1]] <- basis_object

          } else {
              # No basis: Use transformed data directly. Add selectors, then the transformed modulator.
              assign(modulator_var_name, modulator_data_transformed, envir = term_data_env) # Update env
              for(sel_var in event_selector_vars) vars_for_spec[[length(vars_for_spec)+1]] <- rlang::sym(sel_var)
              vars_for_spec[[length(vars_for_spec)+1]] <- rlang::sym(modulator_var_name) # Use symbol
          }

      } else { # hrf, trialwise, nuisance
           for (var_name in term_def$variables) {
                vars_for_spec[[length(vars_for_spec) + 1]] <- rlang::sym(var_name)
           }
      }

      # Construct hrfspec or nuisancespec based on type
      # Need to map contrasts correctly based on which ones apply to this term
      # For now, pass all model contrasts; filtering might be needed in contrast_weights
      term_type <- term_def$type %||% "hrf" # Default to hrf if type missing

      if (term_type == "nuisance") {
          # Create a nuisancespec or similar structure that construct.nuisancespec can handle
          # It should result in a matrix_term when constructed
          # Placeholder: Treat as hrfspec for now, but needs specific handling
          warning("Nuisance term type needs specific implementation in builder/construct.")
           term_spec_obj <- hrfspec(
             vars = vars_for_spec, basis = HRF_SPMG1, # Dummy HRF
             onsets = model_data[[config$spec$events$onset_column]],
             durations = model_data[[config$spec$events$duration_column]],
             subset = subset_expr, contrasts = model_contrasts,
             id = term_name, data_env = term_data_env
          )
          term_spec_obj$.term_type <- "nuisance" # Add flag
      } else {
           # For hrf, parametric, trialwise
           # Need to ensure trialwise generates the right 'vars' (e.g., trial index factor)
           if (term_type == "trialwise") {
               if (length(term_def$variables) != 1) stop("Trialwise term requires exactly one variable (trial identifier).")
               trial_id_var <- term_def$variables[1]
               if (!trial_id_var %in% names(model_data)) stop(sprintf("Trial ID variable '%s' not found.", trial_id_var))
               # Create a unique factor for trials
               unique_trial_factor <- factor(paste0("trial_", seq_len(nrow(model_data)))) # Safer than using potentially non-unique values
               assign(trial_id_var, unique_trial_factor, envir=term_data_env) # Use this factor
               vars_for_spec <- list(rlang::sym(trial_id_var)) # Reset vars to just this factor
           }

           term_spec_obj <- hrfspec(
               vars = vars_for_spec,
               basis = hrf_object,
               onsets = model_data[[config$spec$events$onset_column]],
               durations = model_data[[config$spec$events$duration_column]],
               subset = subset_expr,
               contrasts = model_contrasts, # Pass all model contrasts for now
               id = term_name,
               data_env = term_data_env,
               summate = term_def$hrf_summate %||% TRUE # Add other HRF params if needed
           )
           term_spec_obj$.term_type <- term_type # Store the intended type
      }

      event_model_terms_specs[[term_name]] <- term_spec_obj
  }

  # Create the event_model
  # Need to handle the term types correctly - create_event_model might need updates
  # or we use construct() directly on each spec and combine terms.
  # Using create_event_model for now, assuming it handles different specs.

  event_model <- create_event_model(
      event_terms = event_model_terms_specs, # Pass the list of specs
      events = model_data,
      onsets = model_data[[config$spec$events$onset_column]],
      block = model_data[[config$spec$events$block_column]], # Pass actual block data
      sampling_frame = sampling_frame,
      durations = model_data[[config$spec$events$duration_column]]
  )

  # --- 6. Combine into fmri_model ---
  final_model <- fmri_model(event_model, baseline_model)

  message(sprintf("Successfully built fmri_model for subject %s.", subject_id))
  return(final_model)
}

# --- Add/Ensure Required R Functions ---

# Transformation functions
center <- function(x) scale(x, center = TRUE, scale = FALSE)[, 1]
scale_sd <- function(x) scale(x, center = FALSE, scale = TRUE)[, 1] # Renamed to avoid conflict
zscore <- function(x) scale(x, center = TRUE, scale = TRUE)[, 1]
within_subject <- function(x, subject_ids) {
   if (length(x) != length(subject_ids)) stop("Length mismatch in within_subject.")
   ave(x, subject_ids, FUN = function(v) scale(v, center = TRUE, scale = FALSE)[,1])
}

# Basis function wrappers (if not already present/correctly implemented)
# Ensure Poly, BSpline, Standardized, Ident exist in fmrireg namespace
# and handle parameters as expected by create_parametric_basis_object.
# Example for Poly possibly handling intercept:
# Poly <- function(x, degree, intercept = TRUE) {
#    pres <- stats::poly(x, degree)
#    if (!intercept) {
#       # Logic to remove intercept contribution if needed, or rely on model fitting
#    }
#    # ... rest of Poly constructor ...
# }