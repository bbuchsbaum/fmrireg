#' Create HRF Object from Specification
#'
#' @param hrf_spec A list containing HRF specifications:
#'   - type: Type of HRF (e.g., "HRF_SPMG1", "HRF_GAUSSIAN")
#'   - parameters: Optional list of parameters specific to the HRF type
#'   - definition: Required for custom HRFs, name of function
#'
#' @return An HRF object
#' @examples
#' # Basic canonical HRF
#' create_hrf_object(list(type = "HRF_SPMG1"))
#' 
#' # Gaussian HRF with custom parameters
#' create_hrf_object(list(
#'   type = "HRF_GAUSSIAN",
#'   parameters = list(mean = 5, sd = 1)
#' ))
#' 
#' @export
create_hrf_object <- function(hrf_spec) {
  # Get all available HRF functions
  hrf_funs <- ls(envir = asNamespace("fmrireg"), pattern = "^HRF_")
  
  type <- hrf_spec$type
  if (!type %in% c(hrf_funs, "custom")) {
    stop("Unknown HRF type: ", type, 
         "\nAvailable types: ", paste(c(hrf_funs, "custom"), collapse=", "))
  }
  
  # Handle custom HRFs
  if (type == "custom") {
    if (is.null(hrf_spec$definition)) {
      stop("Custom HRF requires 'definition' field specifying function name")
    }
    
    fun <- try(get(hrf_spec$definition), silent = TRUE)
    if (inherits(fun, "try-error")) {
      stop("Custom HRF function '", hrf_spec$definition, "' not found")
    }
    
    return(do.call(gen_hrf, c(list(fun), hrf_spec$parameters)))
  }
  
  # Get HRF function
  hrf_fun <- get(type, envir = asNamespace("fmrireg"))
  
  # If no parameters, return the default HRF
  if (is.null(hrf_spec$parameters)) {
    return(hrf_fun)
  }
  
  # Validate parameters based on HRF type
  valid_params <- switch(type,
    "HRF_GAUSSIAN" = c("mean", "sd"),
    "HRF_GAMMA" = c("shape", "rate", "scale"),
    "HRF_SPMG1" = character(0),  # No parameters
    "HRF_SPMG2" = character(0),
    "HRF_SPMG3" = character(0),
    "HRF_BSPLINE" = c("knots", "degree"),
    stop("Unknown HRF type: ", type)
  )
  
  # Check for invalid parameters
  invalid_params <- setdiff(names(hrf_spec$parameters), valid_params)
  if (length(invalid_params) > 0) {
    stop("Invalid parameters for ", type, ": ", 
         paste(invalid_params, collapse=", "), 
         "\nValid parameters are: ", paste(valid_params, collapse=", "))
  }
  
  # Create HRF with parameters
  switch(type,
    "HRF_GAUSSIAN" = {
      params <- list(
        mean = hrf_spec$parameters$mean %||% 5,
        sd = hrf_spec$parameters$sd %||% 1
      )
      do.call(gen_hrf, c(list(hrf_gaussian), params))
    },
    
    "HRF_GAMMA" = {
      params <- list(
        shape = hrf_spec$parameters$shape %||% 6,
        rate = hrf_spec$parameters$rate %||% 1,
        scale = hrf_spec$parameters$scale %||% 1
      )
      do.call(gen_hrf, c(list(hrf_gamma), params))
    },
    
    "HRF_BSPLINE" = {
      params <- list(
        knots = hrf_spec$parameters$knots %||% seq(0, 32, by=4),
        degree = hrf_spec$parameters$degree %||% 3
      )
      do.call(gen_hrf, c(list(hrf_bspline), params))
    },
    
    # Return built-in HRFs directly
    hrf_fun
  )
}

#' Validate HRF parameters
#' 
#' @param type HRF type
#' @param params List of parameters
#' @return TRUE if valid, throws error otherwise
#' @keywords internal
validate_hrf_params <- function(type, params) {
  # Parameter validation rules
  rules <- list(
    HRF_GAUSSIAN = list(
      mean = function(x) is.numeric(x) && length(x) == 1,
      sd = function(x) is.numeric(x) && x > 0
    ),
    HRF_GAMMA = list(
      shape = function(x) is.numeric(x) && x > 0,
      rate = function(x) is.numeric(x) && x > 0,
      scale = function(x) is.numeric(x) && x > 0
    ),
    HRF_BSPLINE = list(
      knots = function(x) is.numeric(x) && is.vector(x) && length(x) >= 2,
      degree = function(x) is.numeric(x) && x > 0 && x == round(x)
    )
  )
  
  if (!type %in% names(rules)) {
    return(TRUE)  # No validation for types not listed
  }
  
  type_rules <- rules[[type]]
  
  # Check each parameter
  for (param_name in names(params)) {
    if (!param_name %in% names(type_rules)) {
      stop("Invalid parameter '", param_name, "' for ", type)
    }
    
    validate_fn <- type_rules[[param_name]]
    if (!validate_fn(params[[param_name]])) {
      stop("Invalid value for parameter '", param_name, "' in ", type)
    }
  }
  
  TRUE
}

#' Construct Event Term from Regressor Specification
#'
#' This function constructs an event term based on the regressor specification, events data, and HRF objects.
#'
#' @param regressor_spec A list containing the specification of a single regressor from the YAML configuration, including 'type', 'variables', 'hrf', and optionally 'subset'.
#' @param events_data A data frame containing event information, including onsets, durations, and any variables referenced in the regressor specification.
#' @param hrf_objects A named list of HRF objects, where names correspond to HRF identifiers used in 'regressor_spec$hrf'.
#'
#' @return A list containing the 'event_term' object, the 'hrf' object, and the 'type' of the regressor.
#'
#' @examples
#' # Assume that 'events_data' is a data frame with columns 'onset', 'duration', 'stim'
#' events_data <- data.frame(
#'   onset = c(0, 10, 20),
#'   duration = c(1, 1, 1),
#'   stim = c('A', 'B', 'A')
#' )
#' 
#' # Create HRF objects
#' hrf_spec <- list(type = "HRF_SPMG1")
#' hrf_obj <- create_hrf_object(hrf_spec)
#' hrf_objects <- list(canonical = hrf_obj)
#' 
#' # Regressor specification
#' regressor_spec <- list(
#'   type = "hrf",
#'   variables = c("stim"),
#'   hrf = "canonical",
#'   subset = NULL
#' )
#' 
#' # Construct the event term
#' term_info <- construct_event_term_from_spec(regressor_spec, events_data, hrf_objects)
#' 
#' # Access the event term
#' event_term_obj <- term_info$event_term
#' 
#' @export
construct_event_term_from_spec <- function(regressor_spec, events_data, hrf_objects) {
  regressor_type <- regressor_spec$type
  variables <- regressor_spec$variables
  hrf_name <- regressor_spec$hrf
  hrf_obj <- hrf_objects[[hrf_name]]
  
  # Handle subset if specified
  if (!is.null(regressor_spec$subset)) {
    subset_expr <- parse(text = regressor_spec$subset)
    subs <- eval(subset_expr, envir = events_data)
    events_data <- events_data[subs, , drop = FALSE]
  }
  
  onsets <- events_data$onset
  durations <- events_data$duration
  
  # Build variables list
  var_list <- list()
  
  for (var_name in variables) {
    var_expr <- sym(var_name)
    
    # Apply transformations if specified
    if (!is.null(regressor_spec$transform)) {
      for (trans in regressor_spec$transform) {
        trans_expr <- parse(text = paste0(trans, "(", var_name, ")"))
        var_expr <- expr(!!trans_expr)
      }
    }
    
    # For parametric regressors with basis functions
    if (!is.null(regressor_spec$basis)) {
      basis_expr <- create_parametric_basis_expression(regressor_spec$basis, var_name)
      var_expr <- basis_expr
    }
    
    var_list[[length(var_list) + 1]] <- var_expr
  }
  
  # Create hrfspec using the expressions
  hrfspec_obj <- hrfspec(
    vars = var_list,
    basis = hrf_obj,
    onsets = onsets,
    durations = durations,
    data_env = list2env(as.list(events_data), parent = baseenv())
  )
  
  return(hrfspec_obj)
}

#' Convolve Event Term with HRF
#'
#' This function convolves an event term with the specified HRF object over the sampling frame.
#'
#' @param event_term_obj An event term object created by \code{event_term()}.
#' @param hrf_obj An HRF object created by \code{create_hrf_object()}.
#' @param sampling_frame A data frame defining the time points at which to sample the convolved signal. Should contain a 'time' column.
#'
#' @return A numeric vector containing the convolved signal sampled at the times specified in \code{sampling_frame}.
#'
#' @examples
#' # Assuming 'event_term_obj' and 'hrf_obj' have been created as in previous examples
#' sampling_frame <- data.frame(time = seq(0, 30, by = 1))
#' convolved_signal <- convolve_event_term(event_term_obj, hrf_obj, sampling_frame)
#' 
#' plot(sampling_frame$time, convolved_signal, type = 'l', xlab = 'Time', ylab = 'Signal')
#'
#' @export
convolve_event_term <- function(event_term_obj, hrf_obj, sampling_frame) {
  cterm <- convolve(event_term_obj, hrf_obj, sampling_frame)
  return(cterm)
}

#' Construct Event Model from Convolved Terms
#'
#' This function assembles an event model from a list of convolved terms and a sampling frame.
#'
#' @param convolved_terms A named list of convolved terms (numeric vectors), where names correspond to regressor names.
#' @param sampling_frame A data frame defining the time points at which the signals have been sampled.
#'
#' @return An \code{event_model} object containing the design matrix, sampling frame, and terms.
#'
#' @examples
#' # Assuming 'convolved_terms' is a list of convolved signals
#' # 'sampling_frame' is a data frame with a 'time' column
#' event_model_obj <- construct_event_model(convolved_terms, sampling_frame)
#' 
#' # View the design matrix
#' head(event_model_obj$design_matrix)
#' 
#' @export
construct_event_model <- function(convolved_terms, sampling_frame) {
  design_matrix <- do.call(cbind, lapply(convolved_terms, function(cterm) cterm))
  colnames(design_matrix) <- names(convolved_terms)
  event_model <- list(
    design_matrix = as.data.frame(design_matrix),
    sampling_frame = sampling_frame,
    terms = convolved_terms
  )
  class(event_model) <- "event_model"
  return(event_model)
}

#' Build Event Model from Specification
#'
#' This function builds an entire event model from the given specification, events data, and sampling frame.
#'
#' @param spec A list containing the model specification, including 'hrfs' and 'regressors' definitions.
#' @param events_data A data frame containing event information, including onsets, durations, and any variables referenced in the regressors.
#' @param sampling_frame A data frame defining the time points at which to sample the convolved signals. Should contain a 'time' column.
#'
#' @return An \code{event_model} object containing the design matrix, sampling frame, and terms.
#'
#' @examples
#' # Define the specification
#' spec <- list(
#'   hrfs = list(
#'     canonical = list(type = "HRF_SPMG1"),
#'     gamma_hrf = list(
#'       type = "HRF_GAMMA",
#'       parameters = list(shape = 6, rate = 1)
#'     )
#'   ),
#'   regressors = list(
#'     main_effect = list(
#'       type = "hrf",
#'       variables = c("stim"),
#'       hrf = "canonical"
#'     ),
#'     rt_effect = list(
#'       type = "hrf_parametric",
#'       variables = c("RT"),
#'       hrf = "gamma_hrf"
#'     )
#'   )
#' )
#' 
#' # Events data
#' events_data <- data.frame(
#'   onset = c(0, 10, 20),
#'   duration = c(1, 1, 1),
#'   stim = c('A', 'B', 'A'),
#'   RT = c(0.5, 0.6, 0.55)
#' )
#' 
#' # Sampling frame
#' sampling_frame <- data.frame(time = seq(0, 30, by = 1))
#' 
#' # Build the event model
#' event_model_obj <- build_event_model_from_spec(spec, events_data, sampling_frame)
#' 
#' # View the design matrix
#' head(event_model_obj$design_matrix)
#'
#' @export
build_event_model_from_spec <- function(spec, events_data, sampling_frame) {
  # Build HRF objects
  hrf_objects <- lapply(spec$hrfs, create_hrf_object)
  
  # Build event terms
  event_terms <- lapply(spec$regressors, function(reg_spec) {
    construct_event_term_from_spec(reg_spec, events_data, hrf_objects)
  })
  
  # Create event_model using create_event_model
  event_model <- create_event_model(
    event_terms = event_terms,
    events = events_data,
    onsets = events_data$onset,
    block = events_data$block,
    sampling_frame = sampling_frame,
    durations = events_data$duration
  )
  
  return(event_model)
}

#' Create an HRF specification programmatically
#' 
#' @param variable Name of variable or parametric basis object
#' @param hrf HRF object to use
#' @param onsets Vector of onset times
#' @param durations Vector of durations (optional)
#' @param data Data frame containing variables (for parametric bases)
#' @return An hrfspec object
#' @keywords internal
create_hrfspec <- function(variable, hrf, onsets, durations = NULL, 
                          blockids = NULL, data = NULL) {
  # Handle parametric basis objects
  if (inherits(variable, "ParametricBasis")) {
    basis_obj <- variable
    name <- basis_obj$name
    label <- paste0("hrf(", basis_obj$name, ")")
    
    if (is.null(data)) {
      stop("Data required for parametric basis")
    }
    transformed_data <- predict(basis_obj, data)
    
  } else {
    name <- variable
    label <- paste0("hrf(", variable, ")")
    transformed_data <- NULL
  }
  
  if (is.null(blockids)) {
    blockids <- rep(1, length(onsets))
  }
  
  spec <- list(
    name = name,
    label = label,
    onsets = onsets,
    durations = durations %||% rep(0, length(onsets)),
    blockids = blockids,  # Add block IDs
    hrf = hrf,
    basis = if(inherits(variable, "ParametricBasis")) variable else NULL,
    transformed_data = transformed_data
  )
  
  class(spec) <- c("hrfspec", "list")
  spec
}

#' Create parametric basis from specification
#' 
#' @param basis_spec Basis specification from YAML
#' @param data Data frame containing the variable
#' @return A ParametricBasis object
#' @keywords internal
create_parametric_basis <- function(basis_spec, data) {
  var <- data[[basis_spec$variable]]
  
  switch(basis_spec$type,
    "Poly" = Poly(var, basis_spec$parameters$degree),
    "BSpline" = BSpline(var, basis_spec$parameters$degree),
    "Standardized" = Standardized(var),
    "Ident" = Ident(var),
    stop("Unknown basis type: ", basis_spec$type)
  )
}

create_parametric_basis_expression <- function(basis_spec, variable_name) {
  var_sym <- sym(variable_name)
  
  basis_expr <- switch(basis_spec$type,
    "Poly" = expr(Poly(!!var_sym, degree = !!basis_spec$parameters$degree)),
    "BSpline" = expr(BSpline(!!var_sym, degree = !!basis_spec$parameters$degree)),
    "Standardized" = expr(Standardized(!!var_sym)),
    "Ident" = var_sym,
    stop("Unknown basis type: ", basis_spec$type)
  )
  
  basis_expr
}



#' Build event model from configuration components
#' 
#' @param events_data Event data from bidser::read_events
#' @param events_spec Event specification from YAML
#' @param regressors_spec Regressor specification from YAML
#' @param hrfs List of HRF objects
#' @param sampling_frame Time series sampling frame
#' @return An event_model object
#' @keywords internal
build_event_model <- function(events_data, events_spec, regressors_spec, hrfs, sampling_frame) {
  # Get block information from specified column
  block_vals <- events_data[[events_spec$block]]
  if (is.factor(block_vals)) {
    block_vals <- as.numeric(as.character(block_vals))
  }
  if (!is.unsorted(block_vals)) {
    stop("Block values must be non-decreasing")
  }
  
  # Create block IDs from block values
  block_ids <- rank(block_vals)
  
  # Create hrfspec for each regressor
  var_specs <- list()
  
  for (reg_name in names(regressors_spec)) {
    reg <- regressors_spec[[reg_name]]
    
    # Get variable data and create basis if needed
    if (!is.null(reg$basis)) {
      # Create parametric basis
      basis_obj <- create_parametric_basis(
        list(
          variable = reg$variables[1],  # Assume first variable for basis
          type = reg$basis$type,
          parameters = reg$basis$parameters
        ),
        events_data
      )
      
      var_specs[[reg_name]] <- create_hrfspec(
        variable = basis_obj,
        hrf = hrfs[[reg$hrf]],
        onsets = events_data[[events_spec$onset]],
        durations = events_data[[events_spec$duration]],
        blockids = block_ids,  # Use proper block IDs
        data = events_data
      )
    } else {
      # Regular variable
      var_specs[[reg_name]] <- create_hrfspec(
        variable = reg$variables[1],
        hrf = hrfs[[reg$hrf]],
        onsets = events_data[[events_spec$onset]],
        durations = events_data[[events_spec$duration]],
        blockids = block_ids  # Use proper block IDs
      )
    }
  }
  
  # Create model spec
  model_spec <- list(
    event_table = events_data,
    onsets = events_data[[events_spec$onset]],
    event_spec = var_specs,
    blockvals = block_vals,
    blockids = block_ids,
    durations = events_data[[events_spec$duration]],
    sampling_frame = sampling_frame,
    drop_empty = TRUE,
    precision = 0.3
  )
  
  class(model_spec) <- "event_model_spec"
  construct_model(model_spec)
}

# Improve create_subject_baseline() to handle all confound cases
create_subject_baseline <- function(config, subject_id, sampling_frame) {
  # Get baseline spec with proper defaults
  baseline_spec <- config$spec$baseline %||% list(
    basis = "bspline",
    degree = 3,
    intercept = "runwise"
  )
  
  # Handle confounds with better validation
  confounds_list <- NULL
  if (!is.null(baseline_spec$confounds) || !is.null(config$spec$confounds)) {
    confounds_spec <- baseline_spec$confounds %||% config$spec$confounds
    confound_data <- read_and_validate_confounds(
      config$project,
      subject_id,
      config$tasks,
      confounds_spec
    )
    confounds_list <- process_confounds(confound_data, confounds_spec)
  }
  
  # Create baseline model with all components
  baseline_model(
    basis = baseline_spec$basis,
    degree = baseline_spec$degree,
    sframe = sampling_frame,
    intercept = baseline_spec$intercept,
    nuisance_list = confounds_list
  )
}

# Add helper functions for improved confound handling
read_and_validate_confounds <- function(project, subject_id, tasks, confounds_spec) {
  # Implementation here...
}

process_confounds <- function(confound_data, confounds_spec) {
  # Implementation here...
}

# Add function to load events for a subject
load_subject_events <- function(config, subject_id) {
  # Read events for each task
  events_by_task <- lapply(config$tasks, function(task) {
    events_data <- bidser::read_events(
      config$project,
      subid = subject_id,
      task = task
    )
    
    if (length(events_data) == 0) {
      stop(sprintf("No event files found for subject %s, task %s", 
                  subject_id, task))
    }
    
    # Extract block IDs and ensure they're numeric
    block_col <- config$events_info$mapping$block
    events <- events_data$data[[1]]
    block_vals <- events[[block_col]]
    
    if (is.factor(block_vals)) {
      block_vals <- as.numeric(as.character(block_vals))
    }
    
    if (!is.numeric(block_vals)) {
      stop(sprintf("Block values must be numeric for subject %s, task %s", 
                  subject_id, task))
    }
    
    # Add block IDs to events data
    events$blockids <- rank(block_vals)
    events
  })
  names(events_by_task) <- config$tasks
  events_by_task
}

library(rlang)

create_parametric_basis_expression <- function(basis_spec, variable_name) {
  var_sym <- sym(variable_name)
  
  basis_expr <- switch(basis_spec$type,
    "Poly" = expr(Poly(!!var_sym, degree = !!basis_spec$parameters$degree)),
    "BSpline" = expr(BSpline(!!var_sym, degree = !!basis_spec$parameters$degree)),
    "Standardized" = expr(Standardized(!!var_sym)),
    "Ident" = var_sym,
    stop("Unknown basis type: ", basis_spec$type)
  )
  
  basis_expr
}

