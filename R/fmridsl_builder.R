###############################################################################
#                 Create HRF Object from Specification
###############################################################################

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
  # Get all available HRF functions from your package
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
  
  # Otherwise, we have a built-in HRF in the namespace
  hrf_fun <- get(type, envir = asNamespace("fmrireg"))
  
  # If no parameters, return the default HRF
  if (is.null(hrf_spec$parameters)) {
    return(hrf_fun)
  }
  
  # Validate parameters if needed
  # (calls validate_hrf_params or something if you want)
  
  # Create HRF with parameters
  switch(type,
         "HRF_GAUSSIAN" = {
           params <- list(
             mean = hrf_spec$parameters$mean %||% 5,
             sd   = hrf_spec$parameters$sd   %||% 1
           )
           do.call(gen_hrf, c(list(hrf_gaussian), params))
         },
         
         "HRF_GAMMA" = {
           params <- list(
             shape = hrf_spec$parameters$shape %||% 6,
             rate  = hrf_spec$parameters$rate  %||% 1,
             scale = hrf_spec$parameters$scale %||% 1
           )
           do.call(gen_hrf, c(list(hrf_gamma), params))
         },
         
         "HRF_BSPLINE" = {
           params <- list(
             knots  = hrf_spec$parameters$knots  %||% seq(0, 32, by=4),
             degree = hrf_spec$parameters$degree %||% 3
           )
           do.call(gen_hrf, c(list(hrf_bspline), params))
         },
         
         # If SPMG1/2/3 or other built-in
         hrf_fun
  )
}


###############################################################################
#   Regressor / Event-Model Builders & Parametric Expansions
###############################################################################

#' Validate HRF parameters
#' 
#' (Optional function, if you want to extend param checking)
#'
#' @keywords internal
validate_hrf_params <- function(type, params) {
  # This can be extended as needed
  TRUE
}

#' Construct Event Term from Regressor Specification
#'
#' This function constructs an event term based on the regressor specification,
#' events data, and HRF objects.
#'
#' @param regressor_spec A list containing the specification of a single regressor
#' @param events_data A data frame containing event information
#' @param hrf_objects A named list of HRF objects
#'
#' @return An hrfspec object or similar
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
  
  # Build variables list (placeholder logic)
  var_list <- list()
  
  for (var_name in variables) {
    var_sym <- rlang::sym(var_name)
    var_list[[length(var_list)+1]] <- var_sym
  }
  
  # Create hrfspec using expressions
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
#' This function convolves an event term with the specified HRF object
#' over the sampling frame.
#'
#' @param event_term_obj An event term object
#' @param hrf_obj An HRF object
#' @param sampling_frame A data frame defining the time points, must have 'time' col
#'
#' @return A numeric vector containing the convolved signal
#' @export
convolve_event_term <- function(event_term_obj, hrf_obj, sampling_frame) {
  cterm <- convolve(event_term_obj, hrf_obj, sampling_frame)
  return(cterm)
}


#' Construct Event Model from Convolved Terms
#'
#' @param convolved_terms A named list of convolved terms (numeric vectors)
#' @param sampling_frame A data frame defining time points
#'
#' @return An \code{event_model} object
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
#' This function builds an entire event model from the given specification,
#' events data, and sampling frame.
#'
#' @param spec A list containing 'hrfs' and 'regressors' definitions
#' @param events_data A data frame with event info
#' @param sampling_frame A data frame with 'time' col
#' @return An \code{event_model} object
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
  
  event_model
}

#' Create an HRF specification programmatically
#'
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
    blockids = blockids,
    hrf = hrf,
    basis = if (inherits(variable, "ParametricBasis")) variable else NULL,
    transformed_data = transformed_data
  )
  
  class(spec) <- c("hrfspec", "list")
  spec
}

#' Create parametric basis from specification
#' 
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

#' create_parametric_basis_expression
#'
#' @keywords internal
create_parametric_basis_expression <- function(basis_spec, variable_name) {
  var_sym <- rlang::sym(variable_name)
  
  basis_expr <- switch(basis_spec$type,
                       "Poly" = rlang::expr(Poly(!!var_sym, degree = !!basis_spec$parameters$degree)),
                       "BSpline" = rlang::expr(BSpline(!!var_sym, degree = !!basis_spec$parameters$degree)),
                       "Standardized" = rlang::expr(Standardized(!!var_sym)),
                       "Ident" = var_sym,
                       stop("Unknown basis type: ", basis_spec$type)
  )
  
  basis_expr
}

#' Build event model from configuration components
#'
#' @keywords internal
build_event_model <- function(events_data, events_spec, regressors_spec, hrfs, sampling_frame) {
  block_vals <- events_data[[events_spec$block]]
  if (is.factor(block_vals)) {
    block_vals <- as.numeric(as.character(block_vals))
  }
  
  if (!is.numeric(block_vals)) {
    stop("Block values must be numeric")
  }
  
  if (is.unsorted(block_vals, strictly=FALSE)) {
    stop("Block values must be non-decreasing in build_event_model.")
  }
  
  block_ids <- rank(block_vals)
  
  var_specs <- list()
  for (reg_name in names(regressors_spec)) {
    reg <- regressors_spec[[reg_name]]
    if (!is.null(reg$basis)) {
      basis_obj <- create_parametric_basis(
        list(
          variable = reg$variables[1],
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
        blockids = block_ids,
        data = events_data
      )
    } else {
      var_specs[[reg_name]] <- create_hrfspec(
        variable = reg$variables[1],
        hrf = hrfs[[reg$hrf]],
        onsets = events_data[[events_spec$onset]],
        durations = events_data[[events_spec$duration]],
        blockids = block_ids
      )
    }
  }
  
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


###############################################################################
#           Additional subject baseline & confounds logic (optional)
###############################################################################

create_subject_baseline <- function(config, subject_id, sampling_frame) {
  baseline_spec <- config$spec$baseline %||% list(
    basis = "bspline",
    degree = 3,
    intercept = "runwise"
  )
  
  confounds_list <- NULL
  # If baseline spec or config has confounds
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
  
  # This references a baseline_model function from your code
  baseline_model(
    basis = baseline_spec$basis,
    degree = baseline_spec$degree,
    sframe = sampling_frame,
    intercept = baseline_spec$intercept,
    nuisance_list = confounds_list
  )
}

# Example placeholders for confound logic
read_and_validate_confounds <- function(project, subject_id, tasks, confounds_spec) {
  # Implementation up to you. Possibly read confound tsv for each run in tasks, etc.
  # Return merged or list of confounds
  list()
}

process_confounds <- function(confound_data, confounds_spec) {
  # Example. You might standardize them or just pass them through.
  confound_data
}

# Add function to load events for a subject
load_subject_events <- function(config, subject_id) {
  events_by_task <- lapply(config$tasks, function(task) {
    events_data <- bidser::read_events(config$project, subid = subject_id, task = task)
    
    if (length(events_data) == 0) {
      stop(sprintf("No event files found for subject %s, task %s", subject_id, task))
    }
    
    block_col <- config$events_info$mapping$block
    events <- events_data$data[[1]]
    block_vals <- events[[block_col]]
    if (is.factor(block_vals)) {
      block_vals <- as.numeric(as.character(block_vals))
    }
    
    if (!is.numeric(block_vals)) {
      stop(sprintf("Block values must be numeric for subject %s, task %s", subject_id, task))
    }
    events$blockids <- rank(block_vals)
    events
  })
  names(events_by_task) <- config$tasks
  
  events_by_task
}
