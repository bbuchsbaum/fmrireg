#' @noRd
#' @keywords internal
get_methods <- function(obj) {
  unique(purrr::map_chr(class(obj), ~ methods(class= . )))
}


#' @noRd
#' @keywords internal
with_package <- function(name) {
  if (!requireNamespace(name, quietly=TRUE)) {
    stop(paste("Please install the", name, "package to use this functionality"))
  }
}
  


#' @noRd
#' @keywords internal
as_vectors <- function(x) { function(x, ...) UseMethod("as_vectors") }

#' @importFrom methods setGeneric 
setGeneric("as_vectors") 



#' Construct an event model
#' 
#' This function creates an event-based fMRI regression model, represented as a data structure.
#' 
#' @importFrom lazyeval f_eval f_rhs f_lhs
#' @param x The model specification, typically a `formula`. The formula should have the following format:
#'    response ~ predictor1 + predictor2 + ... + predictorN
#' where `response` is a numeric vector of fMRI signal values, and `predictor1` to `predictorN` are
#' predictor variables. Each predictor variable should be specified as a function of categorical
#' variables (factors) and/or continuous variables. The functions should have the prefix "hrf", and can
#' be defined using the `hrf()` function (see `hrf` documentation for details).
#' @param data A data frame containing the experimental design, with one row per time point and
#' one column per variable used in the model formula. If a categorical variable is used in the formula,
#' it should be a factor in the data frame. The data frame should also contain a column with the fMRI
#' signal values (the response variable).
#' @param block A formula specifying the block structure of the design. This formula should have
#' the format `block_var1 + block_var2 + ...`, where each `block_var` is a categorical variable (factor)
#' used to define blocks of time points. The block structure is used to estimate the baseline fMRI
#' signal level separately for each block.
#' @param sampling_frame A sampling frame defining the temporal and block structure of the design.
#' This should be an object of class `sampling_frame` (see `sampling_frame` documentation for details).
#' @param drop_empty Logical value indicating whether to drop empty factor levels in the model.
#' If `TRUE` (default), any factor levels with no observations will be dropped from the model. If `FALSE`,
#' empty levels will be retained and will receive a coefficient of zero in the model.
#' @param durations A numeric vector specifying the duration (in seconds) of each event in the model.
#' If the model contains block variables, the duration of each block should be specified as well.
#' The length of this vector should be equal to the number of events/blocks in the design.
#' Default value is 0 (no duration).
#' @param ... Additional arguments to be passed to methods. Currently not used.
#' 
#' @export
#' 
#' @return A list containing the following elements:
#' \item{formula}{The formula used to create the model.}
#' \item{design}{The design matrix for the model, with one row per time point and one column per predictor variable.}
#' \item{block_indices}{A list of indices defining the start and end time points of each block.}
#' \item{baseline}{A vector containing the estimated baseline fMRI signal level for each block.}
#' \item{dur}{A vector containing the duration (in seconds) of each event or block in the design.}
#' 
#' @examples 
#' # Create a data frame with experimental design
#' event_data <- data.frame(fac=c("a", "B", "A", "B"), onsets=c(1,10,20,80), run=c(1,1,1,1))
#' 
#' # Create a sampling frame with 50-second blocks and a TR of 2 seconds
#' sframe <- sampling_frame(blocklens=50, TR=2)
#' 
#' # Create an event model using the `onsets` variable as a predictor, 
#' #  with a separate baseline for each run
#' evmodel <- event_model(onsets ~ hrf(onsets), data=event_data, block=~run, sampling_frame=sframe)
#' dmat <- design_matrix(evmodel)
event_model <- function(x, data, block, sampling_frame, drop_empty=TRUE, durations=0, ...) { UseMethod("event_model") }

#' get_data
#' 
#' @param x the dataset
#' @param ... extra args
#' @keywords internal
#' @noRd
get_data <- function(x, ...) UseMethod("get_data")


#' get_data_matrix
#' 
#' @param x the dataset
#' @param ... extra args
#' @keywords internal
#' @noRd
get_data_matrix <- function(x, ...) UseMethod("get_data_matrix")


#' get_mask
#' 
#' get the binary inclusion mask associated with a dataset
#' 
#' @param x the dataset
#' @param ... extra args
#' @keywords internal
#' @noRd
get_mask <- function(x, ...) UseMethod("get_mask")




#' get_formula
#' 
#' @param x the object
#' @param ... extra args
#' @keywords internal
#' @noRd
get_formula <- function(x, ...) UseMethod("get_formula")



#' term_matrices
#' 
#' @param x the object
#' @param ... extra args
#' @keywords internal
#' @noRd
term_matrices <- function(x, ...) UseMethod("term_matrices")


#' extract long names of variable
#' 
#' get the extended names of a set of variable levels
#' 
#' @param x the object
#' @param ... extra args
#' @export
longnames <- function(x, ...) UseMethod("longnames")


#' extract short short names of variable
#' 
#' get the short names of a set of variable levels
#' 
#' @param x the object
#' @param ... extra args
#' @export
shortnames <- function(x, ...) UseMethod("shortnames")


#' design_env
#' 
#' return regression design as a set of matrices stored in an environment
#' 
#' @param x the object
#' @param ... extra args
#' @keywords internal
#' @noRd
design_env <- function(x, ...) UseMethod("design_env")


#' contrast_weights
#'
#' Calculate contrast weights for a given contrast specification and term.
#'
#' @description
#' This function calculates the contrast weights based on the contrast specification
#' provided by the user. It is a generic function that dispatches to the appropriate
#' method depending on the class of the contrast
#'
#' @param x the object
#' @param ... extra args
#' @export
contrast_weights <- function(x, ...) UseMethod("contrast_weights")


#' parent_terms
#' 
#' @param x the object
#' @keywords internal
#' @noRd
parent_terms <- function(x) UseMethod("parent_terms")


#' term_names
#' @param x the object to extra term names from
#' @noRd
#' @keywords internal
term_names <- function(x) UseMethod("term_names")



#' The experimental cells of a design
#' 
#' return the experimental cells that are in a model term as a table
#' 
#' @param x the object
#' @param ... extra args
#' @export
cells <- function(x, ...) UseMethod("cells")



#' Conditions
#' 
#' return the set of condition labels associated with a model term
#' 
#' @param x the model term
#' @param ... extra args
#' @export
#' @family conditions
conditions <- function(x, ...) UseMethod("conditions")



#' convolve
#'
#' Convolve a term `x` with a hemodynamic response function (HRF) over a `sampling_frame`.
#'
#' @description
#' This function convolves an event sequence `x` with a hemodynamic response function `hrf`
#' over a specified time series grid `sampling_frame`. It is a generic function that dispatches
#' to the appropriate method
#'
#' @param x an event sequence
#' @param hrf a hemodynamic response function
#' @param sampling_frame the time series grid over which to sample the function.
#' @param ... extra args
#' @export
convolve <- function(x, hrf, sampling_frame,...) UseMethod("convolve")


#' is_continuous
#' 
#' checks to see if a variable is continuous e.g. numeric/non-categorical
#' 
#' @param x the object 
#' @export
is_continuous <- function(x) UseMethod("is_continuous")



#' is_categorical
#' 
#' checks to see if a variable is cateogircal e.g. factor-based
#' 
#' @param x the object
#' @export
is_categorical <- function(x) UseMethod("is_categorical")

#' levels
#' 
#' extract the levels of a term
#' @param x the term
#' @export
levels <- function(x) UseMethod("levels")



#' columns
#' 
#' return the column labels associated with the elements of a term.
#' 
#' @param x the term
#' @keywords internal
#' @noRd
columns <- function(x) UseMethod("columns")



#' event_table
#'
#' Extract event_table as a `data.frame`.
#'
#' @description
#' This function extracts the event table from the provided term `x` and returns it
#' as a `data.frame`.
#'
#' @param x the term
#' @return A `data.frame` containing the event table.
#' @export
event_table <- function(x) UseMethod("event_table")



#' event_terms
#'
#' Extract a list of 'event_term's from the model object.
#'
#' @description
#' This function extracts a list of 'event_term's from the provided model object `x`.
#'
#' @param x the model object
#' @return A list containing the 'event_term's.
#' @export
event_terms <- function(x) UseMethod("event_terms")



#' baseline_terms
#'
#' Extract a list of 'baseline_term's from the model object.
#'
#' @description
#' This function extracts a list of 'baseline_term's from the provided model object `x`.
#'
#' @param x the model object
#' @return A list containing the 'baseline_term's.
#' @export
baseline_terms <- function(x) UseMethod("baseline_terms")


#' term_indices
#'
#' Get the term indices for the associated design matrix.
#'
#' @description
#' This function retrieves the term indices for the associated design matrix
#' from the provided model or term object `x`.
#'
#' @param x the model/term object
#' @param ... additional arguments
#' @return A vector of term indices.
#' @export
term_indices <- function(x,...) UseMethod("term_indices")




#' run
#'
#' Run a command.
#'
#' @description
#' This function runs a command `x` with the provided extra arguments.
#'
#' @param x the command to run
#' @param ... extra args
#' @export
run <- function(x,...) UseMethod("run")


#' design_matrix
#' 
#' construct a design matrix from the term
#' 
#' @param x the term
#' @param ... additional arguments
#' @export
design_matrix <- function(x, ...) UseMethod("design_matrix")

#' elements
#' 
#' return the ordered elements of a term/variable
#' 
#' @param x the term
#' @param ... additional arguments
#' @export
elements <- function(x, ...) UseMethod("elements")


#' Evaluate a regressor object over a time grid
#' 
#' Generic function to evaluate a regressor object over a specified time grid.
#' Different types of regressors may have different evaluation methods.
#'
#' @param x The regressor object to evaluate
#' @param grid A numeric vector specifying the time points at which to evaluate the regressor
#' @param ... Additional arguments passed to specific methods
#' @return A numeric vector or matrix containing the evaluated regressor values
#' @seealso [single_trial_regressor()], [regressor()]
#' @export
evaluate <- function(x, grid, ...) {
  UseMethod("evaluate")
}


#' fitted_hrf
#'
#' This generic function computes the fitted hemodynamic response function (HRF) for an object.
#' The method needs to be implemented for specific object types.
#'
#' @param x An object for which the fitted HRF should be computed.
#' @param sample_at A vector of time points at which the HRF should be sampled.
#' @param ... Additional arguments to be passed to specific methods.
#' @return A fitted HRF for the object.
#' @export
fitted_hrf <- function(x, sample_at, ...) UseMethod("fitted_hrf")


#' extract regressor set
#' 
#' @param x a model object that contains regressors (or can generate them)
#' @param ... extra args
#' @family regressors
regressors <- function(x, ...) UseMethod("regressors")

#' Shift a time series object
#'
#' This is a generic function to shift time series objects. The goal is to
#' provide a simple way to apply time shifts to various time series objects,
#' such as regressors and time series data.
#'
#' @param x An object representing a time series or a time-based data structure.
#' @param ... extra args
#'
#' @return An object of the same class as the input, shifted by the specified amount.
#'
#' @examples
#' \dontrun{
#' # Shift a regressor object
#' shifted_regressor <- shift(my_regressor, 5)
#'
#' # Shift a time series object
#' shifted_time_series <- shift(my_time_series, -2)
#' }
#'
#' @export
#' @family shift
shift <- function(x, ...) {
  UseMethod("shift")
}




#' return the "global" onsets of an object.
#' 
#' global onsets are defined as cumulative time over runs, i.e. it does not reset to zero for each run.
#' 
#' @export
#' @param x the object
#' @param onsets the relative onset times of the events
#' @param ... extra args
#' @family global_onsets
global_onsets <-  function(x, onsets,...) UseMethod("global_onsets")



#' return number of basis functions associated with hrf.
#' 
#' @export
#' @param x the object
nbasis <-  function(x) UseMethod("nbasis")


 
#' return a set of data chunks
#' 
#' @param x the dataset
#' @param nchunks the number of data chunks
#' @param ... extra args
#' @return an \code{iterator} returning on data chunk per iteration
#' @export
data_chunks <- function(x, nchunks, ...) UseMethod("data_chunks")


#' get event onsets of a variable
#' 
#' return an `onset` vector
#' 
#' @param x the object
#' @export
onsets <-  function(x) UseMethod("onsets")


#' get event durations of a variable
#' 
#' return a `durations` vector
#' 
#' @param x the object
#' @export
durations <-  function(x) UseMethod("durations")

#' get amplitude vector
#' 
#' get the intensity or amplitudes associated with each event
#' 
#' @param x the object
#' @export
amplitudes <-  function(x) UseMethod("amplitudes")

#' extract samples
#' 
#' get the sampling times for a regressor or sampling frame
#' 
#' @param x the object
#' @param ... extra args
#' @export
samples <-  function(x, ...) UseMethod("samples")

#' split_by_block
#' 
#' split variables by block id
#' 
#' @param x the object
#' @param ... extra args
#' @export
split_by_block  <-  function(x, ...) UseMethod("split_by_block")

#' get the block indices
#' 
#' get the block number associated with each image/scan
#' 
#' @param x the object
#' @export
blockids  <-  function(x) UseMethod("blockids")

#' get block lengths
#' 
#' get the number of scans/images in each block
#' 
#' @param x the object
#' @param ... extra args
#' @export
blocklens  <-  function(x, ...) UseMethod("blocklens")

#' Fcontrasts
#' 
#' generate a set of Fcontrasts for a model term
#' 
#' @param x the term
#' @param ... extra args
#' @export
Fcontrasts <- function(x, ...) UseMethod("Fcontrasts")

#' estcon
#' 
#' @param x the object
#' @param fit the model fit
#' @param ... extra args
#' @export
estcon <- function(x, fit, ...) UseMethod("estcon")

#' construct
#' 
#' construct a term given a an hrf spec and model specification
#' 
#' @param x the object
#' @param model_spec the model specification
#' @param ... extra args
#' @export
construct <- function(x, model_spec,...) UseMethod("construct")



#' standard_error
#' 
#' extract standard errors from fitted model
#' 
#' @param x the object
#' @param ... extra arguments
standard_error <- function(x, ...) UseMethod("standard_error")


#' stats
#' 
#' extract statistics from fitted model
#' 
#' @param x the object
#' @param ... extra arguments
#' @export
stats <- function(x, ...) UseMethod("stats")


#' p_values
#' 
#' extract p-values from fitted model
#' 
#' @param x the object
#' @param ... extra arguments
p_values <- function(x, ...) UseMethod("p_values")


#' longnames
#' 
#' @param x the object
#' @param ... extra args
#' @export
longnames <- function(x, ...) UseMethod("longnames")


#' generate an AFNI linear model command from a configuration file
#' 
#' @param x the config file
#' @param ... extra args
#' @export
gen_afni_lm <- function(x, ...) UseMethod("gen_afni_lm")

#' generate a set of AFNI stimuli for '3dDeconvolve'
#' 
#' @param x the term
#' @param ... extra args
#' @keywords internal
#' @noRd
build_afni_stims <- function(x, ...) UseMethod("build_afni_stims")

#' split an onset vector into a list
#' 
#' @param x the object associated with the onset vector
#' @param ... extra args
#' @export
split_onsets <- function(x, ...) UseMethod("split_onsets")



#' estimate trialwise beta coefficients for a dataset
#' 
#' @param x the dataset 
#' @param ... extra args
#' @export
#' @family estimate_betas
estimate_betas <- function(x, ...) UseMethod("estimate_betas")


#' estimate contrast
#' 
#' @param x the contrast
#' @param fit the model fit
#' @param colind the subset of column indices in the design matrix
#' @param ... extra args
#' @export
#' @family estimate_contrast
estimate_contrast <- function(x, fit, colind, ...) UseMethod("estimate_contrast")


#' estimate a linear model sequentially for each "chunk" (a matrix of time-series) of data
#' 
#' @param x the dataset 
#' @param ... extra args
chunkwise_lm <- function(x, ...) UseMethod("chunkwise_lm")

#' neural input
#' 
#' @param x a regressor
#' @param ... extra args
#' @family neural_input
neural_input <- function(x, ...) UseMethod("neural_input")
