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
as_vectors <- function(x) { UseMethod("as_vectors") }

#' @importFrom methods setGeneric 
setGeneric("as_vectors") 


#' Construct an event model
#' 
#' This function creates an event-based fMRI regression model, represented as a data structure.
#' 
#' @importFrom lazyeval f_eval
#' @param formula_or_list The model specification, typically a `formula`. The formula should have the following format:
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
#' \describe{
#'   \item{formula}{The formula used to create the model.}
#'   \item{design}{The design matrix for the model, with one row per time point and one column per predictor variable.}
#'   \item{block_indices}{A list of indices defining the start and end time points of each block.}
#'   \item{baseline}{A vector containing the estimated baseline fMRI signal level for each block.}
#'   \item{dur}{A vector containing the duration (in seconds) of each event or block in the design.}
#' }
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
event_model <- function(formula_or_list, data, block, sampling_frame, drop_empty=TRUE, durations=0, ...) { UseMethod("event_model") }



#' get_data
#' 
#' @param x the dataset
#' @param ... extra args
#' @keywords internal
#' @noRd
get_data <- function(x, ...) UseMethod("get_data")


#' Extract data matrix
#'
#' Return the BOLD time series from a dataset as a numeric matrix. Rows
#' correspond to voxels and columns to scans.
#'
#' @param x The dataset object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A numeric matrix containing the data values.
#'
#' @examples
#' X <- matrix(rnorm(20), 4, 5)
#' dset <- matrix_dataset(X, TR = 1, run_length = 5)
#' mat <- get_data_matrix(dset)
#' dim(mat)
#'
#' @export
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
#' @export
term_matrices <- function(x, ...) UseMethod("term_matrices")



#' design_env
#' 
#' return regression design as a set of matrices stored in an environment
#' 
#' @param x the object
#' @param ... extra args
#' @keywords internal
#' @noRd
design_env <- function(x, ...) UseMethod("design_env")


#' Calculate contrast weights for a given contrast specification and term.
#'
#' @description
#' This function calculates the contrast weights based on the contrast specification
#' provided by the user. It is a generic function that dispatches to the appropriate
#' method depending on the class of the contrast specification (e.g., unit_contrast_spec,
#' pair_contrast_spec, poly_contrast_spec, etc.).
#'
#' @param x The contrast specification object
#' @param ... Extra arguments passed to specific methods
#' @return A list containing:
#' \describe{
#'     \item{term}{The model term the contrast is applied to}
#'     \item{name}{The name of the contrast}
#'     \item{weights}{A matrix of contrast weights}
#'     \item{condnames}{The condition names associated with the weights}
#'     \item{contrast_spec}{The original contrast specification}
#' }
#' @examples
#' # Create a data frame with experimental design
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 10, 20, 80),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create a sampling frame
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' 
#' # Create an event model
#' evmodel <- event_model(
#'   onsets ~ hrf(condition),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Create a contrast comparing conditions A and B
#' con <- pair_contrast(
#'   ~condition == "A",
#'   ~condition == "B",
#'   name = "A_vs_B"
#' )
#' 
#' # Calculate the contrast weights
#' weights <- contrast_weights(con, evmodel)
#' @export
#' @family contrast_weights
#' @seealso [pair_contrast()], [unit_contrast()], [poly_contrast()]
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
#' Return the experimental cells that are in a model term as a table. Experimental cells 
#' represent unique combinations of factor levels in the design. For example, if a design 
#' has factors A (levels: a1, a2) and B (levels: b1, b2), the cells would be: a1:b1, 
#' a1:b2, a2:b1, a2:b2.
#' 
#' @param x The object (typically an event_term or event_model)
#' @param drop.empty Logical; if TRUE, cells with no events are removed (default: TRUE)
#' @param ... Additional arguments passed to methods
#' @return A tibble containing the experimental cells with attributes:
#' \describe{
#'     \item{count}{Number of events in each cell}
#'     \item{rownames}{Cell names when cells have multiple factors}
#' }
#' @examples
#' # Create a simple factorial design
#' evlist <- list(
#'   fac1 = factor(c("A", "B", "A", "B")),
#'   fac2 = factor(c("1", "1", "2", "2"))
#' )
#' 
#' # Create an event term
#' eterm <- event_term(
#'   evlist,
#'   onsets = 1:4,
#'   blockids = rep(1, 4)
#' )
#' 
#' # Get the experimental cells
#' cells(eterm)  # Returns cells: A:1, A:2, B:1, B:2
#' 
#' # Create an event model
#' event_data <- data.frame(
#'   fac = c("a", "B", "A", "B"),
#'   onsets = c(1, 10, 20, 80),
#'   run = c(1, 1, 1, 1)
#' )
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' evmodel <- event_model(
#'   onsets ~ hrf(fac),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Get cells from the model
#' cells(evmodel)
#' @export
#' @family cells
#' @seealso [event_term()], [event_model()]
cells <- function(x, ...) UseMethod("cells")



#' Conditions
#' 
#' Return the set of condition labels associated with a model term. Conditions represent 
#' the unique experimental conditions in the design, typically formed from factor levels 
#' and/or basis functions. For example, a term with factor "stimulus" (levels: face, house) 
#' and two basis functions would have conditions: `"stimulus[face]:basis1"`, `"stimulus[face]:basis2"`, 
#' `"stimulus[house]:basis1"`, `"stimulus[house]:basis2"`.
#' 
#' @param x The model term (typically an event_term, event_model, or convolved_term)
#' @param drop.empty Logical; if TRUE, conditions with no events are dropped (default: TRUE)
#' @param expand_basis Logical; if TRUE, basis function labels are included (default: TRUE)
#' @param ... Additional arguments passed to methods
#' @return A character vector of condition labels
#' @examples
#' # Create a simple event model with a categorical predictor
#' event_data <- data.frame(
#'   stimulus = factor(c("face", "house", "face", "house")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create a sampling frame
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' 
#' # Create an event model with canonical HRF
#' evmodel <- event_model(
#'   onsets ~ hrf(stimulus),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Get condition labels
#' conditions(evmodel)  # Returns: c("stimulus[face]", "stimulus[house]")
#' 
#' # Create model with multiple basis functions
#' evmodel2 <- event_model(
#'   onsets ~ hrf(stimulus, basis = "fourier", nbasis = 2),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Get conditions with basis functions
#' conditions(evmodel2)  # Returns conditions with basis labels
#' @export
#' @family conditions
#' @seealso [cells()], [event_model()], [hrf()]
conditions <- function(x, ...) UseMethod("conditions")



#' Convolve a term with a hemodynamic response function
#'
#' @description
#' This function convolves an event sequence with a hemodynamic response function (HRF) 
#' over a specified time series grid. The convolution models the expected BOLD response 
#' to the events. For event-related designs, each event is convolved with the HRF and 
#' the results are summed. For block designs, the duration of each event is taken into 
#' account during convolution.
#'
#' @param x The event sequence (typically an event_term or event_model)
#' @param hrf The hemodynamic response function to use for convolution
#' @param sampling_frame The time series grid over which to sample the convolved function
#' @param drop.empty Logical; if TRUE, empty events are dropped (default: TRUE)
#' @param summate Logical; if TRUE, sum the convolved HRF over event durations (default: TRUE)
#' @param precision Numeric; precision of HRF sampling (default: 0.3)
#' @param ... Additional arguments passed to methods
#' @return A tibble containing the convolved design matrix, with columns for each condition
#' @examples
#' # Create a simple event-related design
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create a sampling frame
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' 
#' # Create an event term
#' eterm <- event_term(
#'   list(condition = event_data$condition),
#'   onsets = event_data$onsets,
#'   blockids = event_data$run
#' )
#' 
#' # Convolve with canonical HRF
#' convolved <- convolve(eterm, HRF_SPMG1, sframe)
#' 
#' # Convolve with multiple basis functions
#' convolved_fourier <- convolve(
#'   eterm, 
#'   getHRF("fourier", nbasis = 2),
#'   sframe
#' )
#' @export
#' @family convolution
#' @seealso [HRF_SPMG1()], [event_term()], [sampling_frame()]
convolve <- function(x, hrf, sampling_frame, ...) UseMethod("convolve")


#' Check if a variable is continuous
#' 
#' @description
#' Determines if a variable represents continuous (numeric) rather than categorical data.
#' For event terms, continuous variables are those that have numeric values (like amplitudes 
#' or modulators) rather than discrete factor levels. For example, reaction times would be 
#' continuous, while trial types would be categorical.
#' 
#' @param x The object to check (typically an event_term, event_seq, or event_matrix)
#' @return Logical; TRUE if the variable is continuous (numeric), FALSE if categorical
#' @examples
#' # Create event terms with different types
#' 
#' # Categorical event (factor)
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' cat_term <- event_term(
#'   list(condition = event_data$condition),
#'   onsets = event_data$onsets,
#'   blockids = event_data$run
#' )
#' is_continuous(cat_term)  # Returns: FALSE
#' 
#' # Continuous event (numeric)
#' event_data$rt <- c(0.8, 1.2, 0.9, 1.1)  # reaction times
#' cont_term <- event_term(
#'   list(rt = event_data$rt),
#'   onsets = event_data$onsets,
#'   blockids = event_data$run
#' )
#' is_continuous(cont_term)  # Returns: TRUE
#' @export
#' @family variable_type
#' @seealso [is_categorical()], [event_term()]
is_continuous <- function(x) UseMethod("is_continuous")



#' Check if a variable is categorical
#' 
#' @description
#' Determines if a variable represents categorical (factor-based) rather than continuous data.
#' For event terms, categorical variables are those that have discrete factor levels (like 
#' trial types or conditions) rather than numeric values. For example, stimulus types 
#' ("face", "house") would be categorical, while reaction times would be continuous.
#' This function is complementary to [is_continuous()].
#' 
#' @param x The object to check (typically an event_term, event_seq, or event_matrix)
#' @return Logical; TRUE if the variable is categorical (factor-based), FALSE if continuous
#' @examples
#' # Create event terms with different types
#' 
#' # Categorical event (factor)
#' event_data <- data.frame(
#'   condition = factor(c("face", "house", "face", "house")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' cat_term <- event_term(
#'   list(condition = event_data$condition),
#'   onsets = event_data$onsets,
#'   blockids = event_data$run
#' )
#' is_categorical(cat_term)  # Returns: TRUE
#' 
#' # Continuous event (numeric)
#' event_data$intensity <- c(0.8, 1.2, 0.9, 1.1)  # stimulus intensity
#' cont_term <- event_term(
#'   list(intensity = event_data$intensity),
#'   onsets = event_data$onsets,
#'   blockids = event_data$run
#' )
#' is_categorical(cont_term)  # Returns: FALSE
#' @export
#' @family variable_type
#' @seealso [is_continuous()], [event_term()]
is_categorical <- function(x) UseMethod("is_categorical")

# Note: levels is a base R generic from package stats
# We provide methods for it but don't document the generic itself



#' columns
#' 
#' return the column labels associated with the elements of a term.
#' 
#' @param x the term
#' @keywords internal
#' @noRd
longnames <- function(x) UseMethod("longnames")

#' Extract Column Names or Identifiers
#'
#' @description
#' Extract column names or identifiers from an object. For parametric basis objects,
#' this returns tokens representing the type of variables (categorical or continuous).
#'
#' @param x An object (typically a ParametricBasis)
#' @return A character vector of column identifiers
#' @export
columns <- function(x) UseMethod("columns")



#' Extract event table from a term or model
#'
#' @description
#' Extract the event table from a term or model as a data frame. The event table contains 
#' the experimental design information, with one row per event and columns for different 
#' variables (e.g., conditions, onsets, durations). For event terms, this returns the raw 
#' event data. For convolved terms, this includes any basis function expansions.
#'
#' @param x The object to extract events from (typically an event_term, convolved_term, or event_model)
#' @return A tibble containing the event information with columns for:
#'   \itemize{
#'     \item Factor variables (e.g., condition, stimulus type)
#'     \item Continuous variables (e.g., reaction times, intensities)
#'     \item Basis function expansions (if applicable)
#'   }
#' @examples
#' # Create an event term with multiple variables
#' event_data <- data.frame(
#'   condition = factor(c("face", "house", "face", "house")),
#'   rt = c(0.8, 1.2, 0.9, 1.1),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create event term
#' eterm <- event_term(
#'   list(
#'     condition = event_data$condition,
#'     rt = event_data$rt
#'   ),
#'   onsets = event_data$onsets,
#'   blockids = event_data$run
#' )
#' 
#' # Extract event table
#' etable <- event_table(eterm)
#' 
#' # Create and extract from convolved term
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' evmodel <- event_model(
#'   onsets ~ hrf(condition) + hrf(rt),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Get event table with basis expansions
#' model_events <- event_table(evmodel)
#' @export
#' @family events
#' @seealso [event_term()], [event_model()]
event_table <- function(x) UseMethod("event_table")



#' Extract event terms from a model
#'
#' @description
#' Extract the event-related terms from a model object, separating them from baseline 
#' or nuisance terms. Event terms represent the experimental conditions and parametric 
#' modulators in an fMRI design. For example, in a model with both task events 
#' (e.g., stimulus presentations) and baseline components (e.g., drift terms, motion 
#' parameters), this function returns only the task-related terms.
#'
#' @param x The model object (typically an fmri_model)
#' @return A list of event_term objects. Each event_term represents a different 
#' component of the experimental design and contains:
#'   \itemize{
#'     \item varname: Name of the term (e.g., "stimulus", "rt")
#'     \item events: List of event objects (factors or continuous variables)
#'     \item event_table: Data frame of event information
#'     \item onsets: Event onset times in seconds
#'     \item blockids: Run/block identifiers
#'     \item durations: Event durations in seconds
#'   }
#' @examples
#' # Create a model with both event and baseline terms
#' event_data <- data.frame(
#'   stimulus = factor(c("face", "house", "face", "house")),
#'   rt = c(0.8, 1.2, 0.9, 1.1),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' 
#' # Create event model
#' evmodel <- event_model(
#'   onsets ~ hrf(stimulus) + hrf(rt),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Create baseline model for drift
#' bmodel <- baseline_model(
#'   basis = "bs",
#'   degree = 3,
#'   sframe = sframe
#' )
#' 
#' # Combine into full model
#' fmodel <- fmri_model(evmodel, bmodel)
#' 
#' # Extract only the event terms
#' event_terms(fmodel)  # Returns list of stimulus and rt terms
#' @export
#' @family model_components
#' @seealso [baseline_terms()], [fmri_model()]
event_terms <- function(x) UseMethod("event_terms")



#' Extract baseline terms from a model
#'
#' @description
#' Extract the baseline and nuisance terms from a model object, separating them from 
#' experimental event terms. Baseline terms represent non-experimental components of the 
#' fMRI signal, such as:
#' \itemize{
#'   \item Drift terms (modeling scanner drift)
#'   \item Block terms (modeling run-specific baselines)
#'   \item Nuisance terms (e.g., motion parameters, physiological noise)
#' }
#'
#' @param x The model object (typically an fmri_model)
#' @return A list of baseline_term objects. Each baseline_term represents a different 
#' component of the non-experimental signal and contains:
#'   \itemize{
#'     \item varname: Name of the term (e.g., "drift", "block", "motion")
#'     \item design_matrix: Matrix of baseline regressors
#'     \item term_type: Type of baseline term ("drift", "block", or "nuisance")
#'   }
#' @examples
#' # Create a model with both event and baseline terms
#' event_data <- data.frame(
#'   stimulus = factor(c("face", "house", "face", "house")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' 
#' # Create event model
#' evmodel <- event_model(
#'   onsets ~ hrf(stimulus),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Create baseline model with drift and block terms
#' bmodel <- baseline_model(
#'   basis = "bs",    # B-spline basis for drift
#'   degree = 3,      # Cubic drift model
#'   sframe = sframe
#' )
#' 
#' # Combine into full model
#' fmodel <- fmri_model(evmodel, bmodel)
#' 
#' # Extract only the baseline terms
#' baseline_terms(fmodel)  # Returns list with drift and block terms
#' @export
#' @family model_components
#' @seealso [event_terms()], [fmri_model()], [baseline_model()]
baseline_terms <- function(x) UseMethod("baseline_terms")


#' Get term indices from a model or term
#'
#' @description
#' Get the indices that map between model terms and their corresponding columns in the 
#' design matrix. These indices are essential for:
#' \itemize{
#'   \item Extracting coefficients for specific terms
#'   \item Computing contrasts for specific model components
#'   \item Mapping between event terms and baseline terms
#'   \item Identifying which design matrix columns belong to which terms
#' }
#'
#' @param x The model or term object (typically an fmri_model, event_model, or convolved_term)
#' @param ... Additional arguments passed to methods
#' @return A named list where each element contains the column indices in the design matrix 
#' corresponding to that term. For example:
#' \itemize{
#'   \item For event terms: Indices for each experimental condition
#'   \item For baseline terms: Indices for drift and block terms
#'   \item For convolved terms: Indices for each basis function
#' }
#' @examples
#' # Create a model with multiple terms
#' event_data <- data.frame(
#'   stimulus = factor(c("face", "house", "face", "house")),
#'   rt = c(0.8, 1.2, 0.9, 1.1),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' 
#' # Create event model with multiple terms
#' evmodel <- event_model(
#'   onsets ~ hrf(stimulus) + hrf(rt, basis = "fourier", nbasis = 2),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Get indices for each term
#' indices <- term_indices(evmodel)
#' # Returns list with:
#' #  - Indices for stimulus conditions
#' #  - Indices for rt basis functions
#' 
#' # Create full model with baseline
#' bmodel <- baseline_model(basis = "bs", degree = 3, sframe = sframe)
#' fmodel <- fmri_model(evmodel, bmodel)
#' 
#' # Get indices for full model
#' full_indices <- term_indices(fmodel)
#' # Returns indices for both event and baseline terms
#' @export
#' @family model_components
#' @seealso [event_terms()], [baseline_terms()], [design_matrix()]
term_indices <- function(x, ...) UseMethod("term_indices")

#' @export
term_indices.event_model <- function(x, ...) {
  # Extract col_indices from the design matrix attribute
  col_indices <- attr(x$design_matrix, "col_indices")
  if (is.null(col_indices)) {
    warning("Event model design matrix missing 'col_indices' attribute.")
    return(NULL)
  }
  return(col_indices)
}




#' run
#'
#' Run a command or analysis step.
#'
#' @description
#' This function runs a command `x` with the provided extra arguments.
#'
#' @param x the command to run
#' @param ... extra args
#' @noRd
run <- function(x,...) UseMethod("run")


#' design_matrix
#' 
#' Extract or construct the design matrix from a model term or object. The design matrix
#' contains the predictor variables used in the model, with one row per time point and
#' one column per predictor. For event-related designs, the design matrix typically
#' contains the convolved HRF responses. For baseline terms, it contains drift and
#' nuisance regressors.
#' 
#' @details
#' ### Column Naming (Post-Refactor)
#' As of version X.Y.Z, column names in the final design matrix generated by 
#' `event_model` follow the structure: 
#' `term_tag` + `_` + `condition_tag` + `_b##` basis suffix
#' 
#' Refer to `event_model` documentation for details on `term_tag` and `condition_tag` generation.
#' The `design_matrix` methods for individual term types (like `design_matrix.event_term`)
#' return unconvolved predictors whose names may *not* yet follow this final structure.
#' The final naming is applied by `convolve.event_term` and assembled by 
#' `build_event_model_design_matrix`.
#' 
#' @param x The term or model object (typically an event_term, event_model, baseline_model, or fmri_model)
#' @param blockid Numeric vector specifying which blocks/runs to include (optional)
#' @param allrows Logical; if TRUE, return all rows including those with no events (default: FALSE)
#' @param drop.empty Logical; if TRUE, drop empty conditions (default: TRUE)
#' @param ... Additional arguments passed to methods
#' 
#' @return A tibble containing the design matrix, where:
#' \itemize{
#'   \item Rows represent time points (scans)
#'   \item Columns represent predictor variables
#'   \item Column names indicate the condition or regressor (see Details)
#' }
#' @export
#' @rdname design_matrix
#' @family design_matrices
#' @seealso [event_model()], [baseline_model()], [fmri_model()]
design_matrix <- function(x, ...) { UseMethod("design_matrix") }

#' elements
#' 
#' Return the ordered elements of a term or variable.
#' 
#' @description
#' Extract the unique elements from a term or variable in their natural order. For 
#' categorical variables (factors), this returns the factor levels. For continuous 
#' variables, this returns the unique values in ascending order. For event terms with 
#' multiple variables, this returns the combined elements.
#' 
#' @param x The term or variable object (typically an event_term, event_factor, or event_variable)
#' @param what Character; what to extract - "unique" for unique values or "all" for all values (default: "unique")
#' @param transformed Logical; if TRUE, return transformed values (default: FALSE)
#' @param ... Additional arguments passed to methods
#' @return A vector containing the ordered elements:
#'   \itemize{
#'     \item For factors: The factor levels in their defined order
#'     \item For numeric variables: Unique values in ascending order
#'     \item For event terms: Combined elements from all variables
#'   }
#' @examples
#' # Create event terms with different types
#' 
#' # Categorical variable
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "A", "B"), levels = c("B", "A")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' cat_term <- event_term(
#'   list(condition = event_data$condition),
#'   onsets = event_data$onsets,
#'   blockids = event_data$run
#' )
#' elements(cat_term)  # Returns: c("B", "A")
#' 
#' # Continuous variable
#' event_data$rt <- c(1.2, 0.8, 1.1, 0.9)
#' cont_term <- event_term(
#'   list(rt = event_data$rt),
#'   onsets = event_data$onsets,
#'   blockids = event_data$run
#' )
#' elements(cont_term)  # Returns: c(0.8, 0.9, 1.1, 1.2)
#' @export
#' @family term_properties
#' @seealso [levels()], [event_term()]
elements <- function(x, ...) UseMethod("elements")


#' correlation_map
#'
#' Create a correlation heatmap for an fMRI design matrix.
#'
#' @description
#' Generate a correlation heatmap showing the relationships between columns in a design 
#' matrix. This visualization helps identify potential collinearity between regressors 
#' in the model. For event models, it shows correlations between different conditions. 
#' For baseline models, it shows correlations between drift and nuisance terms.
#'
#' @param x The model object (event_model, baseline_model, or fmri_model)
#' @param ... Additional arguments passed to methods. Common arguments include:
#'   \describe{
#'     \item{method}{Correlation method ("pearson" or "spearman")}
#'     \item{half_matrix}{Logical; if TRUE, show only lower triangle}
#'     \item{absolute_limits}{Logical; if TRUE, set color limits to \[-1,1\]}
#'   }
#' @return A ggplot2 object containing the correlation heatmap, where:
#'   \itemize{
#'     \item Rows and columns represent model terms
#'     \item Colors indicate correlation strength (-1 to 1)
#'     \item Darker colors indicate stronger correlations
#'   }
#' @examples
#' # Create event data
#' event_data <- data.frame(
#'   condition = factor(c("face", "house", "face", "house")),
#'   rt = c(0.8, 1.2, 0.9, 1.1),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' 
#' # Create event model
#' evmodel <- event_model(
#'   onsets ~ hrf(condition) + hrf(rt),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Plot correlation map for event model
#' correlation_map(evmodel)
#' 
#' # Create baseline model
#' bmodel <- baseline_model(
#'   basis = "bs",
#'   degree = 3,
#'   sframe = sframe
#' )
#' 
#' # Plot correlation map for baseline model
#' correlation_map(bmodel)
#' 
#' # Create full model and plot combined correlations
#' fmodel <- fmri_model(evmodel, bmodel)
#' correlation_map(fmodel, method = "pearson", half_matrix = TRUE)
#' @export
#' @family visualization
#' @seealso [event_model()], [baseline_model()]
correlation_map <- function(x, ...) {
  UseMethod("correlation_map")
}



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
#' @description
#' Compute and return the fitted hemodynamic response function (HRF) for a model object. 
#' The HRF represents the expected BOLD response to neural activity. For models with 
#' multiple basis functions, this returns the combined HRF shape.
#'
#' @param x An object for which the fitted HRF should be computed
#' @param sample_at A vector of time points at which the HRF should be sampled
#' @param ... Additional arguments passed to methods
#' @return A numeric vector containing the fitted HRF values at the requested time points
#' @examples
#' # Create a simple dataset with two conditions
#' X <- matrix(rnorm(100 * 100), 100, 100)  # 100 timepoints, 100 voxels
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 25, 50, 75),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create dataset and sampling frame
#' dset <- matrix_dataset(X, TR = 2, run_length = 100, event_table = event_data)
#' sframe <- sampling_frame(blocklens = 100, TR = 2)
#' 
#' # Create event model with canonical HRF
#' evmodel <- event_model(
#'   onsets ~ hrf(condition),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Fit model
#' fit <- fmri_lm(
#'   onsets ~ hrf(condition),
#'   block = ~run,
#'   dataset = dset
#' )
#' 
#' # Get fitted HRF at specific timepoints
#' times <- seq(0, 20, by = 0.5)  # Sample from 0-20s every 0.5s
#' hrf_values <- fitted_hrf(fit, sample_at = times)
#' @export
#' @family hrf
#' @seealso [HRF_SPMG1()], [fmri_lm()]
fitted_hrf <- function(x, sample_at, ...) UseMethod("fitted_hrf")


#' Extract regressor set
#' 
#' @description
#' Extract a set of regressors from a model object. Regressors represent the predicted 
#' BOLD response for different experimental conditions or model components. For event-related 
#' designs, each regressor is typically a convolution of event onsets with an HRF. For 
#' baseline terms, regressors might represent drift or nuisance components.
#' 
#' @param x A model object that contains regressors (or can generate them)
#' @param ... Additional arguments passed to methods. Common arguments include:
#'   \describe{
#'     \item{condition}{Character; specific condition to extract regressors for}
#'     \item{block}{Numeric; specific block/run to extract regressors from}
#'     \item{basis}{Character; type of basis functions to use}
#'   }
#' @return A list of regressor objects, where each regressor contains:
#'   \itemize{
#'     \item values: Numeric vector of regressor values over time
#'     \item onsets: Original event onset times
#'     \item condition: Associated experimental condition
#'     \item block: Associated run/block number
#'   }
#' @examples
#' # Create event data with two conditions
#' event_data <- data.frame(
#'   condition = factor(c("face", "house", "face", "house")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' 
#' # Create event model with canonical HRF
#' evmodel <- event_model(
#'   onsets ~ hrf(condition),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Extract all regressors
#' reg_list <- regressors(evmodel)
#' 
#' # Create model with multiple basis functions
#' evmodel2 <- event_model(
#'   onsets ~ hrf(condition, basis = "fourier", nbasis = 2),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Extract regressors with basis functions
#' reg_list2 <- regressors(evmodel2)
#' @export
#' @family regressors
#' @seealso [event_model()], [HRF_SPMG1()], [convolve()]
regressors <- function(x, ...) UseMethod("regressors")

#' Shift a time series object
#'
#' @description
#' Apply a temporal shift to a time series object. This function shifts the values in time 
#' while preserving the structure of the object. Common uses include:
#' \describe{
#'   \item{alignment}{Aligning regressors with different temporal offsets}
#'   \item{derivatives}{Applying temporal derivatives to time series}
#'   \item{correction}{Correcting for timing differences between signals}
#' }
#'
#' @param x An object representing a time series or a time-based data structure
#' @param shift_amount Numeric; amount to shift by (positive = forward, negative = backward)
#' @param ... Additional arguments passed to methods
#' @return An object of the same class as the input, with values shifted in time:
#'   \describe{
#'     \item{Values}{Values are moved by the specified offset}
#'     \item{Structure}{Object structure and dimensions are preserved}
#'     \item{Padding}{Empty regions are filled with padding value}
#'   }
#' @examples
#' # Create a simple time series with events
#' event_data <- data.frame(
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' 
#' # Create regressor from events
#' reg <- regressor(
#'   onsets = event_data$onsets,
#'   sampling_frame = sframe
#' )
#' 
#' # Shift regressor forward by 2 seconds
#' reg_forward <- shift(reg, offset = 2)
#' 
#' # Shift regressor backward by 1 second
#' reg_backward <- shift(reg, offset = -1)
#' 
#' # Evaluate original and shifted regressors
#' times <- seq(0, 50, by = 2)
#' orig_values <- evaluate(reg, times)
#' shifted_values <- evaluate(reg_forward, times)
#' @export
#' @family time_series
#' @seealso [regressor()], [evaluate()]
shift <- function(x, ...) {
  UseMethod("shift")
}




#' Return the global onsets of an object
#' 
#' @description
#' Convert relative onset times to global (cumulative) onset times across runs. Global onsets 
#' are defined as cumulative time over runs, meaning they do not reset to zero for each run. 
#' This is useful for:
#' \itemize{
#'   \item Converting run-specific onsets to experiment-wide timings
#'   \item Aligning events across multiple runs
#'   \item Computing temporal distances between events in different runs
#' }
#'
#' @param x The object containing timing information (typically a sampling_frame)
#' @param onsets A numeric vector of relative onset times within each run/block
#' @param blockids Numeric vector specifying which block/run each onset belongs to
#' @param ... Additional arguments passed to methods
#' @return A numeric vector of global onset times where:
#'   \itemize{
#'     \item Each onset is adjusted by the cumulative duration of previous runs
#'     \item Times are in the same units as the input onsets (typically seconds)
#'     \item NA is returned for onsets that exceed their block duration
#'   }
#' @examples
#' # Create a sampling frame with three runs
#' sframe <- sampling_frame(
#'   blocklens = c(100, 100, 100),  # 100 scans per run
#'   TR = 2                         # 2 seconds per scan
#' )
#' 
#' # Define events in each run
#' run_onsets <- c(10, 20, 30)     # Events at 10s, 20s, 30s
#' run_ids <- c(1, 2, 3)           # One event per run
#' 
#' # Convert to global onsets
#' global_times <- global_onsets(
#'   sframe,
#'   onsets = run_onsets,
#'   blockids = run_ids
#' )
#' # Returns: c(10, 220, 430)
#' # Because:
#' #  - Run 1: 10s
#' #  - Run 2: 20s + (100 scans * 2s) = 220s
#' #  - Run 3: 30s + (200 scans * 2s) = 430s
#' @export
#' @family timing
#' @seealso [sampling_frame()], [event_model()]
global_onsets <- function(x, onsets, ...) UseMethod("global_onsets")



#' Return number of basis functions associated with HRF
#' 
#' @description
#' Get the number of basis functions used in a hemodynamic response function (HRF) or 
#' model term. For canonical HRFs (like SPM's canonical HRF), this returns 1. For 
#' flexible basis sets (like Fourier or B-spline bases), this returns the number of 
#' basis functions used to model the response shape.
#' 
#' @param x The object to query (typically an HRF, hrfspec, or convolved_term)
#' @param ... Additional arguments passed to methods
#' @return An integer indicating the number of basis functions:
#'   \itemize{
#'     \item 1 for canonical HRFs (e.g., SPM gamma)
#'     \item >1 for flexible basis sets (e.g., Fourier, B-spline)
#'     \item For convolved terms: number of basis functions per condition
#'   }
#' @examples
#' # Check basis functions for different HRF types
#' 
#' # Canonical HRF (single basis)
#' canonical_hrf <- HRF_SPMG1
#' nbasis(canonical_hrf)  # Returns: 1
#' 
#' # Fourier basis set
#' fourier_hrf <- getHRF("fourier", nbasis = 3)
#' nbasis(fourier_hrf)  # Returns: 3
#' 
#' # Create event model with multiple basis functions
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' 
#' # Model with Fourier basis
#' evmodel <- event_model(
#'   onsets ~ hrf(condition, basis = "fourier", nbasis = 3),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Get number of basis functions for model term
#' nbasis(evmodel)  # Returns: 3 (basis functions per condition)
#' @export
#' @family hrf
#' @seealso [HRF_SPMG1()], [event_model()]
nbasis <- function(x, ...) UseMethod("nbasis")


#' Generate penalty matrix for regularization
#'
#' @description
#' Generate a penalty matrix for regularizing HRF basis coefficients. The penalty matrix
#' encodes shape priors that discourage implausible or overly wiggly HRF estimates.
#' Different HRF types use different penalty structures:
#' 
#' \itemize{
#'   \item{FIR/B-spline bases: Roughness penalties based on discrete derivatives}
#'   \item{SPM canonical + derivatives: Differential shrinkage of derivative terms}
#'   \item{Fourier bases: Penalties on high-frequency components}
#'   \item{Default: Identity matrix (ridge penalty)}
#' }
#'
#' @param x The HRF object or basis specification
#' @param ... Additional arguments passed to specific methods
#' @return A symmetric positive definite penalty matrix of dimension nbasis(x) × nbasis(x)
#' @details
#' The penalty matrix R is used in regularized estimation as λ * h^T R h, where h are
#' the basis coefficients and λ is the regularization parameter. Well-designed penalty
#' matrices can significantly improve HRF estimation by encoding smoothness or other
#' shape constraints.
#' 
#' @examples
#' # FIR basis with smoothness penalty
#' fir_hrf <- HRF_FIR
#' R_fir <- penalty_matrix(fir_hrf)
#' 
#' # B-spline basis with second-order smoothness
#' bspline_hrf <- HRF_BSPLINE  
#' R_bspline <- penalty_matrix(bspline_hrf, order = 2)
#' 
#' # SPM canonical with derivative shrinkage
#' spmg3_hrf <- HRF_SPMG3
#' R_spmg3 <- penalty_matrix(spmg3_hrf, shrink_deriv = 4)
#' 
#' @export
#' @family hrf
#' @seealso [nbasis()], [HRF_objects]
penalty_matrix <- function(x, ...) UseMethod("penalty_matrix")


#' Combine HRF Basis with Coefficients
#'
#' Create a new HRF by linearly weighting the basis functions of an existing HRF.
#' This is useful for turning estimated basis coefficients into a single
#' functional HRF.
#'
#' @param hrf An object of class `HRF`.
#' @param h   Numeric vector of length `nbasis(hrf)` giving the weights.
#' @param name Optional name for the resulting HRF.
#' @param ... Additional arguments passed to methods.
#'
#' @return A new `HRF` object with `nbasis = 1`.
#' @export
hrf_from_coefficients <- function(hrf, h, ...) { UseMethod("hrf_from_coefficients") }


#' Reconstruction matrix for an HRF basis
#'
#' Returns a matrix \eqn{\Phi} that converts basis coefficients into a
#' sampled HRF shape.
#'
#' @param hrf An object of class `HRF`.
#' @param sframe A `sampling_frame` object or numeric vector of times.
#' @param precision Optional sampling interval in seconds when `sframe`
#'   is a `sampling_frame`. Defaults to the TR of `sframe`.
#' @param ... Additional arguments passed to methods
#' @return A numeric matrix with one column per basis function.
#' @export
reconstruction_matrix <- function(hrf, sframe, ...) { UseMethod("reconstruction_matrix") }


 
#' Return a set of data chunks
#' 
#' @description
#' Split a dataset into manageable chunks for processing. This is particularly useful 
#' for parallel processing of large fMRI datasets. Chunks can be created either by run 
#' (runwise=TRUE) or by dividing the data into a specified number of pieces. Each chunk 
#' contains a subset of the data and metadata about its position in the full dataset.
#' 
#' @param x The dataset to chunk (typically an fmri_dataset or matrix_dataset)
#' @param nchunks Integer; number of chunks to create (ignored if runwise=TRUE)
#' @param ... Additional arguments passed to methods. Common arguments include:
#'   \describe{
#'     \item{runwise}{Logical; if TRUE, create one chunk per run}
#'     \item{parallel}{Logical; if TRUE, prepare chunks for parallel processing}
#'   }
#' @return An iterator object that yields data chunks, where each chunk contains:
#'   \describe{
#'     \item{data}{Matrix of data values for this chunk}
#'     \item{chunk_num}{Index of this chunk}
#'     \item{voxel_ind}{Indices of voxels in this chunk}
#'     \item{row_ind}{Indices of timepoints in this chunk}
#'   }
#' @examples
#' # Create a simple matrix dataset
#' X <- matrix(rnorm(100 * 1000), 100, 1000)  # 100 timepoints, 1000 voxels
#' dset <- matrix_dataset(
#'   X, 
#'   TR = 2,
#'   run_length = c(50, 50)  # Two runs of 50 timepoints each
#' )
#' 
#' # Create chunks by run
#' run_chunks <- data_chunks(dset, runwise = TRUE)
#' 
#' # Process each run chunk
#' foreach::foreach(chunk = run_chunks) %do% {
#'   # chunk$data contains the data for one run
#'   # chunk$row_ind shows which timepoints are included
#'   mean_signal <- colMeans(chunk$data)
#' }
#' 
#' # Create arbitrary number of chunks
#' vox_chunks <- data_chunks(dset, nchunks = 4)
#' 
#' # Process chunks in parallel
#' foreach::foreach(chunk = vox_chunks) %dopar% {
#'   # chunk$data contains subset of voxels
#'   # chunk$voxel_ind shows which voxels are included
#'   apply(chunk$data, 2, sd)
#' }
#' @export
#' @family iterators
#' @seealso [matrix_dataset()], [fmri_dataset()], [foreach::foreach()]
data_chunks <- function(x, nchunks, ...) UseMethod("data_chunks")


#' Get event onsets from an object
#' 
#' @description
#' Extract the onset times of events from a model object. Onsets represent the timing of 
#' experimental events in an fMRI design, typically in seconds from the start of each run. 
#' These times are used to:
#' \itemize{
#'   \item Create regressors by convolving with HRF
#'   \item Verify event timing in the design
#'   \item Analyze temporal patterns of events
#' }
#'
#' @param x The object containing event information (typically an event_term or event_model)
#' @return A numeric vector of onset times in seconds, where:
#'   \itemize{
#'     \item Each value represents the start time of an event
#'     \item Times are relative to the start of each run
#'     \item Order matches the original event sequence
#'   }
#' @examples
#' # Create event data with multiple conditions
#' event_data <- data.frame(
#'   condition = factor(c("face", "house", "face", "house")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' 
#' # Create event term
#' eterm <- event_term(
#'   list(condition = event_data$condition),
#'   onsets = event_data$onsets,
#'   blockids = event_data$run
#' )
#' 
#' # Get onsets from term
#' onset_times <- onsets(eterm)  # Returns: c(1, 10, 20, 30)
#' 
#' # Create and get onsets from event model
#' evmodel <- event_model(
#'   onsets ~ hrf(condition),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' model_onsets <- onsets(evmodel)
#' @export
#' @family timing
#' @seealso [event_term()], [event_model()], [global_onsets()]
onsets <- function(x) UseMethod("onsets")


#' Get event durations from an object
#' 
#' @description
#' Extract the duration of events from a model object. Durations represent how long each 
#' event lasts in an fMRI design, typically in seconds. These are important for:
#' \itemize{
#'   \item Modeling block designs where stimuli have non-zero duration
#'   \item Creating accurate HRF convolutions for extended events
#'   \item Distinguishing between brief and sustained neural activity
#' }
#'
#' @param x The object containing event information (typically an event_term or event_model)
#' @return A numeric vector of durations in seconds, where:
#'   \itemize{
#'     \item Each value represents how long an event lasts
#'     \item Zero values indicate instantaneous events
#'     \item Order matches the corresponding event sequence
#'   }
#' @examples
#' # Create event data with varying durations
#' event_data <- data.frame(
#'   condition = factor(c("block", "event", "block", "event")),
#'   onsets = c(1, 10, 20, 30),
#'   durations = c(8, 0, 8, 0),  # 8s blocks and instantaneous events
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create event term
#' eterm <- event_term(
#'   list(condition = event_data$condition),
#'   onsets = event_data$onsets,
#'   durations = event_data$durations,
#'   blockids = event_data$run
#' )
#' 
#' # Get durations from term
#' dur <- durations(eterm)  # Returns: c(8, 0, 8, 0)
#' @export
#' @family timing
#' @seealso [onsets()], [event_term()]
durations <- function(x) UseMethod("durations")

#' Get event amplitudes from an object
#' 
#' @description
#' Extract the amplitude or intensity values associated with each event. Amplitudes 
#' represent the strength or magnitude of events and can be used to:
#' \itemize{
#'   \item Model parametric modulation of neural responses
#'   \item Weight events by their intensity or importance
#'   \item Create amplitude-modulated regressors
#' }
#'
#' @param x The object containing event information (typically an event_term or event_model)
#' @return A numeric vector of amplitude values, where:
#'   \itemize{
#'     \item Each value represents the intensity of an event
#'     \item Default value of 1 indicates unmodulated events
#'     \item Order matches the corresponding event sequence
#'   }
#' @examples
#' # Create event data with varying amplitudes
#' event_data <- data.frame(
#'   condition = factor(c("stim", "stim", "stim", "stim")),
#'   onsets = c(1, 10, 20, 30),
#'   intensity = c(0.5, 1.0, 1.5, 2.0),  # Parametrically varying intensity
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create event term with amplitudes
#' eterm <- event_term(
#'   list(condition = event_data$condition),
#'   onsets = event_data$onsets,
#'   amplitudes = event_data$intensity,
#'   blockids = event_data$run
#' )
#' 
#' # Get amplitudes from term
#' amp <- amplitudes(eterm)  # Returns: c(0.5, 1.0, 1.5, 2.0)
#' @export
#' @family event_properties
#' @seealso [event_term()], [onsets()]
amplitudes <- function(x) UseMethod("amplitudes")

#' Extract sampling times
#' 
#' @description
#' Get the sampling times for a regressor or sampling frame. These times represent when 
#' fMRI data was acquired and can be either relative (within each run) or global 
#' (cumulative across runs). Sampling times are used to:
#' \itemize{
#'   \item Evaluate regressors at scan acquisition times
#'   \item Align model predictions with data collection
#'   \item Convert between TR-based and time-based representations
#' }
#'
#' @param x The object containing timing information (typically a sampling_frame or regressor)
#' @param blockids Numeric vector specifying which blocks/runs to include (optional)
#' @param global Logical; if TRUE, return cumulative times across runs (default: FALSE)
#' @param ... Additional arguments passed to methods
#' @return A numeric vector of sampling times in seconds, where:
#'   \itemize{
#'     \item Each value represents a scan acquisition time
#'     \item Times account for TR (repetition time) spacing
#'     \item If global=FALSE, times reset at the start of each run
#'     \item If global=TRUE, times accumulate across runs
#'   }
#' @examples
#' # Create a sampling frame with multiple runs
#' sframe <- sampling_frame(
#'   blocklens = c(100, 100, 100),  # 100 scans per run
#'   TR = 2,                        # 2 seconds per scan
#'   start_time = 0                 # Start at time 0
#' )
#' 
#' # Get relative sampling times (reset each run)
#' rel_times <- samples(sframe)
#' # First few times: 0, 2, 4, 6, ... (resets each run)
#' 
#' # Get global sampling times (cumulative)
#' glob_times <- samples(sframe, global = TRUE)
#' # Shows: 0, 2, 4, ..., 198, 200, 202, ..., 598
#' 
#' # Get times for specific runs
#' run2_times <- samples(sframe, blockids = 2)
#' # Times for second run only
#' 
#' # Create regressor and get its sampling times
#' event_data <- data.frame(
#'   onsets = c(1, 10, 20),
#'   run = c(1, 1, 1)
#' )
#' reg <- regressor(
#'   onsets = event_data$onsets,
#'   sampling_frame = sframe
#' )
#' reg_times <- samples(reg)
#' @export
#' @family timing
#' @seealso [sampling_frame()], [regressor()], [global_onsets()]
samples <- function(x, ...) UseMethod("samples")

#' Split variables by block ID
#' 
#' @description
#' Split a vector or matrix of values into separate pieces based on block/run IDs. 
#' This function is useful for:
#' \itemize{
#'   \item Separating data into individual runs
#'   \item Processing blocks independently
#'   \item Analyzing run-specific patterns
#' }
#'
#' @param x The object containing data to split (typically a sampling_frame or dataset)
#' @param ... Additional arguments passed to methods
#' @return A list where each element contains data from one block:
#'   \itemize{
#'     \item List length equals number of blocks
#'     \item Each element contains values from one block
#'     \item Order matches the original block sequence
#'   }
#' @examples
#' # Create a sampling frame with multiple runs
#' sframe <- sampling_frame(
#'   blocklens = c(50, 50, 50),  # 3 runs of 50 scans each
#'   TR = 2
#' )
#' 
#' # Create some example data
#' data_values <- rnorm(150)  # 150 values (50 per run)
#' 
#' # Split data by run
#' run_data <- split_by_block(sframe, data_values)
#' # Returns list with 3 elements, each containing 50 values
#' 
#' # Create matrix dataset
#' X <- matrix(rnorm(150 * 10), 150, 10)  # 150 timepoints, 10 voxels
#' dset <- matrix_dataset(
#'   X,
#'   TR = 2,
#'   run_length = c(50, 50, 50)
#' )
#' 
#' # Split matrix data by run
#' run_matrices <- split_by_block(dset)
#' # Returns list with 3 matrices, each 50 x 10
#' @export
#' @family block_operations
#' @seealso [sampling_frame()], [blockids()], [blocklens()]
split_by_block <- function(x, ...) UseMethod("split_by_block")

#' Get block/run indices
#' 
#' @description
#' Get the block or run number associated with each scan/timepoint in the dataset. 
#' Block indices are used to:
#' \itemize{
#'   \item Track which scans belong to which runs
#'   \item Split data by experimental blocks
#'   \item Align events with their corresponding runs
#'   \item Apply run-specific processing
#' }
#'
#' @param x The object containing block information (typically a sampling_frame or dataset)
#' @return A numeric vector where:
#'   \itemize{
#'     \item Each element is the block/run ID for that scan
#'     \item IDs are sequential integers starting from 1
#'     \item Length matches the total number of scans
#'   }
#' @examples
#' # Create a sampling frame with multiple runs
#' sframe <- sampling_frame(
#'   blocklens = c(50, 75, 50),  # Different length runs
#'   TR = 2
#' )
#' 
#' # Get block IDs for all scans
#' block_ids <- blockids(sframe)
#' # Returns: c(1,1,...,1, 2,2,...,2, 3,3,...,3)
#' # 50 ones, 75 twos, 50 threes
#' 
#' # Create a matrix dataset
#' X <- matrix(rnorm(175 * 10), 175, 10)  # 175 timepoints (50+75+50), 10 voxels
#' dset <- matrix_dataset(
#'   X,
#'   TR = 2,
#'   run_length = c(50, 75, 50)
#' )
#' 
#' # Get block IDs from dataset
#' dataset_blocks <- blockids(dset)
#' 
#' # Use block IDs to split data by run
#' run_data <- split(1:nrow(X), dataset_blocks)
#' # Returns list with indices for each run
#' @export
#' @family block_operations
#' @seealso [blocklens()], [split_by_block()], [sampling_frame()]
blockids <- function(x) UseMethod("blockids")

#' Get block/run lengths
#' 
#' @description
#' Get the number of scans or timepoints in each block/run of the dataset. Block lengths 
#' are used to:
#' \itemize{
#'   \item Define the temporal structure of the experiment by specifying scan counts and timing per run
#'   \item Allocate memory for data matrices by pre-allocating arrays based on scan counts
#'   \item Validate data dimensions across runs by checking against expected lengths
#'   \item Calculate global timing information by computing cumulative timing across runs
#' }
#'
#' @param x The object containing block information (typically a sampling_frame or dataset)
#' @param ... Additional arguments passed to methods
#' @return A numeric vector where:
#' \itemize{
#'   \item Each element is the number of scans in a block or run
#'   \item Length equals the number of blocks/runs 
#'   \item Values are positive integers
#' }
#' @examples
#' # Create a sampling frame with varying run lengths
#' sframe <- sampling_frame(
#'   blocklens = c(100, 150, 100),  # Different length runs
#'   TR = 2
#' )
#' 
#' # Get number of scans per run
#' run_lengths <- blocklens(sframe)  # Returns: c(100, 150, 100)
#' 
#' # Use block lengths to create a dataset
#' total_scans <- sum(run_lengths)  # 350 total timepoints
#' X <- matrix(rnorm(total_scans * 10), total_scans, 10)  # 10 voxels
#' dset <- matrix_dataset(
#'   X,
#'   TR = 2,
#'   run_length = run_lengths
#' )
#' 
#' # Verify block lengths in dataset
#' dset_lengths <- blocklens(dset)
#' 
#' # Use lengths to create time vectors for each run
#' time_vectors <- lapply(run_lengths, function(len) seq(0, by = 2, length.out = len))
#' @export
#' @family block_operations
#' @seealso [blockids()], [split_by_block()], [sampling_frame()]
blocklens <- function(x, ...) UseMethod("blocklens")

#' Generate F-contrasts for a model term
#' 
#' @description
#' Create F-contrasts to test for overall effects of model terms. F-contrasts are used to:
#'\describe{
#'   \item{categorical}{Test for any effect of a categorical predictor}
#'   \item{basis}{Compare multiple basis functions simultaneously}
#'   \item{nonlinear}{Test for nonlinear effects of continuous predictors}
#'   \item{overall}{Evaluate overall significance of model terms}
#'}
#'
#' @param x The model term to generate contrasts for (typically an event_term or event_model)
#' @param ... Additional arguments passed to methods. Common arguments include:
#'\describe{
#'     \item{basis}{Character; type of basis functions used}
#'     \item{nbasis}{Integer; number of basis functions} 
#'     \item{exclude}{Character vector of conditions to exclude}
#'}
#' @return A list of contrast specifications where each contains:
#'\describe{
#'     \item{weights}{Matrix of contrast weights}
#'     \item{term}{The model term being tested}
#'     \item{name}{Descriptive name for the contrast}
#'     \item{df}{Degrees of freedom for the contrast}
#'}
#' @examples
#' # Create event data with multiple conditions
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "C", "A", "B", "C")),
#'   rt = c(0.8, 1.2, 0.9, 1.1, 0.7, 1.3),
#'   onsets = c(1, 10, 20, 30, 40, 50),
#'   run = c(1, 1, 1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame
#' sframe <- sampling_frame(blocklens = 60, TR = 2)
#' 
#' # Create event model with multiple terms
#' evmodel <- event_model(
#'   onsets ~ hrf(condition) + hrf(rt),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Get F-contrast for main effect of condition
#' cond_contrast <- Fcontrasts(evmodel)
#' 
#' # Create model with multiple basis functions
#' evmodel2 <- event_model(
#'   onsets ~ hrf(condition, basis = "fourier", nbasis = 3),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Get F-contrasts testing all basis functions
#' basis_contrasts <- Fcontrasts(evmodel2)
#' @export
#' @family contrasts
#' @seealso [event_model()], [contrast_weights()]
Fcontrasts <- function(x, ...) UseMethod("Fcontrasts")

#' estcon
#' 
#' @param x the object
#' @param fit the model fit
#' @param ... extra args
#' @noRd
#' @keywords internal
estcon <- function(x, fit, ...) UseMethod("estcon")

#' construct
#' 
#' construct a term given a an hrf spec and model specification
#' 
#' @param x the object
#' @param model_spec the model specification
#' @param ... extra args
#' @noRd
#' @keywords internal
construct <- function(x, model_spec,...) UseMethod("construct")




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
#' @rdname build_afni_stims
build_afni_stims <- function(x, ...) UseMethod("build_afni_stims")

#' Split Event Onsets into Lists by Factor Levels or Blocks
#'
#' Split a vector of event onsets into separate lists based on factor levels and/or block IDs.
#' This is useful for:
#' 
#' \describe{
#'   \item{separation}{Separating events by experimental conditions}
#'   \item{organization}{Organizing onsets by scanning runs/blocks}
#'   \item{preparation}{Preparing onset times for AFNI analysis}
#'   \item{analysis}{Analyzing timing patterns within conditions}
#' }
#'
#' @param x The object containing onset information (typically an event_term or event_model)
#' @param sframe A sampling_frame object defining the temporal structure
#' @param global Logical; if TRUE, convert to global onset times (default: FALSE)
#' @param blocksplit Logical; if TRUE, further split by block IDs (default: FALSE)
#' @param ... Additional arguments passed to methods
#' @return A list of numeric vectors where:
#' \describe{
#'     \item{Elements}{Each element contains onsets for one condition/block}
#'     \item{Names}{Names correspond to condition labels}
#'     \item{Nested Structure}{If blocksplit=TRUE, each condition contains a nested list of blocks}
#' }
#' @examples
#' # Create example data with multiple conditions and blocks
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "A", "B", "A", "B")),
#'   onsets = c(1, 10, 30, 40, 70, 80),
#'   run = c(1, 1, 2, 2, 3, 3)
#' )
#' 
#' # Create sampling frame
#' sframe <- sampling_frame(blocklens = c(25, 25, 25), TR = 2)
#' 
#' # Create event term
#' eterm <- event_term(
#'   list(condition = event_data$condition),
#'   onsets = event_data$onsets,
#'   blockids = event_data$run
#' )
#' 
#' # Split onsets by condition
#' split_by_cond <- split_onsets(eterm, sframe)
#' # Returns list with onsets for conditions A and B
#' 
#' # Split by condition and block
#' split_by_block <- split_onsets(eterm, sframe, blocksplit = TRUE)
#' # Returns nested list: conditions -> blocks -> onsets
#' 
#' # Get global onset times
#' split_global <- split_onsets(eterm, sframe, global = TRUE)
#' # Returns onsets adjusted for block timing
#' @family timing_functions
#' @seealso [event_term()], [sampling_frame()], [global_onsets()]
#' @export
split_onsets <- function(x, ...) UseMethod("split_onsets")



#' estimate contrast
#' 
#' @param x the contrast
#' @param fit the model fit
#' @param colind the subset of column indices in the design matrix
#' @param ... extra args
#' @noRd 
#' @keywords internal
estimate_contrast <- function(x, fit, colind, ...) UseMethod("estimate_contrast")


#' estimate a linear model sequentially for each "chunk" (a matrix of time-series) of data
#' 
#' @param x the dataset 
#' @param ... extra args
#' @noRd
#' @keywords internal
chunkwise_lm <- function(x, ...) UseMethod("chunkwise_lm")



#' Extract Standard Errors from a Model Fit
#'
#' Extract standard errors of parameter estimates from a fitted model object.
#' This is part of a family of functions for extracting statistical measures.
#'
#' @param x The fitted model object
#' @param type The type of standard errors to extract: "estimates" or "contrasts" (default: "estimates")
#' @param recon Logical; whether to reconstruct the full matrix representation (default: FALSE)
#' @param ... Additional arguments passed to methods
#' @return A tibble or matrix containing standard errors of parameter estimates
#' @examples
#' # Create example data
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame and dataset
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' dset <- matrix_dataset(
#'   matrix(rnorm(100 * 2), 100, 2),
#'   TR = 2,
#'   run_length = 50,
#'   event_table = event_data
#' )
#' 
#' # Fit model
#' fit <- fmri_lm(
#'   onsets ~ hrf(condition),
#'   block = ~run,
#'   dataset = dset
#' )
#' 
#' # Extract standard errors
#' se <- standard_error(fit)
#' @family statistical_measures
#' @export
standard_error <- function(x, ...) UseMethod("standard_error")

#' Extract Test Statistics from a Model Fit
#'
#' Extract test statistics (e.g., t-statistics, F-statistics) from a fitted model object.
#' This is part of a family of functions for extracting statistical measures.
#'
#' @param x The fitted model object
#' @param type The type of statistics to extract: "estimates", "contrasts", or "F" (default: "estimates")
#' @param ... Additional arguments passed to methods
#' @return A tibble or matrix containing test statistics
#' @examples
#' # Create example data
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame and dataset
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' dset <- matrix_dataset(
#'   matrix(rnorm(100 * 2), 100, 2),
#'   TR = 2,
#'   run_length = 50,
#'   event_table = event_data
#' )
#' 
#' # Fit model
#' fit <- fmri_lm(
#'   onsets ~ hrf(condition),
#'   block = ~run,
#'   dataset = dset
#' )
#' 
#' # Extract test statistics
#' tstats <- stats(fit)
#' @family statistical_measures
#' @export
stats <- function(x, ...) UseMethod("stats")

#' Extract P-values from a Model Fit
#'
#' Extract p-values associated with parameter estimates or test statistics from a fitted model object.
#' This is part of a family of functions for extracting statistical measures.
#'
#' @param x The fitted model object
#' @param ... Additional arguments passed to methods. Common arguments include:
#'   \describe{
#'     \item{type}{The type of p-values to extract (e.g., "estimates" or "contrasts")}
#'   }
#' @return A tibble or matrix containing p-values
#' @examples
#' # Create example data
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame and dataset
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' dset <- matrix_dataset(
#'   matrix(rnorm(100 * 2), 100, 2),
#'   TR = 2,
#'   run_length = 50,
#'   event_table = event_data
#' )
#' 
#' # Fit model
#' fit <- fmri_lm(
#'   onsets ~ hrf(condition),
#'   block = ~run,
#'   dataset = dset
#' )
#' 
#' # Extract p-values
#' pvals <- p_values(fit)
#' @family statistical_measures
#' @export
p_values <- function(x, ...) UseMethod("p_values")

#' Extract Long Names of Variable Levels
#'
#' Get the extended names of variable levels, which include the term prefix and any basis function 
#' information. Long names provide the complete specification of each condition in the model.
#' For example, if a term has conditions "level1" and "level2" with basis functions "basis1" and "basis2",
#' the long names would be "term#level1:basis1", "term#level1:basis2", "term#level2:basis1", "term#level2:basis2".
#'
#' @param x The object to extract names from (typically an event_term, event_model, or convolved_term)
#' @param ... Additional arguments passed to methods. Common arguments include:
#' \describe{
#'     \item{exclude_basis}{Logical; if TRUE, exclude basis function labels from names}
#'     \item{drop_empty}{Logical; if TRUE, drop empty condition levels}
#' }
#' @return A character vector containing the full condition names with term prefixes and basis functions
#' @examples
#' # Create example data with multiple conditions
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "C", "A", "B", "C")),
#'   rt = c(0.8, 1.2, 0.9, 1.1, 0.7, 1.3),
#'   onsets = c(1, 10, 20, 30, 40, 50),
#'   run = c(1, 1, 1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame
#' sframe <- sampling_frame(blocklens = 60, TR = 2)
#' 
#' # Create event model with multiple basis functions
#' evmodel <- event_model(
#'   onsets ~ hrf(condition, basis = "fourier", nbasis = 2),
#'   data = event_data,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#' 
#' # Get long names including basis functions
#' lnames <- longnames(evmodel)
#' # Returns: c("condition#A:basis1", "condition#A:basis2",
#' #           "condition#B:basis1", "condition#B:basis2",
#' #           "condition#C:basis1", "condition#C:basis2")
#' 
#' # Create simple event term
#' eterm <- event_term(
#'   list(condition = event_data$condition),
#'   onsets = event_data$onsets,
#'   blockids = event_data$run
#' )
#' 
#' # Get long names for term
#' term_names <- longnames(eterm)
#' # Returns: c("condition#A", "condition#B", "condition#C")
#' @family variable_names
#' @seealso [shortnames()], [event_model()], [event_term()]
#' @export
longnames <- function(x, ...) UseMethod("longnames")

#' @export
#' @rdname longnames
longnames.event_term <- function(x, ...) {
  # Get the cells (factor level combinations) for this term
  term.cells <- cells(x)
  
  # Create long names by combining variable names with their levels
  # Format: varname#level for each variable, joined with ":"
  apply(as.matrix(sapply(1:ncol(term.cells), 
                         function(i) {
                           paste0(names(term.cells)[i], "#", term.cells[[i]], sep="")
                         })), 1, paste, collapse=":")
}

#' @export
#' @rdname longnames
longnames.event_seq <- function(x, ...) {
  # Delegate to event_term method since event_term inherits from event_seq
  longnames.event_term(x, ...)
}

#' @export
#' @rdname longnames
longnames.afni_hrf_convolved_term <- function(x, ...) {
  # For AFNI terms, use the longnames of the underlying event term
  longnames(x$evterm, ...)
}

#' @export
#' @rdname longnames
longnames.afni_trialwise_convolved_term <- function(x, ...) {
  # For AFNI trialwise terms, use the longnames of the underlying event term
  longnames(x$evterm, ...)
}

#' @export
#' @rdname longnames
longnames.convolved_term <- function(x, ...) {
  # For convolved terms, delegate to the underlying event term
  longnames(x$evterm, ...)
}

#' Estimate Beta Coefficients for fMRI Data
#'
#' @description
#' Estimate beta coefficients (regression parameters) from fMRI data using various methods.
#' This function supports different estimation approaches for:
#' \describe{
#'   \item{single}{Single-trial beta estimation}
#'   \item{effects}{Fixed and random effects}
#'   \item{regularization}{Various regularization techniques}
#'   \item{hrf}{Optional HRF estimation}
#' }
#'
#' @param x The dataset object (fmri_dataset, matrix_dataset, or latent_dataset)
#' @param progress Logical; show progress bar.
#' @param ... Additional arguments passed to specific methods. Common arguments include:
#' \describe{
#'   \item{fixed}{Formula specifying fixed effects (constant across trials)}
#'   \item{ran}{Formula specifying random effects (varying by trial)}
#'   \item{block}{Formula specifying the block/run structure}
#'   \item{method}{Estimation method (e.g., "mixed", "r1", "lss", "pls")}
#'   \item{basemod}{Optional baseline model to regress out}
#'   \item{hrf_basis}{Basis functions for HRF estimation}
#'   \item{hrf_ref}{Reference HRF for initialization}
#' }
#'
#' @return A list of class "fmri_betas" containing:
#' \describe{
#'     \item{betas_fixed}{Fixed effect coefficients}
#'     \item{betas_ran}{Random (trial-wise) coefficients}
#'     \item{design_ran}{Design matrix for random effects}
#'     \item{design_fixed}{Design matrix for fixed effects}
#'     \item{design_base}{Design matrix for baseline model}
#'     \item{method_specific}{Additional components specific to the estimation method used}
#' }
#'
#' @details
#' This is a generic function with methods for different dataset types:
#' \describe{
#'   \item{fmri_dataset}{For volumetric fMRI data}
#'   \item{matrix_dataset}{For matrix-format data}
#'   \item{latent_dataset}{For dimensionality-reduced data}
#' }
#'
#' Available estimation methods include:
#' \describe{
#'   \item{mixed}{Mixed-effects model using rrBLUP}
#'   \item{r1}{Rank-1 GLM with joint HRF estimation}
#'   \item{lss}{Least-squares separate estimation}
#'   \item{pls}{Partial least squares regression}
#'   \item{ols}{Ordinary least squares}
#' }
#'
#' @examples
#' # Create example data
#' event_data <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame and dataset
#' sframe <- sampling_frame(blocklens = 50, TR = 2)
#' dset <- matrix_dataset(
#'   matrix(rnorm(100 * 2), 100, 2),
#'   TR = 2,
#'   run_length = 50,
#'   event_table = event_data
#' )
#' 
#' # Estimate betas using mixed-effects model
#' betas <- estimate_betas(
#'   dset,
#'   fixed = onset ~ hrf(condition),
#'   ran = onset ~ trialwise(),
#'   block = ~run,
#'   method = "mixed"
#' )
#'
#' @references
#' Mumford, J. A., et al. (2012). Deconvolving BOLD activation in event-related designs for multivoxel pattern classification analyses. NeuroImage, 59(3), 2636-2643.
#'
#' Pedregosa, F., et al. (2015). Data-driven HRF estimation for encoding and decoding models. NeuroImage, 104, 209-220.
#'
#' @seealso 
#' \code{\link{fmri_dataset}}, \code{\link{matrix_dataset}}, \code{\link{latent_dataset}}
#' @family model_estimation
#' @export
estimate_betas <- function(x, ...) UseMethod("estimate_betas")

#' Generate Neural Input Function from Event Timing
#'
#' Converts event timing information into a neural input function representing the underlying
#' neural activity before HRF convolution. This function is useful for:
#' 
#' \describe{
#'   \item{stimulus}{Creating stimulus functions for fMRI analysis}
#'   \item{modeling}{Modeling sustained vs. transient neural activity}
#'   \item{inputs}{Generating inputs for HRF convolution}
#'   \item{visualization}{Visualizing the temporal structure of experimental designs}
#' }
#'
#' @param x A regressor object containing event timing information
#' @param ... Additional arguments passed to methods. Common arguments include:
#' \describe{
#'     \item{start}{Numeric; start time of the input function}
#'     \item{end}{Numeric; end time of the input function} 
#'     \item{resolution}{Numeric; temporal resolution in seconds (default: 0.33)}
#' }
#'
#' @return A list containing:
#' \describe{
#'     \item{time}{Numeric vector of time points}
#'     \item{neural_input}{Numeric vector of input amplitudes at each time point}
#' }
#'
#' @examples
#' # Create a regressor with multiple events
#' reg <- regressor(
#'   onsets = c(10, 30, 50),
#'   duration = c(2, 2, 2),
#'   amplitude = c(1, 1.5, 0.8),
#'   hrf = HRF_SPMG1
#' )
#' 
#' # Generate neural input function
#' input <- neural_input(reg, start = 0, end = 60, resolution = 0.5)
#' 
#' # Plot the neural input function
#' plot(input$time, input$neural_input, type = "l",
#'      xlab = "Time (s)", ylab = "Neural Input",
#'      main = "Neural Input Function")
#' 
#' # Create regressor with varying durations
#' reg_sustained <- regressor(
#'   onsets = c(10, 30),
#'   duration = c(5, 10),  # sustained activity
#'   amplitude = c(1, 1),
#'   hrf = HRF_SPMG1
#' )
#' 
#' # Generate and compare neural inputs
#' input_sustained <- neural_input(
#'   reg_sustained,
#'   start = 0,
#'   end = 60,
#'   resolution = 0.5
#' )
#'
#' @family regressor_functions
#' @seealso 
#' \code{\link{regressor}}, \code{\link{evaluate.Reg}}, \code{\link{HRF_SPMG1}}
#' @export
neural_input <- function(x, ...) UseMethod("neural_input")

#' Visualize the entire design matrix as a heatmap
#'
#' Generate a heatmap visualization of a design matrix, showing regressor values over time.
#' This is useful for inspecting the temporal structure of fMRI design matrices.
#'
#' @param x The model object (event_model, baseline_model, or fmri_model)
#' @param ... Additional arguments passed to methods. Common arguments include:
#'   \describe{
#'     \item{rescale_cols}{Logical; if TRUE, columns are rescaled to (-1,1)}
#'     \item{block_separators}{Logical; if TRUE, draw white lines between blocks}
#'     \item{rotate_x_text}{Logical; if TRUE, rotate x-axis labels by 45 degrees}
#'   }
#' @return A ggplot2 object containing the design matrix heatmap
#' @export
#' @family visualization
#' @seealso [correlation_map()], [event_model()], [baseline_model()]
design_map <- function(x, ...) {
  UseMethod("design_map")
}

#' @export
#' @rdname conditions
conditions.afni_hrf_convolved_term <- function(x, ...) {
  # For AFNI terms, use the conditions of the underlying event term
  conditions(x$evterm, ...)
}

#' @export
#' @rdname conditions
conditions.afni_trialwise_convolved_term <- function(x, ...) {
  # For AFNI trialwise terms, use the conditions of the underlying event term
  conditions(x$evterm, ...)
}

#' @export
#' @rdname conditions
conditions.convolved_term <- function(x, ...) {
  # For regular convolved terms, use the conditions of the underlying event term
  conditions(x$evterm, ...)
}

#' @export
#' @rdname event_table
event_table.convolved_term <- function(x) {
  # For convolved terms, delegate to the underlying event term
  event_table(x$evterm)
}

#' @export
#' @rdname nbasis
nbasis.convolved_term <- function(x, ...) {
  # Get nbasis from the HRF object in the hrfspec
  hrfspec <- x$hrfspec
  if (!is.null(hrfspec) && !is.null(hrfspec$hrf)) {
    nbasis(hrfspec$hrf)
  } else {
    1L # Default fallback
  }
}

#' Short Names
#'
#' @description
#' Generate short names for model terms and conditions.
#'
#' @param x The object to generate short names for.
#' @param ... Additional arguments.
#' @return A character vector of short names.
#' @export
shortnames <- function(x, ...) UseMethod("shortnames")

#' @export
#' @rdname shortnames
shortnames.event_term <- function(x, ...) {
  # Get the cells (factor level combinations) for this term
  term.cells <- cells(x)
  
  # Create short names by combining levels with ":" separator (legacy format)
  apply(as.matrix(sapply(1:ncol(term.cells), 
                         function(i) {
                           term.cells[[i]]
                         })), 1, paste, collapse=":")
}

#' @export
#' @rdname shortnames
shortnames.event_model <- function(x, ...) {
  unlist(lapply(terms(x), shortnames))
}

#' @export
#' @rdname longnames
longnames.convolved_term <- function(x, ...) {
  # For convolved terms, delegate to the underlying event term
  longnames(x$evterm, ...)
}

#' @export
#' @rdname design_matrix
design_matrix.convolved_term <- function(x, blockid=NULL, ...) {
  if (is.null(blockid)) {
    x$design_matrix
  } else {
    keep <- blockids(x$sampling_frame) %in% blockid
    x$design_matrix[keep,]
  } 
}

#' @describeIn design_matrix Design matrix for AFNI-convolved terms
#'
#' AFNI-convolved terms rely on AFNI software for design matrix
#' construction. This method stops with an informative message.
#'
#' @export
design_matrix.afni_hrf_convolved_term <- function(x, blockid=NULL, ...) {
  stop("afni_hrf_convolved_term delegates design matrix construction to AFNI")
}
