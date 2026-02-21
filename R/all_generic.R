
#' @noRd
#' @keywords internal
with_package <- function(name) {
  if (!requireNamespace(name, quietly=TRUE)) {
    stop(paste("Please install the", name, "package to use this functionality"))
  }
}
  



# event_model generic is now imported from fmridesign
# See fmridesign-imports.R for re-exports



#' get_data
#' 
#' @param x the dataset
#' @param ... extra args
#' @keywords internal
# get_data <- function(x, ...) UseMethod("get_data") # Now imported from fmridataset


# get_data_matrix <- function(x, ...) UseMethod("get_data_matrix") # Now imported from fmridataset


# get_mask <- function(x, ...) UseMethod("get_mask") # Now imported from fmridataset


#' Retrieve the formula underlying a model object
#'
#' Generic used to expose the original modelling formula for fitted objects.
#'
#' @param x The object to extract a formula from.
#' @param ... Additional arguments passed to methods.
#' @return A formula.
#' @export
get_formula <- function(x, ...) UseMethod("get_formula")


# term_matrices generic is now imported from fmridesign
# term_matrices <- function(x, ...) UseMethod("term_matrices")



#' design_env
#' 
#' return regression design as a set of matrices stored in an environment
#' 
#' @param x the object
#' @param ... extra args
#' @keywords internal
#' @noRd
design_env <- function(x, ...) UseMethod("design_env")

#' Contrast generic
#'
#' Provide a contrast generic that dispatches on the first argument.
#' Falls back to fmridesign::contrast for non-fmri_meta classes.
#' @param x object
#' @param ... passed to methods
#' @return A contrast object with computed contrast weights and statistics
#' @examples
#' meta <- fmrireg:::.demo_fmri_meta()
#' contrast(meta, c("(Intercept)" = 1))
#' @export
contrast <- function(x, ...) UseMethod("contrast")

#' @export
contrast.default <- function(x, ...) fmridesign::contrast(x, ...)

#' Tidy generic
#' 
#' Minimal tidy generic to support fmri_meta tidy() without requiring broom.
#' @param x object
#' @param ... passed to methods
#' @return A tidy data frame with model coefficients and statistics
#' @export
tidy <- function(x, ...) UseMethod("tidy")


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
# contrast_weights generic is now imported from fmridesign
# See fmridesign-imports.R for re-exports





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
# cells generic is now imported from fmridesign
# cells <- function(x, ...) UseMethod("cells")



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
# conditions generic is now imported from fmridesign
# conditions <- function(x, ...) UseMethod("conditions")



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
#'   fmrihrf::gen_hrf("fourier", nbasis = 2),
#'   sframe
#' )
#' @export
#' @family convolution
#' @seealso [HRF_SPMG1()], [event_term()], [sampling_frame()]
# convolve generic is now imported from fmridesign
# convolve <- function(x, hrf, sampling_frame, ...) UseMethod("convolve")


# is_continuous generic is now imported from fmridesign



# is_categorical generic is now imported from fmridesign



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
#' @param ... Additional arguments passed to method-specific implementations.
#' @return A character vector of column identifiers
#' @examples
#' dat <- data.frame(
#'   onsets = c(0, 4, 8),
#'   condition = factor(c("A", "B", "A")),
#'   run = 1
#' )
#' ev <- event_model(
#'   onsets ~ hrf(condition),
#'   data = dat,
#'   block = ~ run,
#'   sampling_frame = fmrihrf::sampling_frame(blocklens = 12, TR = 2)
#' )
#' columns(ev)
#' @export
columns <- fmridesign::columns

#' @rdname columns
#' @export
columns.event_model <- function(x, ...) {
  dm <- fmridesign::design_matrix(x)
  if (is.null(dm)) {
    character(0)
  } else {
    colnames(dm)
  }
}



# event_table generic is now imported from fmridesign



# event_terms generic is now imported from fmridesign



# baseline_terms generic is now imported from fmridesign


# term_indices generic is now imported from fmridesign
# The event_model method was moved to fmridesign





# design_matrix generic is now imported from fmridesign
# design_matrix <- function(x, ...) { UseMethod("design_matrix") }

# elements generic is now imported from fmridesign
# elements <- function(x, ...) UseMethod("elements")


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
#'     \item{\code{method}}{Correlation method: "pearson" (default) or "spearman"}
#'     \item{\code{half_matrix}}{Logical; if TRUE, show only lower triangle (default: FALSE)}
#'     \item{\code{absolute_limits}}{Logical; if TRUE, set color limits to \[-1,1\] (default: TRUE)}
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
#' # Note: To create a full fmri_model and plot combined correlations,
#' # you would need an fmri_dataset object:
#' # fmodel <- fmri_model(evmodel, bmodel, dataset)
#' # correlation_map(fmodel, method = "pearson", half_matrix = TRUE)
#' @export
#' @family visualization
#' @seealso [event_model()], [baseline_model()]
correlation_map <- function(x, ...) {
  UseMethod("correlation_map")
}



# evaluate is now imported from fmrihrf


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
#' dset <- fmridataset::matrix_dataset(X, TR = 2, run_length = 100, event_table = event_data)
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
#' @family hrf
#' @seealso [HRF_SPMG1()], [fmri_lm()]
#' @export
fitted_hrf <- function(x, sample_at, ...) UseMethod("fitted_hrf")





# global_onsets <- function(x, onsets, ...) UseMethod("global_onsets") # Now imported from fmrihrf



# nbasis generic is now imported from fmrihrf package
# See fmrihrf-imports.R for the import statement





 
# data_chunks <- function(x, nchunks, ...) UseMethod("data_chunks") # Now imported from fmridataset


# onsets generic is now imported from fmridesign
# onsets <- function(x) UseMethod("onsets")


# durations generic is now imported from fmridesign
# durations <- function(x) UseMethod("durations")

# samples <- function(x, ...) UseMethod("samples") # Now imported from fmrihrf

# split_by_block generic is now imported from fmridesign
# split_by_block <- function(x, ...) UseMethod("split_by_block")

# blockids is now imported from fmrihrf

# blocklens is now imported from fmrihrf

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
# Fcontrasts generic is now imported from fmridesign

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
# split_onsets is now imported from fmridesign



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



#' Get Available Coefficient Names
#'
#' Return the names of available coefficients from a fitted model object.
#' This helps users discover which coefficient names can be passed to
#' \code{\link{coef_image}} or other extraction functions.
#'
#' @param x A fitted model object
#' @param ... Additional arguments passed to methods
#' @return A character vector of coefficient names
#' @export
#' @family statistical_measures
#' @seealso \code{\link{coef_image}}, \code{\link{coef}}
coef_names <- function(x, ...) UseMethod("coef_names")

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
#' dset <- fmridataset::matrix_dataset(
#'   matrix(rnorm(50 * 2), 50, 2),
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
#' dset <- fmridataset::matrix_dataset(
#'   matrix(rnorm(50 * 2), 50, 2),
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
#' @param type Character string specifying the type of p-values to extract. 
#'   Options typically include "estimates" for parameter estimates and "contrasts" 
#'   for contrast tests. Defaults to "estimates" in most methods.
#' @param ... Additional arguments passed to methods
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
#' dset <- fmridataset::matrix_dataset(
#'   matrix(rnorm(50 * 2), 50, 2),
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
longnames.convolved_term <- function(x, ...) {
  # For convolved terms, delegate to the underlying event term
  longnames(x$evterm, ...)
}

#' @export
#' @rdname longnames
longnames.event_model <- function(x, ...) {
  terms <- try(event_terms(x), silent = TRUE)
  if (inherits(terms, "try-error") || is.null(terms)) {
    return(character(0))
  }
  unname(unlist(lapply(terms, function(tt) longnames(tt, ...))))
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
#'   onset = c(1, 10, 20, 30),
#'   run = c(1, 1, 1, 1)
#' )
#' 
#' # Create sampling frame and dataset
#' sframe <- sampling_frame(blocklens = 100, TR = 2)
#' dset <- fmridataset::matrix_dataset(
#'   matrix(rnorm(100 * 2), 100, 2),
#'   TR = 2,
#'   run_length = 100,
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
# design_map generic moved to fmridesign package



#' @export
#' @rdname conditions
conditions.convolved_term <- function(x, ...) {
  # For regular convolved terms, use the conditions of the underlying event term
  conditions(x$evterm, ...)
}

#' Extract event table from convolved term
#'
#' Extract the event table from a convolved term object.
#'
#' @param x A convolved_term object
#' @param ... Additional arguments passed to the underlying event_table method
#' @return A data.frame containing the event table
#' @method event_table convolved_term
#' @export
event_table.convolved_term <- function(x, ...) {
  # For convolved terms, delegate to the underlying event term
  event_table(x$evterm, ...)
}

#' @export
#' @method nbasis convolved_term
nbasis.convolved_term <- function(x, ...) {
  # Get nbasis from the HRF object in the hrfspec
  hrfspec <- x$hrfspec
  if (!is.null(hrfspec) && !is.null(hrfspec$hrf)) {
    fmrihrf::nbasis(hrfspec$hrf)
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

#' Design Matrix for Convolved Terms
#' 
#' Extract the design matrix from a convolved term object, optionally filtered by block ID.
#' 
#' @param x A convolved_term object
#' @param blockid Optional numeric vector specifying which blocks/runs to include
#' @param ... Additional arguments (not used)
#' @return A matrix containing the convolved design matrix
#' @method design_matrix convolved_term
#' @export
design_matrix.convolved_term <- function(x, blockid=NULL, ...) {
  if (is.null(blockid)) {
    x$design_matrix
  } else {
    keep <- fmrihrf::blockids(x$sampling_frame) %in% blockid
    x$design_matrix[keep,]
  } 
}

#' Block IDs for event_model
#'
#' Return the run/block IDs associated with an event_model's sampling frame.
#'
#' @param x An event_model object
#' @param ... Additional arguments passed through
#' @return Integer vector of block IDs
#' @examples
#' ev <- fmrireg:::.demo_event_model()
#' blockids(ev)
#' @method blockids event_model
#' @export
blockids.event_model <- function(x, ...) {
  fmrihrf::blockids(x$sampling_frame, ...)
}


#' Write Results from fMRI Analysis
#'
#' Generic function to export statistical maps and analysis results from fitted fMRI models
#' to standardized file formats with appropriate metadata.
#'
#' @param x A fitted fMRI model object
#' @param ... Additional arguments passed to methods
#' @return Invisible list of created file paths
#' @export
#' @family result_export
write_results <- function(x, ...) UseMethod("write_results")
