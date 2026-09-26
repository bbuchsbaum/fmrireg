
#' @noRd
#' @keywords internal
with_package <- function(name) {
  if (!requireNamespace(name, quietly=TRUE)) {
    stop(paste("Please install the", name, "package to use this functionality"))
  }
}
  



# event_model generic is now imported from fmridesign
# See fmridesign-imports.R for re-exports



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


# contrast_weights generic is now imported from fmridesign
# See fmridesign-imports.R for re-exports





# cells generic is now imported from fmridesign
# cells <- function(x, ...) UseMethod("cells")



# conditions generic is now imported from fmridesign
# conditions <- function(x, ...) UseMethod("conditions")



# convolve generic is now imported from fmridesign
# convolve <- function(x, hrf, sampling_frame, ...) UseMethod("convolve")


# is_continuous generic is now imported from fmridesign



# is_categorical generic is now imported from fmridesign






# event_table generic is now imported from fmridesign



# event_terms generic is now imported from fmridesign



# baseline_terms generic is now imported from fmridesign


# term_indices generic is now imported from fmridesign
# The event_model method was moved to fmridesign





# design_matrix generic is now imported from fmridesign
# design_matrix <- function(x, ...) { UseMethod("design_matrix") }

# elements generic is now imported from fmridesign
# elements <- function(x, ...) UseMethod("elements")





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
#' dset <- matrix_frame(X, TR = 2, run_length = 100, event_table = event_data)
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

# Fcontrasts generic is now imported from fmridesign





#' generate an AFNI linear model command from a configuration file


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
#' dset <- matrix_frame(
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
#' dset <- matrix_frame(
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
#' dset <- matrix_frame(
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
#' @param x The dataset: an `fmri_frame` (see [matrix_frame()], [neurovec_frame()],
#'   [nifti_frame()], [latent_frame()])
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
#' This is a generic function whose `fmri_frame` method adapts to the frame's
#' feature space: volumetric frames (`volume_space`) return `NeuroVec` betas,
#' while matrix-format (`index_space`) and latent (`basis_space`) frames return
#' coefficient matrices.
#'
#' Available estimation methods include:
#' \describe{
#'   \item{mixed}{Mixed-effects model using ridge/BLUP estimation}
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
#' dset <- matrix_frame(
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
#' \code{\link{matrix_frame}}, \code{\link{neurovec_frame}}, \code{\link{latent_frame}}
#' @family model_estimation
#' @export
estimate_betas <- function(x, ...) UseMethod("estimate_betas")


# design_map generic moved to fmridesign package







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
