
get_methods <- function(obj) {
  unique(purrr::map_chr(class(obj), ~ methods(class= . )))
}


as_vectors <- function(x) { function(x, ...) UseMethod("as_vectors") }

setGeneric("as_vectors") 

#' get_formula
#' 
#' @param x the object
#' @param ... extra args
get_formula <- function(x, ...) UseMethod("get_formula")

#' term_matrices
#' 
#' @param x the object
#' @param ... extra args
term_matrices <- function(x, ...) UseMethod("term_matrices")




#' longnames
#' 
#' @param x the object
#' @param ... extra args
#' @export
longnames <- function(x, ...) UseMethod("longnames")

#' shortnames
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
#' @export
design_env <- function(x, ...) UseMethod("design_env")



#' contrast_weights
#' 
#' @param x the object
#' @param ... extra args
#' @export
contrast_weights <- function(x, ...) UseMethod("contrast_weights")


#' parent_terms
#' 
#' @param x the object
#' @export
parent_terms <- function(x) UseMethod("parent_terms")

#' term_names
#' @param x the object to extra term names from
#' @export
term_names <- function(x) UseMethod("term_names")

#' cells
#' 
#' return the experimental cells that are in a model term
#' 
#' @param x the object
#' @param ... extra args
#' @export
cells <- function(x, ...) UseMethod("cells")

#' conditions
#' 
#' return the set of condition labels associated with a model term
#' 
#' @param x the model term
#' @param ... extra args
#' @export
conditions <- function(x, ...) UseMethod("conditions")

#' convolve
#' 
#' convolve a term \code{x} with a hemodynamic response over a \code{sampling_frame}
#' 
#' @export
#' @param x an event sequence
#' @param hrf a hemodynamic response function
#' @param sampling_frame the time series grid over whcih to sample the fucntion.
#' @param ... extra args
convolve <- function(x, hrf, sampling_frame,...) UseMethod("convolve")

#' is_continuous
#' 
#' checks to see if a variable is continuous e.g. numeric/non-categorical
#' 
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
#' @param x the term
#' @export
columns <- function(x) UseMethod("columns")


#' event_table
#' 
#' extract event_table as \code{data.frame}
#' 
#' @param x the term
#' @export
event_table <- function(x) UseMethod("event_table")


#' event_terms
#' 
#' extract list of 'event_term's from model object
#' 
#' @param x the model object
#' @export
event_terms <- function(x) UseMethod("event_terms")

#' baseline_terms
#' 
#' extract list of 'baseline_term's from model object
#' 
#' @param x the model object
#' @export
baseline_terms <- function(x) UseMethod("baseline_terms")


#' term_indices
#' 
#' the term indices for the associated design matrix
#' 
#' @param x the model/term object
#' @param ... additional arguments
#' @export
term_indices <- function(x,...) UseMethod("term_indices")


#' run a command
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


#' evaluate
#' 
#' evaluate a regression or hemodynamic response function
#' 
#' @param x the object to evaluate
#' @param grid the sampling grid
#' @param ... extra args
#' 
#' @export
evaluate <-  function(x, grid, ...) UseMethod("evaluate")

#' global_onsets
#' 
#' @export
#' @param x the object
#' @param onsets the relative onset times of the events
#' @param ... extra args
global_onsets <-  function(x, onsets,...) UseMethod("global_onsets")

#' nbasis
#' 
#' return number of basis functions associated with hrf.
#' 
#' @export
#' @param x the object
nbasis <-  function(x) UseMethod("nbasis")


#' data_chunks
#' 
#' return the onset vector
#' @param x the dataset
#' @param nchunks the numbe rof data chunks
#' @param ... extra args
#' @return an \code{iterator} returning on data chunk per iteration
#' @export
data_chunks <- function(x, nchunks, ...) UseMethod("data_chunks")


#' onsets
#' 
#' return the onset vector
#' @param x the object
#' @export
onsets <-  function(x) UseMethod("onsets")

#' durations
#' 
#' return the durations vector
#' 
#' @param x the object
#' @export
durations <-  function(x) UseMethod("durations")

#' amplitudes
#' 
#' return the amplitude vector
#' 
#' @param x the object
#' @export
amplitudes <-  function(x) UseMethod("amplitudes")

#' samples
#' 
#' extract samples
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

#' blockids
#' 
#' extract block indices
#' 
#' @param x the object
#' @export
blockids  <-  function(x) UseMethod("blockids")

#' blocklens
#' 
#' extract block lengths
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


#' p_value
#' 
#' extract p-values from fitted model
#' 
#' @param x the object
#' @param ... extra arguments
p_value <- function(x, ...) UseMethod("p_value")


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
#' @export
build_afni_stims <- function(x, ...) UseMethod("build_afni_stims")

#' split an onset vector into a list
#' 
#' @param x the object associated with the onset vector
#' @param ... extra args
#' @export
split_onsets <- function(x, ...) UseMethod("split_onsets")
  

