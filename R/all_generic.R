
#' longnames
#' 
#' @param x the object
#' @param ... extra args
#' @export
longnames <- function(x, ...) UseMethod("longnames")

#' shortnames
#' @param x the object
#' @param ... extra args
#' @export
shortnames <- function(x, ...) UseMethod("shortnames")



#' contrast_weights
#' 
#' @param x the object
#' @export
contrast_weights <- function(x, term, ...) UseMethod("contrast_weights")


#' parent_terms
#' 
#' @param x the object
#' @export
parent_terms <- function(x) UseMethod("parent_terms")


#' cells
#' 
#' return the experimental cells that are in a model term
#' 
#' @export
cells <- function(x, ...) UseMethod("cells")

#' conditions
#' 
#' return the set of condition labels associated with a model term
#' 
#' @param x the model term
#' @param ...
#' @export
conditions <- function(x, ...) UseMethod("conditions")

#' convolve
#' 
#' convolve a term \code{x} with a hemodynamic response over a \code{sampling_frame}
#' 
#' @export
#' @param x a event sequence
#' @param hrf a hemodynamic response function
#' @param sampling_frame the time series grid over whcih to sample the fucntion.
#' @param ...
convolve <- function(x, hrf, samplingFrame,...) UseMethod("convolve")

#' is_continuous
#' 
#' checks to see if a variable is continuous e.g. numeric/non-categorical
#' @export
is_continuous <- function(x) UseMethod("is_continuous")

#' is_categorical
#' 
#' checks to see if a variable is cateogircal e.g. factor-based
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


#' design_matrix
#' 
#' construct a design matrix from the term
#' 
#' @param x the term
#' @export
design_matrix <- function(x, ...) UseMethod("design_matrix")

#' elements
#' 
#' return the ordered elements of a term/variable
#' 
#' @export
elements <- function(x, ...) UseMethod("elements")


#' evaluate
#' 
#' evaluate a regression or hemodynamic response function
#' 
#' @param x
#' @param grid
#' @param ... extra args
#' 
#' @export
evaluate <-  function(x, grid, ...) UseMethod("evaluate")

#' global_onsets
#' 
#' @export
#' @param x the object
#' @param onsets the relative onset times of the events
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
#' @export
durations <-  function(x) UseMethod("durations")

#' amplitudes
#' 
#' return the amplitude vector
#' 
#' @export
amplitudes <-  function(x) UseMethod("amplitudes")

#' samples
#' 
#' extract samples
#' 
#' @export
samples <-  function(x, ...) UseMethod("samples")

#' split_by_block
#' 
#' split variables by block id
#' 
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
#' @export
blocklens  <-  function(x, ...) UseMethod("blocklens")

#' Fcontrasts
#' 
#' generate a set of Fcontrasts for a model term
#' 
#' @param x the term
#' @export
Fcontrasts <- function(x, ...) UseMethod("Fcontrasts")

#' estcon
#' 
#' @export
estcon <- function(x, fit, ...) UseMethod("estcon")

#' construct
#' 
#' construct a term given a an hrf spec and model specification
#' 
#' @param x the object
#' @param model_spec the model specification
#' @export
construct <- function(x, model_spec,...) UseMethod("construct")


