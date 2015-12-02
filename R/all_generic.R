
#' @export
parentTerms <- function(x) UseMethod("parentTerms")

#' @export
cells <- function(x, ...) UseMethod("cells")

#' conditions
#' 
#' return the set of condition labels associated with a model term
#' @param x the model term
#'@param ...
#' @export
conditions <- function(x, ...) UseMethod("conditions")

#' convolve
#' 
#' convolve a term \code{x} with a hemodynamic response over a \code{samplingFrame}
#' 
#' @export
#' @param x a event sequence
#' @param hrf a hemodynamic response function
#' @param samplingFrame the time series grid over whcih to sample the fucntion.
#' @param ...
convolve <- function(x, hrf, samplingFrame,...) UseMethod("convolve")

#' isContinuous
#' 
#' checks to see if a variable is continuous e.g. numeric/non-categorical
#' @export
isContinuous <- function(x) UseMethod("isContinuous")

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

#' designMatrix
#' 
#' construct a design matrix from the term
#' @param x the term
#' @export
designMatrix <- function(x, ...) UseMethod("designMatrix")

#' elements
#' 
#' return the ordered elements of a term/variable
#' 
#' @export
elements <- function(x, ...) UseMethod("elements")

#' evaluate
#' 
#' @export
evaluate <-  function(x, samplingGrid, ...) UseMethod("evaluate")

#' nbasis
#' @export
globalOnsets <-  function(x, onsets,...) UseMethod("globalOnsets")

#' @export
nbasis <-  function(x) UseMethod("nbasis")

#' onsets
#' @export
onsets <-  function(x) UseMethod("onsets")

#' durations
#' @export
durations <-  function(x) UseMethod("durations")

#' amplitudes
#' @export
amplitudes <-  function(x) UseMethod("amplitudes")

#' samples
#' @export
samples <-  function(x, ...) UseMethod("samples")