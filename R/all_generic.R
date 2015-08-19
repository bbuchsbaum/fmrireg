
#' @export
parentTerms <- function(x) UseMethod("parentTerms")

#' @export
cells <- function(x, ...) UseMethod("cells")

#' @export
conditions <- function(x, ...) UseMethod("conditions")

#' @export
convolve <- function(x, hrf, samplingFrame,...) UseMethod("convolve")

#' @export
isContinuous <- function(x) UseMethod("isContinuous")

#' @export
levels <- function(x) UseMethod("levels")

#' @export
columns <- function(x) UseMethod("columns")

#' @export
designMatrix <- function(x, ...) UseMethod("designMatrix")

#' @export
elements <- function(x, ...) UseMethod("elements")

#' @export
evaluate <-  function(x, samplingGrid, ...) UseMethod("evaluate")

#' @export
globalOnsets <-  function(x, onsets,...) UseMethod("globalOnsets")

#' @export
nbasis <-  function(x) UseMethod("nbasis")

#' @export
onsets <-  function(x) UseMethod("onsets")

#' @export
durations <-  function(x) UseMethod("durations")

#' @export
amplitudes <-  function(x) UseMethod("amplitudes")

#' @export
samples <-  function(x, ...) UseMethod("samples")