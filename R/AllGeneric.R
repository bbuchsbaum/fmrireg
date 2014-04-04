

#' Generic function to extract onsets
#' @param x the object with onsets slot
#' @export 
#' @rdname onsets-methods
setGeneric("onsets",  function(x) standardGeneric("onsets"))

#' Generic function to extract amplitudes
#' @param x the object with amplitude slot
#' @export 
#' @rdname amplitude-methods
setGeneric("amplitude", function(x) standardGeneric("amplitude"))

#' Generic function to extract durations
#' @param x the object with duration slot
#' @export 
#' @rdname duration-methods
setGeneric("duration",  function(x) standardGeneric("duration"))


#' Generic function to evaluate a regressor over a discrate sampling grid
#' @param x the object to evaluate
#' @param grid the sampling grid
#' @param ... additional arguments
#' @export 
#' @rdname evaluate-methods
setGeneric("evaluate",  function(x, grid, ...) standardGeneric("evaluate"))