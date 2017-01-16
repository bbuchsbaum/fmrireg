

#' @param event_model
#' @param baseline_model
#' @param event_table
#' @param basis
#' @param durations
#' @param sampling_frame
#' @export
fmri_model <- function(event_model, baseline_model) {
  assert_that(inherits(event_model, "event_model"))
  assert_that(inherits(baseline_model, "baseline_model"))
  
  fmodel <- list(event_model=event_model, baseline_model=baseline_model)
  class(fmodel) <- "fmri_model"
  fmodel
}


#' @importFrom tibble as_tibble
design_matrix.fmri_model <- function(x) {
  tibble::as_tibble(cbind(design_matrix(x$event_model), design_matrix(x$baseline_model)))
}



#' @export
terms.fmri_model <- function(x) {
  terms(x$event_model)
}

#' @export
terms.event_model <- function(x) {
  x$terms
}



#' @export
conditions.fmri_model <- function(x) {
  unlist(lapply(terms(x), function(t) conditions(t)), use.names=FALSE)
}

#' @export
print.fmri_model <- function(object) {
  cat("fmri_model", "\n")
  cat(" ", "Event Model:  ", Reduce(paste, deparse(object$model_spec$formula)), "\n")
  cat(" ", "Baseline Model:  ", Reduce(paste, deparse(object$model_spec$baseline_formula)), "\n")
  cat(" ", "Num Terms", length(terms(object)), "\n")
  cat(" ", "Num Events: ", nrow(object$model_spec$event_table), "\n")
  cat(" ", "Num Columns: ", length(conditions(object)), "\n")
  cat(" ", "Num Blocks: ", length(object$model_spec$blocklens), "\n")
  cat(" ", "Length of Blocks: ", paste(object$model_spec$blocklens, collapse=", "), "\n")
  for (i in 1:length(terms(object))) {
    t <- terms(object)[[i]]
    cat("Term:", i, " ")
    print(t)
    cat("\n")
  }
  
}


