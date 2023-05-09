
#' Construct an fMRI regression model
#'
#' This function constructs an fMRI regression model consisting of an event model
#' and a baseline model. The resulting model can be used for the analysis of fMRI data.
#'
#' @param event_model An object of class "event_model" representing the event-related part of the fMRI regression model.
#' @param baseline_model An object of class "baseline_model" representing the baseline-related part of the fMRI regression model.
#' @return An object of class "fmri_model" containing the event and baseline models.
#' @export
#' @seealso event_model, baseline_model
fmri_model <- function(event_model, baseline_model) {
  assert_that(inherits(event_model, "event_model"))
  assert_that(inherits(baseline_model, "baseline_model"))
  
  fmodel <- list(event_model=event_model, baseline_model=baseline_model)
  class(fmodel) <- "fmri_model"
  fmodel
}


#' @importFrom tibble as_tibble
#' @param blockid the block id to extract
#' @export
#' @rdname design_matrix
design_matrix.fmri_model <- function(x, blockid=NULL, ...) {
  suppressMessages(tibble::as_tibble(cbind(design_matrix(x$event_model, blockid), 
                                           design_matrix(x$baseline_model, blockid)),.name_repair="check_unique"))
}

#' @importFrom tibble as_tibble
#' @keywords internal
design_env.fmri_model <- function(x, blockid=NULL) {
  stop("not implemented")
  
}


#' @export
terms.fmri_model <- function(x,...) {
  c(terms(x$event_model), terms(x$baseline_model))
}

#' @export
blocklens.fmri_model <- function(x,...) {
  blocklens(x$event_model)
}


#' @export
event_terms.fmri_model <- function(x) {
  terms(x$event_model)
}

#' @export
baseline_terms.fmri_model <- function(x) {
  terms(x$baseline_model)
}

#' @export
contrast_weights.fmri_model <- function(x, ...) {
  contrast_weights.event_model(x$event_model,...)
}


#' @export
conditions.fmri_model <- function(x, ...) {
  unlist(lapply(terms(x), function(t) conditions(t)), use.names=FALSE)
}

#' @export
conditions.baseline_model <- function(x, ...) {
  unlist(lapply(terms(x), function(t) conditions(t)), use.names=FALSE)
}


#' @importFrom cowplot plot_grid
#' @export
plot.fmri_model <- function(x,...) {
  p1 <- plot(x$event_model) + ggplot2::ggtitle("Event Model")
  p2 <- plot(x$baseline_model) + ggplot2::ggtitle("Baseline Model")
  cowplot::plot_grid(p1,p2,nrow=2, align="h")
}



#' @export
print.fmri_model <- function(x,...) {
  cat("fmri_model", "\n")
  cat(" ", "Event Model:  ", Reduce(paste, deparse(x$event_model$model_spec$formula)), "\n")
  cat(" ", "Baseline Model:  ", x$baseline_model$drift_term$varname, "\n")
  cat(" ", "Num Terms", length(terms(x)), "\n")
  cat(" ", "Num Events: ", nrow(x$model_spec$event_table), "\n")
  cat(" ", "Num Columns: ", length(conditions(x)), "\n")
  cat(" ", "Num Blocks: ", length(x$event_model$model_spec$sampling_frame$blocklens), "\n")
  cat(" ", "Length of Blocks: ", paste(x$event_model$model_spec$sampling_frame$blocklens, collapse=", "), "\n")
  for (i in 1:length(terms(x))) {
    cat("\n")
    t <- terms(x)[[i]]
    cat("Term:", i, " ")
    print(t)
  }
  
}


