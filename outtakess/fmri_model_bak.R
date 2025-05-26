#' Create an fMRI Model
#'
#' This function creates an fMRI model consisting of an event model and a baseline model.
#'
#' @param formula The model formula for experimental events.
#' @param block The model formula for block structure.
#' @param baseline_model (Optional) A \code{baseline_model} object. Default is \code{NULL}.
#' @param dataset An \code{fmri_dataset} object containing the time-series data.
#' @param drop_empty Logical. Whether to remove factor levels with zero size. Default is \code{TRUE}.
#' @param durations A vector of event durations. Default is \code{0}.
#' @return An \code{fmri_model} object.
#' @keywords internal
#' @noRd
create_fmri_model <- function(formula, block, baseline_model = NULL, dataset, drop_empty = TRUE, durations = 0) {
  assert_that(is.formula(formula), msg = "'formula' must be a formula")
  assert_that(is.formula(block), msg = "'block' must be a formula")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset'")
  assert_that(is.numeric(durations), msg = "'durations' must be numeric")
  
  if (is.null(baseline_model)) {
    baseline_model <- baseline_model(
      basis = "bs",
      degree = max(ceiling(median(dataset$sampling_frame$blocklens) / 100), 3),
      sframe = dataset$sampling_frame
    )
  } else {
    assert_that(inherits(baseline_model, "baseline_model"),
                msg = "'baseline_model' must have class 'baseline_model'")
  }
  
  ev_model <- event_model(
    x = formula,
    block = block,
    data = dataset$event_table,
    sampling_frame = dataset$sampling_frame,
    drop_empty = drop_empty,
    durations = durations
  )
  
  fmri_model(ev_model, baseline_model)
}


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

#' @noRd
prediction_matrix <- function(x) {
  #model$event_model
  ## --> tells us whether this term needs expanding. nbasis(model$event_model$terms[[1]]$evterm)
  ## iterate over every term. 
  ## get factors and discover continuos variables.
  ## get the range of each continuous variable, get the levels of each factor.
  ## cross all elements with expand.grid, with basis functions if necessary.
  ## do not include redundant variables.
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
#' @noRd
design_env.fmri_model <- function(x, blockid=NULL) {
  stop("not implemented")
  
}


#' @export
terms.fmri_model <- function(x,...) {
  c(terms(x$event_model), terms(x$baseline_model))
}

#' @export
#' @autoglobal
cells.fmri_model <- function(x, ...) {
  c1 <- cells(x$event_model) %>% dplyr::mutate(type="event")
  c2 <- cells(x$baseline_model) %>% dplyr::mutate(type="baseline")
  rbind(c1,c2) %>% dplyr::relocate(index, type)
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


#' @export
plot.fmri_model <- function(x,...) {
  with_package("cowplot")
  p1 <- plot(x$event_model) + ggplot2::ggtitle("Event Model")
  p2 <- plot(x$baseline_model) + ggplot2::ggtitle("Baseline Model")
  cowplot::plot_grid(p1,p2,nrow=2, align="h")
}



#' @export
print.fmri_model <- function(x, ...) {
  # Header with fancy border
  cat("\n╔══════════════════════════════════════════╗")
  cat("\n║             fMRI Model                   ║")
  cat("\n╠══════════════════════════════════════════╣")
  
  # Event Model Section
  cat("\n║ Event Model                              ║")
  cat("\n╟──────────────────────────────────────────╢")
  cat("\n║ Formula:", crayon::cyan(Reduce(paste, deparse(x$event_model$model_spec$formula))))
  
  # Event Model Summary
  cat("\n║ Summary:")
  cat("\n║   • Terms:", crayon::yellow(length(terms(x$event_model))))
  cat("\n║   • Events:", crayon::yellow(nrow(x$event_model$model_spec$event_table)))
  cat("\n║   • Design Columns:", crayon::yellow(length(conditions(x$event_model))))
  cat("\n║   • Blocks:", crayon::yellow(length(unique(x$event_model$blockids))))
  
  # Baseline Model Section (if present)
  if (!is.null(x$baseline_model)) {
    cat("\n╟──────────────────────────────────────────╢")
    cat("\n║ Baseline Model                           ║")
    cat("\n║ Components:")
    
    # Drift term info
    if (!is.null(x$baseline_model$drift_term)) {
      drift_name <- x$baseline_model$drift_term$varname
      basis_type <- x$baseline_model$drift_spec$basis
      degree <- x$baseline_model$drift_spec$degree
      drift_cols <- ncol(design_matrix(x$baseline_model$drift_term))
      cat("\n║   • Drift:", crayon::magenta(drift_name))
      cat("\n║     - Type:", crayon::blue(basis_type))
      cat("\n║     - Degree:", crayon::blue(degree))
      cat("\n║     - Columns:", crayon::yellow(drift_cols))
    }
    
    # Block term info
    if (!is.null(x$baseline_model$block_term)) {
      const_cols <- ncol(design_matrix(x$baseline_model$block_term))
      cat("\n║   • Block Terms:", crayon::yellow(const_cols), "columns")
    }
    
    # Nuisance term info
    if (!is.null(x$baseline_model$nuisance_term)) {
      nuis_cols <- ncol(design_matrix(x$baseline_model$nuisance_term))
      cat("\n║   • Nuisance Terms:", crayon::yellow(nuis_cols), "columns")
    }
  }
  
  # Total Model Summary
  cat("\n╟──────────────────────────────────────────╢")
  cat("\n║ Total Model                              ║")
  total_cols <- ncol(design_matrix(x))
  cat("\n║   • Total Design Columns:", crayon::yellow(total_cols))
  
  # Footer
  cat("\n╚══════════════════════════════════════════╝\n")
}


