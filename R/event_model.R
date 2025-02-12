
###############################################################################
## event_model.R
##
## This file constructs event models for fMRI analysis. Two primary interfaces 
## are provided: one for when the model is specified as a list of HRF specifications 
## (event_model.list) and one when specified as a formula (event_model.formula). 
## Additional helper functions are provided for constructing the model specification,
## extracting terms, parsing variables, building design matrices, computing contrasts,
## and printing/plotting the event model.
##
## Changes made:
##  - Replaced lazyeval calls with rlang equivalents.
##  - Added additional error checking and clarifying inline comments.
##  - Removed extraneous debugging prints.
###############################################################################


.pre_eval_hrf_calls <- function(expr, fenv) {
  # If it's NOT a call, just return as-is
  if (!rlang::is_call(expr)) {
    return(expr)
  }
  
  # If the function name is 'hrf' or 'trialwise', evaluate the entire call 
  # in the formula environment. Then return that result (which should be an hrfspec).
  fun_name <- rlang::call_name(expr)
  if (fun_name %in% c("hrf", "trialwise", "afni_hrf", "afni_trialwise")) {
    # Evaluate the call in the formula’s environment
    evald <- rlang::eval_tidy(expr, env=fenv)
    # 'evald' should be an hrfspec (or related).
    return(evald)
  }
  
  # Otherwise, we recursively walk arguments of the call
  # so if the user has nested calls, we can catch them too.
  expr_args <- as.list(expr)
  # first element is the function name
  for (i in seq_along(expr_args)[-1]) {
    expr_args[[i]] <- .pre_eval_hrf_calls(expr_args[[i]], fenv)
  }
  # Rebuild the call
  as.call(expr_args)
}

## ============================================================================
## Section 1: Event Model Construction Functions
## ============================================================================
#' event_model.list
#'
#' Construct an event model from a list of HRF specifications.
#'
#' @description
#' Constructs an event model from a list of hemodynamic response function (HRF)
#' specifications (`x`), along with additional parameters (data, block, sampling frame,
#' durations, etc.). All elements in `x` must be of type `hrfspec`.
#'
#' @param x A list of HRF specifications.
#' @param data The dataset (a data.frame) containing event information.
#' @param block The block variable; either a formula or a vector of block values.
#' @param sampling_frame The time series grid over which to sample the function.
#' @param drop_empty Logical indicating whether to drop empty events. Default is TRUE.
#' @param durations A numeric vector of event durations. Default is 0 for all events.
#' @param precision Numeric value indicating the precision for HRF sampling. Default is 0.3.
#' @param ... Additional arguments.
#' @return An event model object.
#' @examples 
#' # Example for event_model.list:
#' # Assume HRF_SPMG1 is a predefined hrfspec object.
#' df <- data.frame(onset = seq(1, 100, by = 10),
#'                  run = rep(1:2, each = 5))
#' # Create a list of HRF specifications:
#' hrf_list <- list(HRF_SPMG1)
#' # Create a sampling frame (assume sampling_frame() is defined)
#' sframe <- sampling_frame(blocklens = c(50, 50), TR = 2)
#' # Construct the event model
#' ev_model <- event_model.list(x = hrf_list, data = df, block = ~ run, sampling_frame = sframe)
#' @export
event_model.list <- function(x, data, block, sampling_frame, drop_empty = TRUE, durations = 0, 
                             precision = 0.3, ...) {
  assert_that(all(sapply(x, function(a) inherits(a, "hrfspec"))),
              msg = "event_model.list: all `x` elements must be of type `hrfspec`")
  
  if (rlang::is_formula(block)) {
    ## TODO: check that the block variable exists in data and warn if onsets are off.
    block_rhs <- rlang::f_rhs(block)
    blockvals <- rlang::eval_tidy(block_rhs, data)
  } else {
    blockvals <- block
  }
  
  if (is.factor(blockvals)) {
    blockvals <- as.integer(as.character(blockvals))
  }
  
  assert_that(is.increasing(blockvals), msg = "'blockvals' must consist of strictly increasing integers")
  assert_that(length(blockvals) == nrow(data))
  
  blocks <- unique(blockvals)
  ranked_blocks <- rank(blocks)
  blockids <- ranked_blocks[match(blockvals, blocks)]
  
  if (missing(durations)) {
    ## assume zero-duration impulse for all events
    durations <- rep(0, nrow(data))
  }
  
  form <- paste("~", paste(sapply(x, function(z) z$label), collapse = "+"))
  
  model_spec <- list(formula = form,
                     event_table = data, 
                     onsets = x[[1]]$onsets, 
                     event_spec = x, 
                     blockvals = blockvals,
                     blockids = blockids, 
                     durations = durations, 
                     sampling_frame = sampling_frame,
                     drop_empty = drop_empty,
                     precision = precision)
  
  class(model_spec) <- c("event_model_spec", "list")
  fmodel <- construct_model(model_spec)
  fmodel
}


#' event_model.formula
#'
#' Construct an event model from a formula and data.
#'
#' @description
#' Constructs an event model using a formula `x`, a dataset `data`, a block variable,
#' and additional parameters (sampling frame, durations, precision, etc.). The right-hand
#' side of the formula must consist of HRF terms (of class "hrfspec").
#'
#' @param x A formula specifying the event model.
#' @param data A data.frame containing the event information.
#' @param block The block variable; either a formula or a vector of block values.
#' @param sampling_frame The time series grid over which to sample the function.
#' @param drop_empty Logical indicating whether to drop empty events. Default is TRUE.
#' @param durations A numeric vector of event durations. Default is 0 for all events.
#' @param precision Numeric value indicating the precision of the model. Default is 0.3.
#' @param ... Additional arguments.
#' @return An event model object.
#' @examples 
#' # Example for event_model.formula:
#' df <- data.frame(onset = seq(1, 100, by = 10),
#'                  run = rep(1:2, each = 5),
#'                  x = rnorm(10),
#'                  y = rnorm(10))
#' # Create a sampling frame (assume sampling_frame() is defined)
#' sframe <- sampling_frame(blocklens = c(50, 50), TR = 2)
#' # Construct an event model using a formula.
#' # Here the left-hand side represents onsets and the right-hand side contains HRF terms.
#' ev_model <- event_model(x = onset ~ hrf(x) + hrf(y), data = df, block = ~ run, sampling_frame = sframe)
#' @export
#' @importFrom purrr map_lgl
event_model.formula <- function(x, data, block, sampling_frame, drop_empty = TRUE, durations = 0, precision = 0.3, ...) {
  formula <- x
  #environment(formula) <- environment()
  stopifnot(inherits(formula, "formula"))
  assert_that(inherits(data, "data.frame"), msg = "`data` must be a data.frame")
  
  if (rlang::is_formula(block)) {
    ## TODO: check that block exists in data and warn if onsets are way off.
    block_rhs <- rlang::f_rhs(block)
    blockvals <- rlang::eval_tidy(block_rhs, data)
  } else {
    blockvals <- block
  }
  
  if (is.factor(blockvals)) {
    blockvals <- as.integer(as.character(blockvals))
  }
  
  assert_that(is.increasing(blockvals), msg = "'blockvals' must consist of strictly increasing integers")
  assert_that(length(blockvals) == nrow(data))
  
  blocks <- unique(blockvals)
  ranked_blocks <- rank(blocks)
  blockids <- ranked_blocks[match(blockvals, blocks)]
  
  if (missing(durations)) {
    durations <- rep(0, nrow(data))
  }
  
  formspec <- function(formula, table) {
    
    vterms <- extract_terms(formula, table)
    resp <- attr(vterms, "response")
    variables <- extract_variables(formula, table, vterms)
    
    lhs <- if (resp != 0) variables[[resp]] else NULL
    rhs <- variables[(resp + 1):length(variables)]
    
    ret <- list(lhs = lhs, rhs = rhs)
    class(ret) <- "formula_extraction"
    ret
  }
  
  #browser()
  
  event_spec <- formspec(formula, data)
  assert_that(all(map_lgl(event_spec$rhs, inherits, "hrfspec")),
              msg = "all terms on right-hand side must be 'hrfspec' terms")
  
  model_spec <- list(
    formula = formula, 
    event_table = data, 
    onsets = event_spec$lhs, 
    event_spec = event_spec$rhs, 
    blockvals = blockvals,
    blockids = blockids, 
    durations = durations, 
    sampling_frame = sampling_frame,
    drop_empty = drop_empty,
    contrasts = contrasts,
    precision = precision
  )
  
  class(model_spec) <- c("event_model_spec", "list")
  fmodel <- try(construct_model(model_spec))
  fmodel
}

#' @export
blocklens.event_model <- function(x, ...) {
  blocklens(x$sampling_frame)
}

#' @export
blockids.event_model <- function(x) {
  x$blockids
}


## ============================================================================
## Section 2: Model Specification and Construction Helpers
## ============================================================================

#' Construct an Event Model.
#'
#' Given a model_spec list (constructed by event_model.list or event_model.formula),
#' this function constructs the final event model.
#'
#' @param x A model_spec list.
#' @return An event_model object.
#' @keywords internal
construct_model <- function(x) {
  # Ensure onsets is numeric.
  assert_that(is.numeric(x$onsets))
  
  # Get term names from event_spec (each should have a "name" element).
  term_names <- sapply(x$event_spec, "[[", "name")
  term_names <- .sanitizeName(term_names)
  
  # If there are duplicates, append an index.
  dups <- sum(duplicated(term_names)) > 0
  if (dups) {
    dup_ids <- ave(term_names, term_names, FUN = seq_along)
    term_names <- paste0(term_names, "_", dup_ids)
  }
  
  # Construct each term by calling construct() on each hrfspec.
  terms <- lapply(x$event_spec, function(m) construct(m, x))
  names(terms) <- term_names
  
  term_lens <- sapply(lapply(terms, conditions), length)
  spans <- c(0, cumsum(term_lens))
  
  term_indices <- lapply(1:(length(spans) - 1), function(i) {
    seq(spans[i] + 1, spans[i + 1])
  })
  names(term_indices) <- term_names
  
  ret <- list(
    term_indices = term_indices,
    terms = terms,
    blockids = x$blockids,
    sampling_frame = x$sampling_frame,
    contrasts = x$contrasts,
    model_spec = x
  )
  class(ret) <- c("event_model", "list")
  ret
}

#' @keywords internal
#' @noRd
extract_terms <- function(formula, data) {
  if (!inherits(formula, "terms")) {
    terms(formula, data = data)
  } else {
    formula
  }
}

#' @keywords internal
#' @noRd
is_parametric_basis <- function(obj) { inherits(obj, "ParametricBasis") }

#' @keywords internal
#' @noRd
extract_variables <- function(form, data, .terms = NULL) {
  if (is.null(.terms)) {
    .terms <- extract_terms(form, data)
  }
  
  env <- environment(.terms)
  varnames <- sapply(attr(.terms, "variables"), deparse, width.cutoff = 500)[-1]
  # Remove debugging print.
  #print(varnames)
  #browser()
  #variables <- try(eval(attr(.terms, "variables"), data))
  #browser()
  variables <- eval(attr(.terms, "variables"), envir=data, enclos=environment(form))
  
  names(variables) <- varnames
  variables
}

#' @keywords internal
#' @noRd
parse_term <- function(vars, ttype) {
  dim <- length(vars) # number of variables
  term <- deparse(vars[[1]], backtick = TRUE) # first covariate
  if (dim > 1) {
    for (i in 2:dim) term[i] <- deparse(vars[[i]], backtick = TRUE)
  }
  for (i in 1:dim) term[i] <- attr(terms(reformulate(term[i])), "term.labels")
  
  full.call <- paste(ttype, "(", term[1], sep = "")
  if (dim > 1) {
    for (i in 2:dim) full.call <- paste(full.call, ",", term[i], sep = "")
  }
  label <- paste(full.call, ")", sep = "")
  list(term = term, label = label)
}

## ============================================================================
## Section 3: Design Matrix and Contrast Functions for the Event Model
## ============================================================================

#' @export
#' @importFrom tibble as_tibble
design_matrix.event_model_spec <- function(x, ...) {
  termlist <- lapply(x$varspec, function(m) construct(m, x))
  ret <- lapply(termlist, design_matrix)
  vnames <- unlist(lapply(ret, colnames))
  dmat <- suppressMessages(tibble::as_tibble(do.call(cbind, ret), .name_repair = "check_unique"))
  names(dmat) <- vnames
  dmat
}

#' @importFrom tibble as_tibble
#' @export
#' @rdname design_matrix
design_matrix.event_model <- function(x, blockid = NULL, ...) {
  ret <- lapply(x$terms, design_matrix, blockid)
  vnames <- unlist(lapply(ret, names))
  dmat <- suppressMessages(tibble::as_tibble(do.call(cbind, ret), .name_repair = "check_unique"))
  names(dmat) <- vnames
  dmat
}

#' @export
terms.event_model <- function(x, ...) {
  x$terms
}

#' @export
conditions.event_model <- function(x, ...) {
  unlist(lapply(terms(x), function(t) conditions(t)), use.names = FALSE)
}

#' @export
#' @rdname contrast_weights
contrast_weights.convolved_term <- function(x, ...) {
  lapply(x$contrasts, function(cspec) {
    if (!is.null(cspec))
      contrast_weights(cspec, x)
  })
}

#' @export
#' @rdname Fcontrasts
Fcontrasts.convolved_term <- function(x, ...) {
  Fcontrasts(x$evterm)
}

#' @export
#' @rdname contrast_weights
contrast_weights.fmri_model <- function(x, ...) { 
  contrast_weights(x$event_model) 
}

#' @export
term_names.event_model <- function(x) {
  xt <- terms(x)
  unlist(lapply(xt, function(term) term$varname))
}

#' @export
#' @rdname contrast_weights
contrast_weights.event_model <- function(x, ...) {
  tnames <- term_names(x)
  tind <- x$term_indices
  ncond <- length(conditions(x))
  ret <- lapply(seq_along(terms(x)), function(i) {
    cwlist <- contrast_weights(terms(x)[[i]])
    len <- length(conditions(terms(x)[[i]]))
    
    if (!is.null(cwlist) && length(cwlist) > 0) {
      ret <- lapply(cwlist, function(cw) {
        out <- matrix(0, ncond, ncol(cw$weights))
        out[tind[[i]], ] <- cw$weights
        attr(cw, "term_indices") <- as.vector(tind[[i]])
        attr(cw, "offset_weights") <- out
        cw$offset_weights <- out
        cw
      })
      
      cnames <- sapply(cwlist, function(cw) cw$name)
      names(ret) <- cnames
      ret
    }
  })
  names(ret) <- tnames
  ret
}

#' @export
Fcontrasts.event_model <- function(x, ...) {
  tind <- x$term_indices
  tnames <- names(terms(x))
  ret <- unlist(lapply(seq_along(terms(x)), function(i) {
    eterm <- terms(x)[[i]]$evterm
    len <- length(conditions(eterm))
    cwlist <- Fcontrasts(terms(x)[[i]])
    if (!is.null(cwlist)) {
      ret <- lapply(cwlist, function(cw) {
        out <- matrix(0, len, ncol(cw))
        ti <- tind[[i]]
        out <- as.matrix(cw)
        attr(out, "term_indices") <- as.vector(ti)
        row.names(out) <- row.names(cw)
        out
      })
      cnames <- names(cwlist)
      prefix <- tnames[i]
      names(ret) <- paste0(prefix, "#", cnames)
      ret
    }
  }), recursive = FALSE)
  ret
}

#' @export
#' @rdname design_matrix
design_matrix.convolved_term <- function(x, blockid = NULL, ...) {
  if (is.null(blockid)) {
    x$design_matrix
  } else {
    keep <- blockids(x$sampling_frame) %in% blockid
    x$design_matrix[keep, ]
  }
}

#' @export
design_matrix.afni_hrf_convolved_term <- function(x, blockid = NULL, ...) {
  stop("afni_hrf_convolved_term delegates design matrix construction to AFNI")
}

#' matrix_term
#'
#' Create a matrix_term object.
#'
#' @description
#' Creates a matrix_term object, which is a set of regression variables stored as
#' a numeric matrix.
#'
#' @param varname The name of the variable.
#' @param mat The matrix of values.
#' @return A matrix_term object.
#' @importFrom tibble as_tibble
#' @examples 
#' mat <- matrix(rnorm(100 * 10), 100, 10)
#' mterm <- matrix_term("mterm", mat)
#' @export
matrix_term <- function(varname, mat) {
  stopifnot(is.matrix(mat))
  ret <- list(varname = varname, design_matrix = suppressMessages(tibble::as_tibble(mat, .name_repair = "minimal")))
  class(ret) <- c("matrix_term", "fmri_term", "list")
  ret
}

#' @export
#' @rdname design_matrix
design_matrix.matrix_term <- function(x, ...) {
  if (is.null(names(x$design_matrix))) {
    cnames <- paste0(x$varname, "_", 1:ncol(x$design_matrix))
    names(x$design_matrix) <- cnames
  }
  x$design_matrix
}

#' @export
#' @rdname event_table
event_table.convolved_term <- function(x) event_table(x$evterm)

#' @export
#' @rdname nbasis
nbasis.convolved_term <- function(x) nbasis(x$hrf)

#' @export
cells.event_model <- function(x, ...) {
  ret <- lapply(x$terms, function(z) cells(z, ...))
  tnames <- names(ret)
  out <- do.call(rbind, lapply(seq_along(tnames), function(i) {
    dplyr::tibble(term = tnames[i], level = ret[[i]][, 1], basis = ret[[i]][, 2])
  })) %>% dplyr::mutate(index = 1:dplyr::n())
}

#' @export
#' @rdname longnames
longnames.convolved_term <- function(x, ...) {
  term.cells <- cells(x)
  apply(as.matrix(sapply(1:ncol(term.cells), function(i) {
    paste0(names(term.cells)[i], "#", term.cells[[i]], sep = "")
  })), 1, paste, collapse = ":")
}

#' @export
#' @rdname longnames
longnames.afni_hrf_convolved_term <- function(x, ...) {
  term.cells <- cells(x, exclude_basis = TRUE)
  apply(as.matrix(sapply(1:ncol(term.cells), function(i) {
    paste0(names(term.cells)[i], "#", term.cells[[i]], sep = "")
  })), 1, paste, collapse = ":")
}

#' @export
#' @rdname longnames
longnames.event_model <- function(x, ...) {
  unlist(lapply(terms(x), longnames))
}

#' @export
#' @rdname longnames
longnames.event_term <- function(x, ...) {
  term.cells <- cells(x)
  apply(as.matrix(sapply(1:ncol(term.cells), function(i) {
    paste0(names(term.cells)[i], "#", term.cells[[i]], sep = "")
  })), 1, paste, collapse = ":")
}

#' @export
#' @rdname shortnames
shortnames.event_model <- function(x, ...) {
  unlist(lapply(terms(x), shortnames))
}

#' @export
#' @rdname shortnames
shortnames.convolved_term <- function(x, ...) {
  term.cells <- cells(x)
  apply(as.matrix(sapply(1:ncol(term.cells), function(i) {
    term.cells[[i]]
  })), 1, paste, collapse = ":")
}


#' @export
#' @rdname shortnames
shortnames.event_term <- function(x, ...) {
  term.cells <- cells(x)
  apply(as.matrix(sapply(1:ncol(term.cells), function(i) {
    term.cells[[i]]
  })), 1, paste, collapse = ":")
}

#' @export
#' @rdname shortnames
shortnames.matrix_term <- function(x, ...) {
  colnames(x$design_matrix)
}

#' @export
#' @rdname longnames
longnames.matrix_term <- function(x, ...) {
  paste0(x$name, "#", colnames(design_matrix(x)))
}

#' @export
print.event_model <- function(x, ...) {
  cat("\n╔══════════════════════════════════════════╗")
  cat("\n║           fMRI Event Model               ║")
  cat("\n╚══════════════════════════════════════════╝\n")
  
  cat("\n Model Formula:\n")
  cat("  ", crayon::cyan(Reduce(paste, deparse(x$model_spec$formula))), "\n")
  
  cat("\n📈 Model Summary:\n")
  cat("  • Number of Terms:", crayon::yellow(length(terms(x))), "\n")
  cat("  • Total Events:", crayon::yellow(nrow(x$model_spec$event_table)), "\n")
  cat("  • Design Matrix Columns:", crayon::yellow(length(conditions(x))), "\n")
  cat("  • Number of Blocks:", crayon::yellow(length(unique(x$blockids))), "\n")
  
  if (nrow(x$model_spec$event_table) > 0) {
    cat("\n📋 Event Table Preview:\n")
    cat("  • Variables:", crayon::green(paste(names(x$model_spec$event_table), collapse = ", ")), "\n")
    if (nrow(x$model_spec$event_table) > 3) {
      print(head(x$model_spec$event_table, 3))
      cat("  ... (", nrow(x$model_spec$event_table) - 3, " more rows )\n")
    } else {
      print(x$model_spec$event_table)
    }
  }
  
  cat("\n🔍 Model Terms:\n")
  for (i in seq_along(terms(x))) {
    term <- terms(x)[[i]]
    cat("\n  Term", crayon::blue(i), ":", crayon::magenta(names(terms(x))[i]), "\n")
    
    term_conditions <- conditions(term)
    if (length(term_conditions) > 0) {
      cat("    • Conditions:", crayon::green(paste(term_conditions, collapse = ", ")), "\n")
    }
    
    if (!is.null(term$hrfspec)) {
      cat("    • HRF Type:", crayon::yellow(attr(term$hrfspec$hrf, "name")), "\n")
      if (!is.null(attr(term$hrfspec$hrf, "nbasis"))) {
        cat("    • Basis Functions:", crayon::yellow(attr(term$hrfspec$hrf, "nbasis")), "\n")
      }
    }
    
    if (!is.null(term$contrasts) && length(term$contrasts) > 0) {
      cat("    • Contrasts:", crayon::cyan(paste(names(term$contrasts), collapse = ", ")), "\n")
    }
  }
  
  cat("\n")
}

#' Plot an event_model object
#'
#' @description
#' Creates a detailed plot of an event_model object by visualizing its design matrix.
#' The design matrix is first converted into a long-format tibble and then plotted over time,
#' with separate panels for each block. Several customization options are available.
#'
#' @param x An event_model object.
#' @param y Unused.
#' @param term_name Optional: Name of the term to plot. If NULL, the first term is plotted.
#' @param longnames Logical flag; if TRUE, use long condition names; if FALSE, use short names. Default is TRUE.
#' @param title Optional title for the plot. If not provided, a default title is generated.
#' @param xlab Label for the x-axis. Default is "Time".
#' @param ylab Label for the y-axis. Default is "Value".
#' @param line_size Numeric value for the line thickness in the plot. Default is 1.
#' @param color_palette The name of a ColorBrewer palette to use. Default is "Set1".
#' @param facet_ncol Number of columns for facet_wrap. Default is 1.
#' @param theme_custom A ggplot2 theme to apply. Default is theme_bw with base_size 14.
#' @param hide_legend_threshold If the number of unique conditions exceeds this threshold, the legend is hidden. Default is 25.
#' @param ... Additional arguments passed to ggplot2::geom_line().
#'
#' @importFrom ggplot2 ggplot aes_string geom_line facet_wrap labs theme_bw scale_color_brewer guides
#' @importFrom tidyr gather
#' @return A ggplot2 plot object.
#' @export
plot.event_model <- function(x, y, term_name = NULL, longnames = TRUE,
                             title = NULL, xlab = "Time", ylab = "Value",
                             line_size = 1, color_palette = "Set1", facet_ncol = 1,
                             theme_custom = ggplot2::theme_bw(base_size = 14),
                             hide_legend_threshold = 25, ...) {
  
  # Extract event terms and their names.
  all_terms <- terms(x)
  term_names <- sapply(all_terms, "[[", "varname")
  
  # Retrieve the sampling frame and ensure it has the required fields.
  sframe <- x$sampling_frame
  if (!("time" %in% names(sframe)) || !("blockids" %in% names(sframe))) {
    stop("The sampling_frame must contain 'time' and 'blockids' components.")
  }
  
  # For each event term, convert its design matrix to long format.
  dflist <- lapply(all_terms, function(term) {
    dm <- suppressMessages(tibble::as_tibble(design_matrix(term), .name_repair = "check_unique"))
    if (!longnames) {
      names(dm) <- shortnames(term)
    }
    # Add block and time information from the sampling frame.
    dm$.block <- sframe$blockids
    dm$.time <- sframe$time
    # Convert from wide to long format.
    tidyr::gather(dm, key = "condition", value = "value", - .time, - .block)
  })
  names(dflist) <- term_names
  
  # Select the term to plot.
  if (is.null(term_name)) {
    plot_term <- term_names[1]
    dfx <- dflist[[plot_term]]
  } else {
    if (!(term_name %in% term_names)) {
      stop("The specified term_name was not found among event terms.")
    }
    plot_term <- term_name
    dfx <- dflist[[term_name]]
  }
  
  # Build the ggplot.
  p <- ggplot2::ggplot(dfx, ggplot2::aes_string(x = ".time", y = "value", colour = "condition")) +
    ggplot2::geom_line(size = line_size, ...) +
    ggplot2::facet_wrap(~ .block, ncol = facet_ncol) +
    ggplot2::labs(title = if (!is.null(title)) title else paste("Event Model:", plot_term),
                  x = xlab, y = ylab, colour = "Condition") +
    ggplot2::scale_color_brewer(palette = color_palette) +
    theme_custom
  
  # Hide the legend if there are too many unique conditions.
  if (length(unique(dfx$condition)) > hide_legend_threshold) {
    p <- p + ggplot2::guides(colour = "none")
  }
  
  p
}

#' Detect variables inside exprs(...) calls and copy those columns into `events`.
#' 
#' @param event_terms The list of event term specifications
#' @param events A data.frame of event variables
#' @return The updated `events` data frame, ensuring all variables found in expressions exist as columns
#' @keywords internal
.extract_expr_variables_into_events <- function(event_terms, events) {
  
  for (et in event_terms) {
    vars <- et$variables
    
    # If expression, parse out the symbols inside.
    for (v in vars) {
      if (rlang::is_quosure(v)) {
        # Convert to call first
        vcall <- rlang::get_expr(v)
        needed <- .all_symbols_in_call(vcall)
      } else if (rlang::is_call(v)) {
        needed <- .all_symbols_in_call(v)
      } else {
        # If it's just a symbol or character, that is simpler:
        needed <- character(0)
        if (rlang::is_symbol(v)) {
          needed <- as.character(v)
        } else if (is.character(v)) {
          needed <- v
        }
      }
      # For each variable name in 'needed', ensure events has that column
      # If missing, throw an error or do something minimal:
      for (nm in needed) {
        if (! nm %in% names(events)) {
          stop(sprintf(
            "create_event_model: The expression references variable '%s' which is not a column in `events`.", 
            nm
          ))
        }
      }
    }
  }
  
  events
}

#' Recursively find all symbol names used inside a call/expression.
#' 
#' For example, given `Poly(modulator, 2)`, returns `c("modulator")`.
#' 
#' @param call An R call object or symbol
#' @return A character vector of bare variable names
#' @keywords internal
.all_symbols_in_call <- function(call) {
  out <- character(0)
  if (rlang::is_symbol(call)) {
    return(as.character(call))
  }
  if (!rlang::is_call(call)) {
    return(character(0))
  }
  for (arg in as.list(call)[-1]) {
    out <- c(out, .all_symbols_in_call(arg))
  }
  unique(out)
}

#' Create an event model directly from components.
#'
#' Constructs an event model from its components. Variables can be specified as 
#' character strings or expressions (captured using rlang::exprs). The function creates 
#' HRF specifications (hrfspec objects), builds a formula string, and then assembles 
#' the final event model.
#'
#' @param event_terms A list of event term specifications. Each term is a list with components:
#'        - `variables`: Character vector of variable names or expressions.
#'        - `hrf`: HRF specification for the term.
#' @param events A data frame of event variables (all variables must be present).
#' @param onsets A numeric vector of event onset times (in seconds).
#' @param block A vector of block IDs (must be strictly increasing integers).
#' @param sampling_frame The time series grid over which to sample the function.
#' @param durations A numeric vector of event durations (default is 0 for all events).
#' @param drop_empty Logical indicating whether to drop empty events (default is TRUE).
#' @param precision Numeric value indicating the precision of HRF sampling (default is 0.3).
#'
#' @return An event_model object.
#'
#' @examples
#' library(fmrireg)
#' library(rlang)
#'
#' # Example with variables as character strings:
#' event_terms <- list(
#'   list(
#'     variables = c("x", "y"),
#'     hrf = "spmg1"
#'   )
#' )
#'
#' # Example with variables as expressions:
#' event_terms <- list(
#'   list(
#'     variables = exprs(x, Poly(y, 2)),
#'     hrf = "spmg1"
#'   )
#' )
#'
#' @export
create_event_model <- function(event_terms,
                               events,
                               onsets,
                               block,
                               sampling_frame,
                               durations = 0,
                               drop_empty = TRUE,
                               precision = 0.3) {
  
  # Input validation.
  assert_that(is.list(event_terms), msg = "event_terms must be a list of term specifications")
  assert_that(length(onsets) == nrow(events), msg = "onsets length must match number of events")
  
  # Process block IDs.
  if (is.factor(block)) {
    block <- as.integer(as.character(block))
  }
  assert_that(is.increasing(block), msg = "'block' must consist of strictly increasing integers")
  assert_that(length(block) == nrow(events))
  
  blocks <- unique(block)
  ranked_blocks <- rank(blocks)
  blockids <- ranked_blocks[match(block, blocks)]
  
  # Process durations.
  if (length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  
  #browser()
  # MINIMAL FIX: ensure that *all* variables in each hrfspec 
  # are actually in `events` data frame. If an hrfspec has
  # `variables = exprs(Poly(modulator, 2))`, we look at 
  # the expression to find "modulator".
  events <- .extract_expr_variables_into_events(event_terms, events)
  
  # Create a data environment from the event data.
  data_env <- list2env(as.list(events), parent = as.environment("package:fmrireg"))
  
  # Create hrfspec objects from event_terms.
  hrfspecs <- lapply(event_terms, function(term) {
    variables <- term$variables
    hrf_spec <- term$hrf
    
    # Prepare variables: if characters, convert to symbols, etc.
    if (is.character(variables)) {
      variables <- lapply(variables, rlang::sym)
    } else if (is.expression(variables) || rlang::is_quosure(variables)) {
      variables <- list(variables)
    } else if (!is.list(variables)) {
      stop("variables in event_terms must be a character vector or expressions captured by rlang::exprs")
    }
    
    # Prepare HRF basis.
    if (is.character(hrf_spec)) {
      hrf <- getHRF(hrf_spec)
    } else if (inherits(hrf_spec, "HRF")) {
      hrf <- hrf_spec
    } else if (is.list(hrf_spec)) {
      hrf <- do.call(getHRF, c(list(hrf_spec$hrf), hrf_spec$parameters))
    } else {
      stop("Invalid HRF specification")
    }
    
    hrfspec(
      vars = variables,
      basis = hrf,
      onsets = onsets,
      durations = durations,
      precision = precision,
      data_env = data_env
    )
  })
  
  # Construct a formula string from event_terms.
  formula_terms <- lapply(event_terms, function(term) {
    variables <- term$variables
    hrf_spec <- term$hrf
    var_strings <- sapply(variables, function(var) {
      if (is.symbol(var) || is.language(var)) deparse(var)
      else if (is.character(var)) var
      else stop("variables must be symbols, expressions, or character strings")
    })
    hrf_type <- if (is.character(hrf_spec)) {
      hrf_spec
    } else if (inherits(hrf_spec, "HRF")) {
      attr(hrf_spec, "name")
    } else if (is.list(hrf_spec)) {
      hrf_spec$hrf
    }
    var_string <- paste(var_strings, collapse = ", ")
    paste0("hrf(", var_string, ", basis='", hrf_type, "')")
  })
  
  formula_str <- paste("onsets ~", paste(formula_terms, collapse = " + "))
  model_formula <- as.formula(formula_str)
  
  # model_spec <- list(
  #   formula = model_formula,
  #   event_table = events,
  #   onsets = onsets,
  #   event_spec = hrfspecs,
  #   blockvals = block,
  #   blockids = blockids,
  #   durations = durations,
  #   sampling_frame = sampling_frame,
  #   drop_empty = drop_empty,
  #   precision = precision
  # )
  # 
  # class(model_spec) <- c("event_model_spec", "list")
  # 
  # construct_model(model_spec)
  
  form <- formula(formula_str)
  event_model(form, events, block, sampling_frame, drop_empty=drop_empty, durations=durations, precision=precision)
}

plotly.event_model <- function(x, 
                               term_name = NULL,
                               title = NULL,
                               xlab = "Time",
                               ylab = "Value",
                               line_size = 2,
                               hide_legend_threshold = 25,
                               ...) {
  library(dplyr)
  library(tidyr)
  library(plotly)
  
  all_terms  <- terms(x)
  term_names <- sapply(all_terms, `[[`, "varname")
  sframe     <- x$sampling_frame
  stopifnot("time" %in% names(sframe), "blockids" %in% names(sframe))
  
  # Helper: convert one event_term -> long tibble with .Term label
  to_long <- function(term, lbl) {
    dm <- design_matrix(term)
    df <- as_tibble(dm, .name_repair="unique") %>%
      mutate(.time  = sframe$time,
             .block = sframe$blockids) %>%
      pivot_longer(cols = -c(.time, .block), 
                   names_to = "condition", 
                   values_to = "value") %>%
      mutate(.Term = lbl)
    df
  }
  
  # If a specific term is requested, just plot it (no animation needed).
  if (!is.null(term_name)) {
    if (!(term_name %in% term_names))
      stop("term_name not found among: ", paste(term_names, collapse=", "))
    idx <- match(term_name, term_names)
    df  <- to_long(all_terms[[idx]], term_name)
    
    p <- plot_ly(df, 
                 x = ~.time, y = ~value,
                 split = ~condition,       # each condition -> separate trace
                 type="scatter", mode="lines",
                 line = list(width=line_size),
                 ...) %>%
      layout(title = title %||% paste("Event Model:", term_name),
             xaxis = list(title=xlab),
             yaxis = list(title=ylab))
    if (length(unique(df$condition)) > hide_legend_threshold) {
      p <- p %>% layout(showlegend = FALSE)
    }
    return(p)
  }
  
  # Otherwise gather all terms so that we can animate from frame to frame
  big_df <- purrr::map2_dfr(all_terms, term_names, to_long)
  
  # Convert .Term and condition to factor
  # (We want them consistent across frames.)
  all_terms_char <- unique(big_df$.Term)
  all_conds_char <- unique(big_df$condition)
  
  big_df <- big_df %>%
    mutate(.Term     = factor(.Term,     levels = all_terms_char),
           condition = factor(condition, levels = all_conds_char))
  
  # "Complete" so each frame has every .time × condition combo:
  # fill missing with NA
  big_df <- big_df %>%
    complete(.Term, .time, condition, fill = list(value = NA, .block=NA))
  
  # Build the plotly with frame=~.Term for animation
  p <- plot_ly(
    data = big_df,
    x = ~.time,
    y = ~value,
    split = ~condition,
    frame = ~.Term,
    type  = "scatter",
    mode  = "lines",
    line  = list(width=line_size),
    ...
  ) %>%
    layout(
      title = title %||% "Event Model (All Terms)",
      xaxis = list(title=xlab),
      yaxis = list(title=ylab)
    ) %>%
    # Force redraw so old lines don't persist
    animation_opts(frame=1, transition=0, redraw=TRUE)
  
  if (length(all_conds_char) > hide_legend_threshold) {
    p <- p %>% layout(showlegend = FALSE)
  }
  p
}


#' correlation_map.event_model
#'
#' @description
#' Generates a correlation heatmap of the columns in an \code{event_model}'s design matrix.
#'
#' @param x An \code{event_model}.
#' @param method Correlation method, "pearson" or "spearman".
#' @param half_matrix Logical; if TRUE, mask out the upper-triangle.
#' @param absolute_limits Logical; if TRUE, fix color scale to -1..+1.
#' @param ... Further arguments passed to \code{\link{.correlation_map_common}}.
#' @export
correlation_map.event_model <- function(x,
                                        method          = c("pearson", "spearman"),
                                        half_matrix     = FALSE,
                                        absolute_limits = TRUE,
                                        ...) {
  DM <- as.matrix(design_matrix(x))
  .correlation_map_common(DM, method=method, half_matrix=half_matrix,
                          absolute_limits=absolute_limits, ...)
}


#' Visualize the entire design matrix of an event_model as a heatmap
#'
#' @description
#' Produces a single heatmap of all columns in the design matrix for an
#' `event_model` object, with rows corresponding to scan number and columns
#' corresponding to regressors. A white horizontal line separates each block/run.
#'
#' @param x An `event_model` object.
#' @param block_separators Logical; if `TRUE`, draw white horizontal lines between blocks.
#' @param rotate_x_text Logical; if `TRUE`, rotate x-axis labels by 45 degrees.
#' @param fill_midpoint Numeric or `NULL`; if not `NULL`, passed as the `midpoint`
#'   argument to [ggplot2::scale_fill_gradient2()] to center the color scale. If
#'   your design matrix has both positive and negative values, you may want
#'   `fill_midpoint = 0`.
#' @param fill_limits Numeric vector of length 2 or `NULL`; passed to `limits=`
#'   in the fill scale. Useful to clip color range.
#' @param ... Other arguments passed on to `geom_tile()`.
#'
#' @import ggplot2
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#'
#' @return A ggplot2 plot object.
#' @export
#' Visualize the entire design matrix of an event_model as a heatmap, with optional rescaling
#'
#' @param x An `event_model` object.
#' @param rescale_cols Logical; if TRUE, columns are rescaled to (-1,1) before plotting.
#' @inheritParams design_matrix_heatmap.event_model
#' @export
design_map.event_model <- function(x,
                                         rescale_cols    = TRUE,
                                         block_separators = TRUE,
                                         rotate_x_text    = TRUE,
                                         fill_midpoint    = NULL,
                                         fill_limits      = NULL,
                                         ...) {
  # 1) Extract the design matrix
  DM <- design_matrix(x)
  
  # 2) Optionally rescale each regressor column to (-1,1)
  if (rescale_cols) {
    DM <- rescale_columns_01_to_neg1_pos1(DM)
  }
  
  # 3) Convert to long format
  n_scans <- nrow(DM)
  df_long <- tibble::as_tibble(DM, .name_repair = "unique")
  df_long$scan_number <- seq_len(n_scans)  
  df_long <- tidyr::pivot_longer(
    df_long,
    cols      = -scan_number,
    names_to  = "Regressor",
    values_to = "Value"
  )
  
  # 4) Build ggplot
  plt <- ggplot(df_long, aes(x = Regressor, y = scan_number, fill = Value)) +
    geom_tile(...)
  
  # Reverse y-axis so that scan_number=1 is at top
  plt <- plt + scale_y_reverse()
  
  # Color scale: either normal gradient or centered
  if (is.null(fill_midpoint)) {
    plt <- plt + scale_fill_gradientn(colours = c("navy","white","firebrick"),
                                      limits  = fill_limits)
  } else {
    plt <- plt + scale_fill_gradient2(midpoint = fill_midpoint,
                                      low      = "navy",
                                      mid      = "white",
                                      high     = "firebrick",
                                      limits   = fill_limits)
  }
  
  # Draw white lines for block boundaries
  if (block_separators && !is.null(x$blockids)) {
    block_ids    <- x$blockids
    runs         <- rle(block_ids)
    row_breaks   <- cumsum(runs$lengths)
    n_regressors <- ncol(DM)
    for (rb in row_breaks[-length(row_breaks)]) {
      plt <- plt + 
        annotate("segment",
                 x    = 0.5, 
                 xend = n_regressors + 0.5,
                 y    = rb + 0.5,
                 yend = rb + 0.5,
                 color = "white", size = 1)
    }
  }
  
  # Rotation of x-axis text
  plt <- plt + 
    theme_minimal(base_size = 14) +
    labs(x = "Regressors", y = "Scan Number", fill = "Value") +
    theme(
      panel.grid  = element_blank(),
      axis.text.x = if (rotate_x_text) element_text(angle = 45, hjust = 1) else element_text()
    )
  
  plt
}


#' Rescale columns of a design matrix to (-1,1).
#'
#' @param DM A numeric matrix or data frame (e.g. from `design_matrix(event_model)`).
#' @param handle_constant Whether to force constant (zero-range) columns to zero, or
#'        raise an error, or skip them. Default "to_zero" sets them all to 0.
#' @return A matrix the same size as `DM` with each column scaled to (-1,1).
#' @export
rescale_columns_01_to_neg1_pos1 <- function(DM, handle_constant = c("to_zero","error","skip")) {
  handle_constant <- match.arg(handle_constant)
  
  DM_scaled <- as.data.frame(DM)  # safe copy
  for (j in seq_len(ncol(DM_scaled))) {
    colvals <- DM_scaled[[j]]
    rng     <- range(colvals, na.rm = TRUE)
    minv    <- rng[1]
    maxv    <- rng[2]
    span    <- maxv - minv
    
    if (abs(span) < 1e-15) {
      # It's effectively constant
      if (handle_constant == "to_zero") {
        DM_scaled[[j]] <- 0
      } else if (handle_constant == "error") {
        stop(sprintf("Column %d is constant; cannot scale to (-1,1).", j))
      } else if (handle_constant == "skip") {
        # leave it as is
      }
    } else {
      # map [minv, maxv] -> [-1,1]
      DM_scaled[[j]] <- 2 * (colvals - minv) / span - 1
    }
  }
  as.matrix(DM_scaled)
}
###############################################################################
## End of event_model.R
###############################################################################