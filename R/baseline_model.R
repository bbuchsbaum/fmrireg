###############################################################################
## BASELINE_MODEL.R
##
## This file implements baseline model construction for fMRI analyses.
## It includes helper routines for:
##   - Extracting column indices from a list of matrices.
##   - Constructing baseline models that account for drift, block‐wise
##     intercepts, and nuisance regressors.
##   - Creating design matrices and terms for baseline models.
##   - Specifying block and nuisance terms.
##   - Printing and plotting baseline models.
##
## The file has been organized into sections with additional inline comments
## for clarity. Functionality is preserved exactly as in the original file.
###############################################################################

## ============================================================================
## Section 1: Helper Functions
## ============================================================================

#' Get Column Indices from a List of Matrices
#'
#' Given a list of matrices, returns a list in which each element is a sequence
#' of column indices corresponding to each matrix.
#'
#' @param Xlist A list of matrices.
#' @return A list of integer sequences representing column indices.
#' @importFrom purrr map_int
#' @noRd
#' @keywords internal
get_col_inds <- function(Xlist) {
  ncols <- purrr::map_int(Xlist, ncol)
  csum <- cumsum(ncols)
  csum1 <- c(0, csum[-length(csum)])
  m <- as.matrix(cbind(csum1 + 1, csum))
  
  # Return a list of sequences from the start to the end index for each matrix.
  lapply(1:nrow(m), function(i) seq(m[i, 1], m[i, 2]))
  # Alternative: apply(m, 1, function(x) seq(x[1], x[2]))
}


## ============================================================================
## Section 2: Baseline Model Construction and Specification
## ============================================================================

#' Construct a Baseline Model.
#'
#' Builds a baseline model to account for noise and non–event-related variance.
#' This model may include a drift term, a block intercept term, and nuisance regressors.
#'
#' @param basis Character; type of basis function ("constant", "poly", "bs", or "ns").
#' @param degree Integer; degree of the spline/polynomial function.
#' @param sframe A sampling_frame object.
#' @param intercept Character; whether to include an intercept ("runwise", "global", or "none").
#'   (Automatically set to FALSE when basis == "constant".)
#' @param nuisance_list Optional list of nuisance matrices (one matrix per fMRI block).
#'
#' @return An object of class "baseline_model".
#'
#' @examples 
#' sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)
#' bmod <- baseline_model(basis = "bs", degree = 3, sframe = sframe)
#' bmod_global <- baseline_model(basis = "bs", degree = 3, sframe = sframe, intercept = "global")
#' bmod_nointercept <- baseline_model(basis = "bs", degree = 3, sframe = sframe, intercept = "none")
#' stopifnot(ncol(design_matrix(bmod)) == 8)
#' @export
baseline_model <- function(basis = c("constant", "poly", "bs", "ns"), degree = 1, sframe, 
                           intercept = c("runwise", "global", "none"), nuisance_list = NULL) {
  
  basis <- match.arg(basis)
  intercept <- match.arg(intercept)
  
  if (basis %in% c("bs", "ns")) {
    assert_that(degree > 2, msg ="'bs' and 'ns' bases must have degree >= 3")
  }
  
  # Construct the drift term using a baseline specification.
  drift_spec <- baseline(degree = degree, basis = basis, intercept = intercept)
  drift_term <- construct(drift_spec, sframe)
  
  # Construct a block intercept term if applicable.
  block_term <- if (intercept != "none" && basis != "constant") {
    construct_block_term("constant", sframe, intercept)
  }
  
  # Construct a nuisance term if a nuisance_list is provided.
  nuisance_term <- if (!is.null(nuisance_list)) {
    total_len <- sum(purrr::map_int(nuisance_list, nrow))
    assertthat::assert_that(total_len == length(blockids(sframe)), 
                            msg = "number of rows of `nuisance_list` must match number of scans in sampling_frame")
    assertthat::assert_that(length(nuisance_list) == length(blocklens(sframe)), 
                            msg = "number of `nuisance_list` elements must match number of blocks in sampling_frame")
    
    colind <- get_col_inds(nuisance_list)
    rowind <- split(1:length(blockids(sframe)), blockids(sframe))
    
    for (i in seq_along(nuisance_list)) {
      nmat <- nuisance_list[[i]]
      colnames(nmat) <- paste0("nuisance#", i, "#", 1:ncol(nmat))
      # Use unclass() because Matrix::bdiag may not work with additional classes.
      nuisance_list[[i]] <- unclass(as.matrix(nmat))
    }
    
    cnames <- unlist(lapply(nuisance_list, colnames))
    nmat <- as.matrix(Matrix::bdiag(nuisance_list))
    colnames(nmat) <- cnames
    baseline_term("nuisance", nmat, colind, rowind)
  }
  
  ret <- list(drift_term = drift_term, 
              drift_spec = drift_spec, 
              block_term = block_term, 
              nuisance_term = nuisance_term, 
              sampling_frame = sframe)
  
  class(ret) <- c("baseline_model", "list")
  ret
}

#' Create a Baseline Specification.
#'
#' Generates a baselinespec for modeling low-frequency drift in fMRI time series.
#'
#' @param degree Number of basis terms per image block (ignored for "constant").
#' @param basis Type of basis ("constant", "poly", "bs", or "ns").
#' @param name Optional name for the term.
#' @param intercept Type of intercept to include ("runwise", "global", or "none").
#'
#' @return A baselinespec list instance.
#' @export
baseline <- function(degree = 1, basis = c("constant", "poly", "bs", "ns"), name = NULL, 
                     intercept = c("runwise", "global", "none")) {
  
  basis <- match.arg(basis)
  intercept <- match.arg(intercept)
  
  if (basis == "constant") {
    degree <- 1
  }
  
  bfun <- switch(basis,
                 bs = splines::bs,
                 ns = splines::ns,
                 poly = poly,
                 constant = function(x, degree) { matrix(rep(1, length(x))) })
  
  if (is.null(name)) {
    name <- paste0("baseline_", basis, "_", degree)
  }
  
  ret <- list(
    degree = degree,
    basis = basis,
    fun = bfun,
    intercept = intercept,
    name = name
  )
  
  class(ret) <- c("baselinespec", "nuisancespec")
  ret
}


## ============================================================================
## Section 3: Design Matrix and Term Functions for Baseline Models
## ============================================================================

#' Construct a Design Matrix for a Baseline Model.
#'
#' Combines the drift term, block intercept term, and nuisance term (if any)
#' into a complete design matrix.
#'
#' @param x A baseline_model object.
#' @param blockid Optional block ID to extract a subset.
#' @param allrows Logical; if TRUE, returns all rows for the block.
#' @return A tibble representing the design matrix.
#' @export
design_matrix.baseline_model <- function(x, blockid = NULL, allrows = FALSE, ...) {
  if (is.null(x$nuisance_term) && is.null(x$block_term)) {
    suppressMessages(tibble::as_tibble(design_matrix(x$drift_term, blockid, allrows), .name_repair = "check_unique"))
  } else if (is.null(x$nuisance_term) && !is.null(x$block_term)) {
    suppressMessages(tibble::as_tibble(cbind(design_matrix(x$block_term, blockid, allrows), 
                                             design_matrix(x$drift_term, blockid, allrows)), .name_repair = "check_unique"))
  } else if (!is.null(x$nuisance_term) && is.null(x$block_term)) {
    suppressMessages(tibble::as_tibble(cbind(design_matrix(x$drift_term, blockid, allrows),
                                             design_matrix(x$nuisance_term, blockid, allrows)), .name_repair = "check_unique"))
  } else {
    suppressMessages(tibble::as_tibble(cbind(design_matrix(x$block_term, blockid, allrows), 
                                             design_matrix(x$drift_term, blockid, allrows), 
                                             design_matrix(x$nuisance_term, blockid, allrows)), .name_repair = "check_unique"))
  }
}

#' Return the Terms of a Baseline Model.
#'
#' @param x A baseline_model object.
#' @return A list of terms (drift, block, nuisance).
#' @export
#' @importFrom purrr map_lgl
terms.baseline_model <- function(x, ...) {
  ret <- list(x$block_term, x$drift_term, x$nuisance_term)
  ret <- ret[!purrr::map_lgl(ret, is.null)]
  names(ret) <- unlist(lapply(ret, "[[", "varname"))
  ret
}

#' Retrieve Cells of a Baseline Model.
#'
#' Combines the cells from all baseline terms into a tibble with an index.
#'
#' @param x A baseline_model object.
#' @return A tibble with columns: term, level, basis, and index.
#' @export
#' @importFrom dplyr mutate
cells.baseline_model <- function(x, ...) {
  .terms <- terms(x)
  do.call(rbind, lapply(.terms, function(t) {
    ncond <- length(conditions(t))
    zstr <- paste0(rep("0", ceiling(log10(ncond + 1e-6))), collapse = "")
    dplyr::tibble(term = t$varname, level = conditions(t), basis = paste0("basis", zstr, 1:length(conditions(t))))
  })) %>% dplyr::mutate(index = 1:dplyr::n())
}

## ============================================================================
## Section 4: Block and Nuisance Specification Helpers
## ============================================================================

#' Create a Block Variable.
#'
#' Returns a block variable that is constant over the span of a scanning run.
#'
#' @param x The block variable.
#' @return An object of class "blockspec".
#' @export
block <- function(x) {
  varname <- substitute(x)
  pterm <- parse_term(as.list(substitute(x)), "block")
  ret <- list(
    name = varname,
    label = pterm$label
  )
  class(ret) <- "blockspec"
  ret
}

#' Construct a Baseline Specification from a Sampling Frame.
#'
#' Given a baselinespec object and a sampling frame, constructs the baseline
#' covariates for each block by applying the baseline function to each block.
#'
#' @param x A baselinespec object.
#' @param model_spec A model specification containing (or serving as) the sampling frame.
#' @return A baseline term object.
#' @noRd
construct.baselinespec <- function(x, model_spec, ...) {
  sampling_frame <- if (!is.null(model_spec$sampling_frame)) model_spec$sampling_frame else model_spec
  
  # Compute baseline covariates for each block.
  ret <- lapply(sampling_frame$blocklens, function(bl) x$fun(seq(1, bl), degree = x$degree))
  
  if (x$basis == "constant" && x$intercept == "global") {
    colind <- lapply(ret, function(x) 1:1)
    rowind <- lapply(seq_along(ret), function(i) {
      rowstart <- sum(sampling_frame$blocklens[1:i]) - sampling_frame$blocklens[i] + 1
      rowend <- rowstart + sampling_frame$blocklens[i] - 1
      rowstart:rowend
    })
    ret <- do.call(rbind, ret)
    cnames <- paste0("base_", x$basis)
    colnames(ret) <- cnames
    baseline_term(x$name, ret, colind, rowind)
  } else {
    nc <- sum(purrr::map_int(ret, ncol))
    nc_per_block <- ncol(ret[[1]])
    mat <- matrix(0, sum(sampling_frame$blocklens), nc)
    
    colind <- vector("list", length(ret))
    rowind <- vector("list", length(ret))
    
    for (i in seq_along(ret)) {
      rowstart <- sum(sampling_frame$blocklens[1:i]) - sampling_frame$blocklens[i] + 1
      rowend <- rowstart + sampling_frame$blocklens[i] - 1
      colstart <- (nc_per_block) * (i - 1) + 1
      colend <- colstart + nc_per_block - 1
      mat[rowstart:rowend, colstart:colend] <- ret[[i]]
      colind[[i]] <- colstart:colend
      rowind[[i]] <- rowstart:rowend
    }
    
    cnames <- apply(expand.grid(paste0("base_", x$basis, 1:nc_per_block), 
                                paste0("block_", 1:length(sampling_frame$blocklens)), 
                                stringsAsFactors = TRUE), 1, paste, collapse = "_")
    colnames(mat) <- cnames
    baseline_term(x$name, mat, colind, rowind)	
  }
}

#' Construct a Baseline Term.
#'
#' Creates a baseline_term object given a covariate matrix and its associated
#' column and row indices.
#'
#' @param varname The name of the term.
#' @param mat A matrix (or data frame) of covariates.
#' @param colind A list of column indices.
#' @param rowind A list of row indices.
#' @return A baseline_term object.
#' @importFrom tibble as_tibble
#' @noRd
#' @keywords internal
baseline_term <- function(varname, mat, colind, rowind) {
  stopifnot(inherits(mat, "matrix") || is.data.frame(mat) || inherits(mat, "Matrix"))
  ret <- list(varname = varname, 
              design_matrix = suppressMessages(tibble::as_tibble(as.matrix(mat), .name_repair = "minimal")), 
              colind = colind, 
              rowind = rowind)
  class(ret) <- c("baseline_term", "matrix_term", "fmri_term", "list")
  ret
}

#' @export
design_matrix.baseline_term <- function(x, blockid = NULL, allrows = FALSE, ...) {
  if (is.null(blockid)) {
    x$design_matrix
  } else {
    if (!allrows) {
      x$design_matrix[unlist(x$rowind[blockid]), unlist(x$colind[blockid]), drop = FALSE]
    } else {
      x$design_matrix[, unlist(x$colind[blockid]), drop = FALSE]
    }
  }
}

#' Construct a Block Term.
#'
#' Constructs a constant block intercept term based on block IDs.
#'
#' @param vname The name of the block variable.
#' @param sframe A sampling_frame object.
#' @param intercept Type of intercept ("global" or "runwise").
#' @return A block_term object.
#' @noRd
construct_block_term <- function(vname, sframe, intercept = c("global", "runwise")) {
  intercept <- match.arg(intercept)
  blockids <- blockids(sframe)
  blockord <- sort(unique(blockids))
  
  # Create an expanded ordered factor for block IDs.
  expanded_blockids <- ordered(rep(blockord, blocklens(sframe)))
  
  if (length(blockord) == 1 || intercept == "global") {
    mat <- matrix(rep(1, length(expanded_blockids)), nrow = length(expanded_blockids), ncol = 1)
    colnames(mat) <- paste0(vname, "_global")
  } else {
    mat <- model.matrix(~ expanded_blockids - 1)
    colnames(mat) <- paste0(vname, "_", blockord)
  }
  
  block_term(vname, blockids = blockids, expanded_blockids = expanded_blockids, mat, type = intercept)
}

#' Construct a Block Term Object.
#'
#' Constructs a constant term for block intercepts.
#'
#' @param varname The name of the term.
#' @param blockids The ordered block IDs.
#' @param expanded_blockids The expanded block IDs (by run).
#' @param mat The matrix of covariates.
#' @param type The type ("runwise" or "global").
#' @importFrom tibble as_tibble
#' @noRd
#' @keywords internal
block_term <- function(varname, blockids, expanded_blockids, mat, type = c("runwise", "global")) {
  type <- match.arg(type)
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(nrow(mat) == length(blockids))
  
  if (type == "runwise") {
    assertthat::assert_that(ncol(mat) == length(unique(blockids)))
  } else {
    assertthat::assert_that(ncol(mat) == 1)
  }
  
  ret <- list(varname = varname, 
              blockids = blockids, 
              expanded_blockids = expanded_blockids, 
              design_matrix = suppressMessages(tibble::as_tibble(mat, .name_repair = "check_unique")), 
              nblocks = ncol(mat),
              colind = as.list(1:ncol(mat)), 
              rowind = split(1:length(blockids), blockids))
  class(ret) <- c("block_term", "matrix_term", "fmri_term", "list")
  ret
}

#' @export
design_matrix.block_term <- function(x, blockid = NULL, allrows = FALSE, ...) {
  if (is.null(blockid)) {
    x$design_matrix
  } else {
    if (!allrows) {
      x$design_matrix[unlist(x$rowind[blockid]), unlist(x$colind[blockid]), drop = FALSE]
    } else {
      x$design_matrix[, unlist(x$colind[blockid]), drop = FALSE]
    }
  }
}

#' @export
term_names.baseline_model <- function(x) {
  xt <- terms(x)
  unlist(lapply(xt, function(term) term$varname))
}

#' Create a Nuisance Specification.
#'
#' Returns a nuisance term specification from a numeric matrix.
#'
#' @param x A matrix.
#' @return An object of class "nuisancespec".
#' @export
nuisance <- function(x) {
  varname <- substitute(x)
  ret <- list(name = varname)
  class(ret) <- "nuisancespec"
  ret
}

#' @export
construct.nuisancespec <- function(x, model_spec, ...) {
  mat <- eval(parse(text = x$varname), envir = model_spec$aux_data, enclos = parent.frame())
  matrix_term(x$name, mat)
}

#' @export
construct.blockspec <- function(x, model_spec, ...) {
  construct_block_term(x$name, model_spec$sampling_frame)
}


## ============================================================================
## Section 4: Print and Plot Methods for Baseline Models
## ============================================================================

#' Print a Baseline Model.
#'
#' Displays key information about the baseline model components and a preview
#' of the design matrix.
#'
#' @param x A baseline_model object.
#' @param ... Additional arguments (ignored).
#' @export
print.baseline_model <- function(x, ...) {
  # Extract component information.
  drift_name <- x$drift_term$varname
  basis_type <- x$drift_spec$basis
  degree <- x$drift_spec$degree
  drift_cols <- ncol(design_matrix(x$drift_term))
  const_cols <- if (!is.null(x$block_term)) ncol(design_matrix(x$block_term)) else 0
  nuis_cols <- if (!is.null(x$nuisance_term)) ncol(design_matrix(x$nuisance_term)) else 0
  total_cols <- ncol(design_matrix(x))
  
  # Print header.
  cat("╔══════════════════════════════════════════╗\n")
  cat("║           Baseline Model                 ║\n")
  cat("╠══════════════════════════════════════════╣\n")
  
  # Drift term info.
  cat("║ Drift Components                         ║\n")
  cat(sprintf("║   • %-34s ║\n", paste("Name:", drift_name)))
  cat(sprintf("║   • %-34s ║\n", paste("Basis type:", basis_type)))
  cat(sprintf("║   • %-34s ║\n", paste("Degree:", degree)))
  cat(sprintf("║   • %-34s ║\n", paste("Drift columns:", drift_cols)))
  
  # Additional components.
  cat("║                                          ║\n")
  cat("║ Additional Components                    ║\n")
  cat(sprintf("║   • %-34s ║\n", paste("Constant columns:", const_cols)))
  cat(sprintf("║   • %-34s ║\n", paste("Nuisance columns:", nuis_cols)))
  
  # Summary.
  cat("║                                          ║\n")
  cat("║ Model Summary                            ║\n")
  cat(sprintf("║   • %-34s ║\n", paste("Total columns:", total_cols)))
  
  # Preview design matrix.
  cat("║                                          ║\n")
  cat("║ Design Matrix Preview                    ║\n")
  dm <- head(design_matrix(x), 3)
  for (i in 1:min(3, nrow(dm))) {
    row_preview <- paste(sprintf("%6.3f", as.numeric(dm[i, 1:min(4, ncol(dm))])), collapse = " ")
    if (ncol(dm) > 4) row_preview <- paste0(row_preview, " ...")
    cat(sprintf("║   %-36s ║\n", row_preview))
  }
  if (nrow(dm) > 3) cat("║   ...                                   ║\n")
  
  cat("╚══════════════════════════════════════════╝\n")
}

#' Plot a Baseline Model.
#'
#' Creates a detailed ggplot2 visualization of the baseline model design matrix.
#' Each non-constant term is plotted over time. The plot includes separate panels
#' for each block and supports customization of titles, axis labels, line size, and color palette.
#'
#' @param x A baseline_model object.
#' @param term_name Optional term name (a character string) specifying which term to plot.
#'   If omitted, the first non-constant term is plotted.
#' @param title Optional title for the plot. If not provided, a default title is generated.
#' @param xlab Label for the x-axis (default: "Time").
#' @param ylab Label for the y-axis (default: "Design Matrix Value").
#' @param line_size Numeric value for line thickness (default: 1).
#' @param color_palette A palette name for the line colors (default: "Set1").
#' @param ... Additional arguments passed to ggplot2::geom_line.
#'
#' @importFrom ggplot2 ggplot aes_string geom_line facet_wrap labs theme_minimal scale_color_brewer
#' @importFrom tidyr gather
#' @autoglobal
#' @export
plot.baseline_model <- function(x, term_name = NULL, title = NULL, 
                                xlab = "Time", ylab = "Design Matrix Value",
                                line_size = 1, color_palette = "Set1", ...) {
  # Extract terms and term names from the baseline model.
  all_terms <- terms(x)
  term_names <- names(all_terms)
  
  # Remove constant terms from plotting.
  cidx <- grep("constant", term_names, ignore.case = TRUE)
  if (length(cidx) > 0) {
    all_terms <- all_terms[-cidx]
    term_names <- term_names[-cidx]
  }
  
  # Check if any non-constant terms remain.
  if (length(all_terms) == 0) {
    stop("No non-constant baseline terms available for plotting.")
  }
  
  # Extract the sampling frame.
  sframe <- x$sampling_frame
  if (is.null(sframe$time) || is.null(sframe$blockids)) {
    stop("The sampling_frame in the baseline model must contain 'time' and 'blockids'.")
  }
  
  # For each term, convert its design matrix into a long-format tibble.
  dflist <- lapply(all_terms, function(term) {
    dm <- design_matrix(term)
    dm_tib <- suppressMessages(tibble::as_tibble(dm, .name_repair = "check_unique"))
    dm_tib$.block <- sframe$blockids
    dm_tib$.time <- sframe$time
    tidyr::gather(dm_tib, key = "condition", value = "value", - .time, - .block)
  })
  names(dflist) <- term_names
  
  # Select the term to plot.
  if (is.null(term_name)) {
    dfx <- dflist[[term_names[1]]]
    plot_term <- term_names[1]
  } else {
    if (!term_name %in% term_names) {
      stop("The specified term_name is not found in the baseline model terms.")
    }
    dfx <- dflist[[term_name]]
    plot_term <- term_name
  }
  
  # Create the ggplot.
  p <- ggplot2::ggplot(dfx, ggplot2::aes_string(x = ".time", y = "value", colour = "condition")) +
    ggplot2::geom_line(size = line_size, ...) +
    ggplot2::facet_wrap(~ .block, ncol = 1) +
    ggplot2::labs(title = if (!is.null(title)) title else paste("Baseline Model:", plot_term),
                  x = xlab, y = ylab, colour = "Condition") +
    ggplot2::scale_color_brewer(palette = color_palette) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   axis.title = ggplot2::element_text(face = "bold"))
  
  p
}

#' correlation_map.baseline_model
#'
#' @description
#' Generates a correlation heatmap of the columns in a \code{baseline_model}'s design matrix.
#'
#' @param x A \code{baseline_model}.
#' @inheritParams correlation_map.event_model
#' @export
correlation_map.baseline_model <- function(x,
                                           method          = c("pearson", "spearman"),
                                           half_matrix     = FALSE,
                                           absolute_limits = TRUE,
                                           ...) {
  DM <- as.matrix(design_matrix(x))
  .correlation_map_common(DM, method=method, half_matrix=half_matrix,
                          absolute_limits=absolute_limits, ...)
}


#' Heatmap visualization of the baseline_model design matrix
#'
#' @description
#' Produces a heatmap of all columns in the design matrix for a `baseline_model` object,
#' with rows corresponding to scans and columns corresponding to regressors. By default,
#' it draws horizontal lines separating runs (blocks), and rotates the column labels diagonally.
#'
#' @param x A `baseline_model` object.
#' @param block_separators Logical; if `TRUE`, draw white horizontal lines between blocks.
#' @param rotate_x_text Logical; if `TRUE`, rotate x-axis labels by 45 degrees.
#' @param fill_midpoint Numeric or `NULL`; if not `NULL`, used as the `midpoint` in
#'   [ggplot2::scale_fill_gradient2()] to center the color scale (for example at 0).
#' @param fill_limits Numeric vector of length 2 or `NULL`; passed to the fill scale
#'   `limits` argument. Can clip or expand the color range.
#' @param ... Additional arguments forwarded to [ggplot2::geom_tile()].
#'
#' @import ggplot2
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @return A ggplot2 plot object.
#' @export
design_map.baseline_model <- function(x,
                                      block_separators = TRUE,
                                      rotate_x_text    = TRUE,
                                      fill_midpoint    = NULL,
                                      fill_limits      = NULL,
                                      ...) {
  # 1) Extract the design matrix
  DM <- design_matrix(x)
  n_scans <- nrow(DM)
  
  # 2) Convert to long format
  df_long <- tibble::as_tibble(DM, .name_repair = "unique")
  df_long$scan_number <- seq_len(n_scans)
  df_long <- tidyr::pivot_longer(
    df_long,
    cols      = -scan_number,
    names_to  = "Regressor",
    values_to = "Value"
  )
  
  # 3) Build the base ggplot
  plt <- ggplot(df_long, aes(x = Regressor, y = scan_number, fill = Value)) +
    geom_tile(...)
  
  # 4) Reverse the y-axis so that scan #1 is at top
  plt <- plt + scale_y_reverse()
  
  # 5) Decide on color scale
  #    - If fill_midpoint is set, use scale_fill_gradient2 to center the scale
  #    - Otherwise, use a default 3-color gradient
  if (is.null(fill_midpoint)) {
    plt <- plt + scale_fill_gradientn(
      colours = c("navy", "white", "firebrick"),
      limits  = fill_limits
    )
  } else {
    plt <- plt + scale_fill_gradient2(
      midpoint = fill_midpoint,
      low      = "navy",
      mid      = "white",
      high     = "firebrick",
      limits   = fill_limits
    )
  }
  
  # 6) Optionally draw white horizontal lines to separate blocks
  sframe <- x$sampling_frame
  if (block_separators && !is.null(sframe$blockids)) {
    block_ids  <- sframe$blockids
    run_info   <- rle(block_ids)             # lengths of each block
    row_breaks <- cumsum(run_info$lengths)   # boundary after each block
    ncols      <- ncol(DM)
    
    # Add horizontal lines
    for (rb in row_breaks[-length(row_breaks)]) {
      plt <- plt + 
        annotate("segment",
                 x    = 0.5,
                 xend = ncols + 0.5,
                 y    = rb + 0.5,
                 yend = rb + 0.5,
                 color = "white", size = 1)
    }
  }
  
  # 7) Clean up theme
  plt <- plt + 
    theme_minimal(base_size = 14) +
    labs(x = "Regressors", y = "Scan Number", fill = "Value") +
    theme(
      panel.grid  = element_blank(),
      axis.text.x = if (rotate_x_text) element_text(angle = 45, hjust = 1) else element_text()
    )
  
  plt
}


## ============================================================================
## Section 5: End of File
###############################################################################
