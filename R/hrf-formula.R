#' @importFrom rlang enquos enexpr syms is_formula is_quosure is_call as_label %||% quo_get_expr is_symbol
#' @importFrom assertthat assert_that
#' @importFrom stats setNames


#' @keywords internal
#' @noRd
make_hrf <- function(basis, lag, nbasis=1) {
  if (!is.numeric(lag) || length(lag) > 1) {
    stop("hrf: 'lag' must be a numeric scalar")
  }
  
  if (is.character(basis)) {
    # Map character names to HRF constants
    hrf_map <- list(
      "spmg1" = fmrihrf::HRF_SPMG1,
      "spmg2" = fmrihrf::HRF_SPMG2,
      "spmg3" = fmrihrf::HRF_SPMG3,
      "gamma" = fmrihrf::HRF_GAMMA,
      "gaussian" = fmrihrf::HRF_GAUSSIAN,
      "bspline" = fmrihrf::HRF_BSPLINE,
      "fir" = fmrihrf::HRF_FIR
    )
    
    if (basis %in% names(hrf_map)) {
      base_hrf_obj <- hrf_map[[basis]]
      # Check if nbasis needs to be customized for bases that support it
      if (basis == "bspline" && nbasis != fmrihrf::nbasis(base_hrf_obj)) {
        # Create custom bspline HRF function with specified nbasis
        span <- attr(base_hrf_obj, "span")
        if (is.null(span)) span <- 24
        degree <- attr(base_hrf_obj, "params")$degree
        if (is.null(degree)) degree <- 3
        custom_bspline_fn <- function(t) {
          fmrihrf::hrf_bspline(t, span = span, N = nbasis, degree = degree)
        }
        base_hrf_obj <- fmrihrf::as_hrf(custom_bspline_fn, name = "bspline", nbasis = nbasis, span = span)
      } else if (basis == "fir" && nbasis != fmrihrf::nbasis(base_hrf_obj)) {
        # Create custom FIR HRF with specified nbasis
        span <- attr(base_hrf_obj, "span")
        if (is.null(span)) span <- 24
        custom_fir_fn <- function(t) {
          fmrihrf::hrf_tent_generator(span = span, nbasis = nbasis)(t)
        }
                base_hrf_obj <- fmrihrf::as_hrf(custom_fir_fn, name = "fir", nbasis = nbasis, span = span)
      }
    } else if (basis == "fourier") {
      # Create Fourier HRF with specified nbasis
      span <- 24  # Default span for Fourier
      custom_fourier_fn <- function(t) {
        fmrihrf::hrf_fourier(t, span = span, nbasis = nbasis)
      }
      base_hrf_obj <- fmrihrf::as_hrf(custom_fourier_fn, name = "fourier", nbasis = nbasis, span = span)
    } else {
      # For unknown basis names, just use the default (this will likely fail, but matches old behavior)
      stop("Unknown HRF basis name: ", basis, ". Available options: ", paste(c(names(hrf_map), "fourier"), collapse = ", "))
    }
    # Apply lag using gen_hrf
    final_hrf <- fmrihrf::gen_hrf(base_hrf_obj, lag = lag)

  } else if (inherits(basis, "HRF")) {
    # If it's already an HRF object, apply lag using gen_hrf
    final_hrf <- fmrihrf::gen_hrf(basis, lag = lag)
    
  } else if (is.function(basis)) {
    # If it's a raw function, gen_hrf will handle conversion via as_hrf and apply lag
    final_hrf <- fmrihrf::gen_hrf(basis, lag = lag)

  } else {
    stop("invalid basis function: must be 1) character string indicating hrf type, e.g. 'gamma' 2) a function or 3) an object of class 'HRF': ", basis)
  }
  
  return(final_hrf)
}

#### TODO character variables need an "as.factor"

#' hemodynamic regressor specification function for model formulas.
#' 
#' This function is to be used in formulas for fitting functions, e.g. onsets ~ hrf(fac1,fac2) ...
#' It captures the variables/expressions provided and packages them with HRF/contrast 
#' information into an `hrfspec` object, which is then processed by `event_model`.
#' 
#' @param ... One or more variable names (bare or character) or expressions involving variables 
#'            present in the `data` argument of `event_model`.
#' @param basis the impulse response function or the name of a pre-supplied function, 
#'        one of: "gamma", "spmg1", "spmg2", "spmg3", "bspline", "gaussian", "tent", "bs". 
#'        Can also be an `HRF` object.
#' @param onsets optional onsets override. If missing, onsets will be taken from the LHS of the main model formula.
#' @param durations optional durations override. If missing, durations argument from `event_model` is used.
#' @param prefix a character string that is prepended to the variable names and used to identify the term. 
#'               Can be used to disambiguate two \code{hrf} terms with the same variable(s) but different onsets or basis functions.
#' @param subset an expression indicating the subset of 'onsets' to keep.
#' @param precision sampling precision in seconds.
#' @param nbasis number of basis functions -- only used for hemodynamic response functions (e.g. bspline) that take a variable number of bases.
#' @param contrasts one or more \code{contrast_spec} objects created with the \code{contrast}, `pair_contrast` etc. functions. 
#'                  Must be NULL, a single contrast spec, or a *named* list of contrast specs.
#' @param id a unique \code{character} identifier used to refer to term, otherwise will be determined from variable names.
#' @param name Optional human-readable name for the term.
#' @param lag a temporal offset in seconds which is added to onset before convolution
#' @param summate whether impulse amplitudes sum up when duration is greater than 0.
#' @examples 
#' 
#' ## 'hrf' is typically used in the context of \code{formula}s passed to `event_model`.
#' 
#' # Simple model with one factor
#' form1 <- onsets ~ hrf(condition, basis="spmg1")
#' 
#' # Model with factor and continuous modulator, using default SPMG1 for both terms
#' form2 <- onsets ~ hrf(condition) + hrf(RT)
#' 
#' # Model with interaction term and SPMG3 basis
#' form3 <- onsets ~ hrf(condition, RT, basis="spmg3")
#' 
#' # Model with an expression and contrasts
#' library(rlang)
#' con1 <- pair_contrast(~ condition == "A", ~ condition == "B", name="AvB")
#' form4 <- onsets ~ hrf(condition, Poly(RT, 2), contrasts=con1)
#' 
#' @export
#' @importFrom rlang enquos enexpr syms is_formula is_quosure is_call as_label %||%
#' @return an \code{hrfspec} instance
hrf <- function(..., basis="spmg1", onsets=NULL, durations=NULL, prefix=NULL, subset=NULL, precision=.3, 
                nbasis=1, contrasts=NULL, id=NULL, name=NULL, lag=0, summate=TRUE) {
  
  vars <- rlang::enquos(...) # Capture variables/expressions as quosures
  
  # --- Handle special named arguments within ... ---
  var_names <- names(vars)
  # Remove any quosures whose name matches hrf formals (these are control args, not vars)
  hrf_formals <- names(formals(hrf))
  var_indices <- ! (var_names %in% hrf_formals)
  vars <- vars[var_indices]
  
  if (length(vars) == 0) {
      stop("`hrf` must have at least one variable or expression specified in `...`")
  }
  
  # --- Determine Term ID/Name (Prioritize id, then name, then auto) --- 
  final_id <- id # Prioritize explicit id
  if (is.null(final_id) && !is.null(name)) { # Use name if id is missing
      final_id <- name
  }
  # final_id remains NULL if neither id nor name provided, hrfspec will generate default
  
  # --- Check contrasts argument --- 
  if (!is.null(contrasts)) {
    if (!is.list(contrasts)) {
        if (inherits(contrasts, "contrast_spec")) {
            contrast_list <- list(contrasts)
            cname <- contrasts$name %||% "contrast1"
            names(contrast_list) <- cname
            contrasts <- contrast_list
        } else {
             stop("`contrasts=` argument must be a single contrast_spec object or a list of them.\n",
                  "  Detected type: ", class(contrasts)[1], "\n",
                  "  Hint: Use functions like contrast(), pair_contrast(), contrast_set() to define contrasts.", 
                  call. = FALSE)
        }
    } else {
        is_spec <- vapply(contrasts, inherits, logical(1), "contrast_spec")
        if (!all(is_spec)) {
             stop("If `contrasts=` is a list, all elements must be contrast_spec objects.\n",
                  "  Hint: Use functions like contrast(), pair_contrast(), contrast_set() to define contrasts.", 
                  call. = FALSE)
        }
        if (is.null(names(contrasts)) || any(names(contrasts) == "")) {
            cnames <- sapply(contrasts, function(cs) cs$name %||% paste0("contrast", which(sapply(contrasts, identical, cs))))
             if (any(duplicated(cnames))) {
                stop("If `contrasts=` is a list, it must be a *named* list (or auto-naming must yield unique names).", call. = FALSE)
            }
            names(contrasts) <- cnames
        }
    }
  }
  # -----------------------------
  
  basis_obj <- make_hrf(basis, lag, nbasis=nbasis)
  
  # Call the internal constructor, passing quosures directly
  ret <- hrfspec(
    vars = vars, # Pass list of quosures captured by enquos
    basis = basis_obj,         
    onsets = onsets,       
    durations = durations, 
    prefix = prefix,       
    subset = rlang::enexpr(subset), # Capture subset expr unevaluated
    precision = precision,
    contrasts = contrasts, ## Pass validated list of contrast specs
    summate = summate,
    id = final_id # Pass the determined ID 
    )

  ret
}


#' Internal constructor for hrfspec objects
#' 
#' Creates the hrfspec list structure. Called by `hrf()`.
#' Generates termname and label from the input variables/expressions.
#'
#' @param vars List of quosures representing variables/expressions.
#' @param label Optional label for the term (if NULL, generated automatically).
#' @param basis An `HRF` object.
#' @param ... Other arguments passed from `hrf()` (onsets, durations, prefix, subset, etc.)
#' @return An `hrfspec` instance (list with class `hrfspec`).
#' @importFrom rlang as_label is_symbol is_call quo_get_expr
#' @noRd
#' @keywords internal
hrfspec <- function(vars, label=NULL, basis=fmrihrf::HRF_SPMG1, ...) {

  assert_that(inherits(basis, "HRF"))
  
  # Generate varnames and termname from quosures
  var_labels <- sapply(vars, rlang::as_label)
  # Simple heuristic for names: use symbol directly, otherwise make.names on label
  varnames <- sapply(vars, function(q) {
       expr <- rlang::quo_get_expr(q)
       if (rlang::is_symbol(expr)) as.character(expr) else make.names(rlang::as_label(expr))
  })
  
  termname <- paste0(varnames, collapse=":")
  
  # Generate label if not provided
  if (is.null(label)) {
      label <- paste0("hrf(", paste0(var_labels, collapse=","), ")")
  }
  
  # Capture other arguments passed via ... from hrf()
  other_args <- list(...)
  
  ret <- list(
    name = termname, 
    label = label,
    id = other_args$id,
    vars = vars,
    varnames = varnames, 
    hrf = basis,
    onsets = other_args$onsets,
    durations = other_args$durations,
    prefix = other_args$prefix,
    subset = other_args$subset, # Should be an expression
    precision = other_args$precision %||% 0.3,
    contrasts = other_args$contrasts, # Already validated list or NULL
    summate = other_args$summate %||% TRUE
    # data_env is NOT stored here, added later during evaluation
  )
  
  # If user supplied an explicit id (or via name= alias), use it as the primary term name
  if (!is.null(other_args$id)) {
      ret$name <- other_args$id
  }
  
  class(ret) <- c("hrfspec", "list")
  ret
}



#' @export
nbasis.hrfspec <- function(x, ...) {
  fmrihrf::nbasis(x$hrf)
}

#' @export
construct.hrfspec <- function(x, model_spec, ...) {
  ons <- if (!is.null(x$onsets)) x$onsets else model_spec$onsets
  et <- construct_event_term(x, model_spec)
  
  # Set the term_tag attribute on the event_term before returning
  # This ensures that when convolution happens later, the correct column names are generated
  term_tag <- x$id %||% x$name
  attr(et, "term_tag") <- term_tag
  
  # DON'T convolve here - let build_event_model_design_matrix handle convolution
  # Just return the event_term with hrfspec attached
  # The hrfspec is already attached in construct_event_term
  
  # Handle add_sum flag if present in the hrfspec (set by trialwise)
  # Store this as an attribute for later processing during convolution
  if (isTRUE(x$add_sum)) {
    attr(et, "add_sum") <- TRUE
    attr(et, "add_sum_label") <- x$id %||% x$name
  }

  # Return the event_term directly, not a convolved_term
  et
}

#' evaluate.hrfspec
#'
#' This function evaluates a hemodynamic response function (HRF) specified by an hrfspec object for a given set of time points (grid) and other parameters.
#' It is a wrapper function that calls the evaluate.HRF function with the HRF function contained in the hrfspec object.
#'
#' @param x The hrfspec object containing the HRF function.
#' @param grid A vector of time points.
#' @param amplitude The scaling value for the event (default: 1).
#' @param duration The duration of the event (default: 0).
#' @param precision The temporal resolution used for computing summed responses when duration > 0 (default: 0.1).
#' @param ... Additional arguments to be passed to the evaluate.HRF function.
#' @return A vector of HRF values at the specified time points.
#' @noRd
evaluate.hrfspec <- function(x, grid, amplitude=1, duration=0, precision=.1, ...) {
  fmrihrf::evaluate(x$hrf, grid, amplitude, duration, precision)
}

#' trialwise
#'
#' Generate one regressor per trial (plus an optional grand-mean column)
#' by delegating everything to `hrf()`.
#'
#' Use it **only on the RHS** of an event-model formula:
#'
#'     onset ~ trialwise(basis = "spmg1", add_sum = TRUE)
#'
#' @param basis,lag,nbasis Passed straight to `hrf()`.
#' @param add_sum If TRUE, append a column that is the average of all
#'                trialwise columns (useful as a conventional main effect).
#' @param label Term label / prefix for the generated columns.
#' @export
trialwise <- function(basis   = "spmg1",
                      lag     = 0,
                      nbasis  = 1,
                      add_sum = FALSE,
                      label   = "trial") {

  # Create an expression that will evaluate .trial_factor(length(onsets)) 
  # when the onsets variable is available at evaluation time
  trial_expr <- rlang::expr(.trial_factor(length(onsets)))

  term <- hrf(!!trial_expr,  # Use !! to inject the expression
              basis = basis,
              lag   = lag,
              nbasis = nbasis,
              id     = label) # Use id argument for naming

  term$add_sum <- isTRUE(add_sum) # flag for construct() to act upon
  term
}

#' Internal helper for generating trial factors
#'
#' @param n Length of the factor to generate.
#' @keywords internal
.trial_factor <- function(n) {
  pad <- nchar(as.character(n))
  factor(sprintf(paste0("%0", pad, "d"), seq_len(n)))
}
