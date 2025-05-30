#' Turn any function into an HRF object
#'
#' This is the core constructor for creating HRF objects in the refactored system.
#' It takes a function `f(t)` and attaches standard HRF attributes.
#'
#' @param f The function to be turned into an HRF object. It must accept a single argument `t` (time).
#' @param name The name for the HRF object. Defaults to the deparsed name of `f`.
#' @param nbasis The number of basis functions represented by `f`. Defaults to 1L.
#' @param span The nominal time span (duration in seconds) of the HRF. Defaults to 24.
#' @param params A named list of parameters associated with the HRF function `f`. Defaults to an empty list.
#' @param ... Additional arguments to pass to the function f.
#' @return A new HRF object.
#' @keywords internal
#' @export
as_hrf <- function(f, name = deparse(substitute(f)), nbasis = 1L, span = 24,
                   params = list()) {
  assertthat::assert_that(is.function(f))
  assertthat::assert_that(is.character(name), length(name) == 1)
  assertthat::assert_that(is.numeric(nbasis), length(nbasis) == 1)
  assertthat::assert_that(is.numeric(span), length(span) == 1)
  assertthat::assert_that(is.list(params))

  structure(
    f,
    class        = c("HRF", "function"),
    name         = name,
    nbasis       = as.integer(nbasis),
    span         = span,
    param_names  = names(params),
    params       = params
  )
}


#' Bind HRFs into a Basis Set
#'
#' Combines multiple HRF objects into a single multi-basis HRF object.
#' The resulting function evaluates each input HRF at time `t` and returns the results column-bound together.
#'
#' @param ... One or more HRF objects created by `as_hrf` or other HRF constructors/decorators.
#'
#' @return A new HRF object representing the combined basis set.
#'
#' @keywords internal
#' @export
#' @importFrom assertthat assert_that
bind_basis <- function(...) {
  xs <- list(...)
  assertthat::assert_that(length(xs) > 0, msg = "bind_basis requires at least one HRF object.")
  assertthat::assert_that(all(sapply(xs, inherits, "HRF")), msg = "All inputs to bind_basis must be HRF objects.")

  # Handle single HRF case explicitly
  if (length(xs) == 1) {
    return(xs[[1]])
  }

  # Calculate combined attributes
  combined_nbasis <- sum(vapply(xs, attr, 0L, "nbasis"))
  combined_span <- max(vapply(xs, attr, 0, "span"))
  combined_name <- paste(sapply(xs, attr, "name"), collapse = " + ")

  # Create the combined function
  combined_func <- function(t) {
    do.call(cbind, lapply(xs, function(f) f(t)))
  }

  # Use as_hrf to create the new HRF object
  as_hrf(
    f = combined_func,
    name = combined_name,
    nbasis = combined_nbasis,
    span = combined_span,
    params = list() # Params usually don't combine meaningfully
  )
}


#' Construct an HRF Instance using Decorators
#' 
#' @description
#' `gen_hrf` takes a base HRF function or object and applies optional lag,
#' blocking, and normalization decorators based on arguments.
#'
#' @param hrf A function `f(t)` or an existing `HRF` object.
#' @param lag Optional lag in seconds. If non-zero, applies `lag_hrf`.
#' @param width Optional block width in seconds. If non-zero, applies `block_hrf`.
#' @param precision Sampling precision for block convolution (passed to `block_hrf`). Default is 0.1.
#' @param summate Whether to summate within blocks (passed to `block_hrf`). Default is TRUE.
#' @param normalize If TRUE, applies `normalise_hrf` at the end. Default is FALSE.
#' @param name Optional name for the *final* HRF object. If NULL (default), a name is generated based on the base HRF and applied decorators.
#' @param span Optional span for the *final* HRF object. If NULL (default), the span is determined by the base HRF and decorators.
#' @param ... Extra arguments passed to the *base* HRF function if `hrf` is a function.
#'
#' @return A final `HRF` object, potentially modified by decorators.
#' 
#' @examples 
#' # Lagged SPMG1
#' grf_lag <- gen_hrf(HRF_SPMG1, lag=3)
#' # Blocked Gaussian
#' grf_block <- gen_hrf(hrf_gaussian, width=5, precision=0.2)
#' # Lagged and Blocked, then Normalized
#' grf_both_norm <- gen_hrf(HRF_SPMG1, lag=2, width=4, normalize=TRUE)
#'
#' @export
gen_hrf <- function(hrf, lag=0, width=0, precision=.1, 
                    summate=TRUE, normalize=FALSE, name=NULL, span=NULL, ...) {

  # 1. Ensure we start with an HRF object
  if (is.function(hrf) && !inherits(hrf, "HRF")) {
    # If it's a plain function, convert it using as_hrf
    # Determine nbasis by evaluating the function
    test_t <- 1:10 # A small sample range
    test_val <- try(hrf(test_t, ...), silent = TRUE)
    determined_nbasis <- if (!inherits(test_val, "try-error") && !is.null(test_val)) {
      if (is.matrix(test_val)) ncol(test_val) else 1L
    } else {
      warning(paste("Could not determine nbasis for function", deparse(substitute(hrf)), "- defaulting to 1. Evaluation failed."))
      1L
    }
    
    # Pass extra args (...) here if they are meant for the base function construction
    base_hrf <- as_hrf(f = function(t) hrf(t, ...),
                       name = deparse(substitute(hrf)),
                       nbasis = determined_nbasis) # Pass determined nbasis
                       # Let as_hrf determine default span, params
  } else if (inherits(hrf, "HRF")) {
    # If already an HRF object, use it directly
    base_hrf <- hrf
    if (length(list(...)) > 0) {
      warning("Ignoring extra arguments (...) because 'hrf' is already an HRF object.")
    }
  } else {
    stop("'hrf' must be a function or an HRF object.")
  }

  # Apply decorators conditionally
  decorated_hrf <- base_hrf

  # Apply width decorator first if needed
  if (width != 0) {
    # Check positivity *before* applying
    stopifnot(width > 0)
    # Note: block_hrf handles normalize=FALSE internally by default
    decorated_hrf <- block_hrf(decorated_hrf, width=width, precision=precision,
                               summate=summate, normalize=FALSE)
  }

  # Apply lag decorator if needed
  if (lag != 0) {
    decorated_hrf <- lag_hrf(decorated_hrf, lag=lag)
  }

  # Apply normalization decorator last if needed
  if (normalize) {
    decorated_hrf <- normalise_hrf(decorated_hrf)
  }

  # Override name and span if provided by user
  if (!is.null(name)) {
    attr(decorated_hrf, "name") <- name
  }
  if (!is.null(span)) {
    attr(decorated_hrf, "span") <- span
  }

  # Return the final (potentially decorated) HRF object
  return(decorated_hrf)
}


#' Generate an Empirical Hemodynamic Response Function
#' 
#' @description
#' `empirical_hrf` generates an empirical HRF using provided time points and values.
#' 
#' @param t Time points.
#' @param y Values of HRF at time `t[i]`.
#' @param name Name of the generated HRF.
#' @return An instance of type `HRF`.
#' @export
empirical_hrf <- function(t, y, name = "empirical_hrf") {
  as_hrf(stats::approxfun(t, y, yright = 0, yleft = 0),
         name = name, nbasis = 1L, span = max(t, na.rm = TRUE))
}

#' @export
#' @rdname empirical_hrf
#' @keywords internal
gen_empirical_hrf <- function(...) {
  .Deprecated("empirical_hrf")
  empirical_hrf(...)
}


#' Generate an HRF Basis Set
#' 
#' @description
#' `hrf_set` constructs an HRF basis set from one or more component HRF objects.
#'
#' @param ... One or more HRF objects.
#' @param name The name for the combined HRF set.
#' @return A combined HRF object.
#' @export
hrf_set <- function(..., name = "hrf_set") {
  combined_hrf <- bind_basis(...)
  attr(combined_hrf, "name") <- name
  combined_hrf
}

#' @export
#' @rdname hrf_set
#' @keywords internal
gen_hrf_set <- function(...) {
  .Deprecated("hrf_set")
  hrf_set(...)
}


#' Generate an HRF library from a parameter grid
#'
#' @description
#' `hrf_library` applies a base HRF generating function to each row of a parameter grid.
#'
#' @param fun A function that generates an HRF, given a set of parameters.
#' @param pgrid A data frame where each row is a set of parameters.
#' @param ... Additional arguments passed to `fun`.
#' @return A combined HRF object representing the library.
#' @importFrom purrr pmap partial
#' @export
hrf_library <- function(fun, pgrid, ...) {
  # Ensure fun returns an HRF object
  hrf_list <- purrr::pmap(pgrid, function(...) {
      params <- list(...)
      # Assuming 'fun' is designed to take these params and return an HRF
      # If fun itself is just a base function, we need to wrap it
      # Let's assume fun already produces an HRF or can be wrapped by as_hrf
      # This part might need adjustment based on typical usage of 'fun'
      do.call(fun, c(params, list(...))) # Pass original dots as well?
      # Safest might be if 'fun' is expected to return an HRF object directly
      # Example: fun = function(lag) HRF_SPMG1 |> lag_hrf(lag)
  })
  # Bind the generated HRFs
  do.call(bind_basis, hrf_list)
}

#' @export
#' @rdname hrf_library
#' @keywords internal
gen_hrf_library <- function(...) {
  .Deprecated("hrf_library")
  hrf_library(...)
}


#' HRF Constructor Function
#'
#' The `HRF` function creates an object representing a hemodynamic response function (HRF). It is a class constructor for HRFs.
#'
#' @param fun A function representing the hemodynamic response, mapping from time to BOLD response.
#' @param name A string specifying the name of the function.
#' @param nbasis An integer representing the number of basis functions, e.g., the columnar dimension of the HRF. Default is 1.
#' @param span A numeric value representing the span in seconds of the HRF. Default is 24.
#' @param param_names A character vector containing the names of the parameters for the HRF function.
#'
#' @return An HRF object with the specified properties.
#'
#' @details
#' The package provides several pre-defined HRF types that can be used in modeling fMRI responses:
#'
#' **Canonical HRFs:**
#' * `"spmg1"` or `HRF_SPMG1`: SPM's canonical HRF (single basis function)
#' * `"spmg2"` or `HRF_SPMG2`: SPM canonical + temporal derivative (2 basis functions)
#' * `"spmg3"` or `HRF_SPMG3`: SPM canonical + temporal and dispersion derivatives (3 basis functions)
#' * `"gaussian"` or `HRF_GAUSSIAN`: Gaussian-shaped HRF with peak around 5-6s
#' * `"gamma"` or `HRF_GAMMA`: Gamma function-based HRF with longer tail
#'
#' **Flexible basis sets:**
#' * `"bspline"` or `"bs"` or `HRF_BSPLINE`: B-spline basis for flexible HRF modeling
#' * `"tent"`: Tent (triangular) basis functions for flexible HRF modeling
#' * `"daguerre"` or `HRF_DAGUERRE`: Daguerre basis functions
#'
#' To see a complete list of available HRF types with details, use the `list_available_hrfs()` function.
#'
#' @examples
#' hrf <- HRF(hrf_gamma, "gamma", nbasis=1, param_names=c("shape", "rate"))
#' resp <- evaluate(hrf, seq(0, 24, by=1))
#'
#' # List all available HRF types
#' list_available_hrfs(details = TRUE)
#'
#' @export
#' @rdname HRF-class
HRF <- function(fun, name, nbasis=1, span=24, param_names=NULL) {
  vals <- try(fun(seq(0, span)), silent = TRUE)

  peak <- if (!inherits(vals, "try-error") && !is.null(vals)) {
    if (nbasis == 1) {
      max(vals, na.rm = TRUE)
    } else if (is.matrix(vals)) {
      max(apply(vals, 2, max, na.rm = TRUE))
    } else {
      NA # Unable to determine peak
    }
  } else {
    NA # Error during evaluation or null result
  }

  scale_factor <- if (!is.na(peak) && peak != 0) {
    1 / peak
  } else {
    NA # Cannot compute scale_factor if peak is NA or zero
  }
  
  structure(fun, name=name, 
            nbasis=as.integer(nbasis), 
            span=span,
            param_names=param_names, 
            scale_factor=scale_factor, 
            class=c("HRF", "function"))
  
}

#' @rdname nbasis
#' @export
nbasis.HRF <- function(x,...) attr(x, "nbasis")



#' @keywords internal
#' @noRd
makeDeriv <- function(HRF, n=1) {
  #with_package("numDeriv")
  if (n == 1) {
    function(t) numDeriv::grad(HRF, t)
  } else {
    Recall(function(t) numDeriv::grad(HRF,t), n-1)
  }
}


#' Generate a Lagged HRF Function
#'
#' @description
#' The `gen_hrf_lagged` function takes an HRF function and applies a specified lag to it. This can be useful for modeling time-delayed hemodynamic responses.
#'
#' @param hrf A function representing the underlying HRF to be shifted.
#' @param lag A numeric value specifying the lag or delay in seconds to apply to the HRF. This can also be a vector of lags, in which case the function returns an HRF set.
#' @param normalize A logical value indicating whether to rescale the output so that the maximum absolute value is 1. Defaults to `FALSE`.
#' @param ... Extra arguments supplied to the `hrf` function.
#'
#' @return A function representing the lagged HRF. If `lag` is a vector of lags, the function returns an HRF set.
#' @family gen_hrf
#' @examples
#' \donttest{
#' hrf_lag5 <- gen_hrf_lagged(HRF_SPMG1, lag=5)
#' hrf_lag5(0:20)
#' }
#'
#' @export
gen_hrf_lagged <- function(hrf, lag=2, normalize=FALSE, ...) {
  force(hrf)
  # TODO deal with nbasis arg in ...
  if (length(lag)>1) {
    do.call(gen_hrf_set, lapply(lag, function(l) gen_hrf_lagged(hrf, l,...)))
  } else {
    function(t) {
      ret <- hrf(t-lag,...)
      if (normalize) {
        ret <- ret/max(abs(ret))
      } 
      
      ret
    }
  }
}

#' @export
#' @describeIn gen_hrf_lagged alias for gen_hrf_lagged
#' @family gen_hrf
#' @return an lagged hrf function
hrf_lagged <- gen_hrf_lagged


#' Generate a Blocked HRF Function
#'
#' @description
#' The `gen_hrf_blocked` function creates a blocked HRF by convolving the input HRF with a boxcar function. This can be used to model block designs in fMRI analysis.
#'
#' @param hrf A function representing the hemodynamic response function. Default is `hrf_gaussian`.
#' @param width A numeric value specifying the width of the block in seconds. Default is 5.
#' @param precision A numeric value specifying the sampling resolution in seconds. Default is 0.1.
#' @param half_life A numeric value specifying the half-life of the exponential decay function, used to model response attenuation. Default is `Inf`, which means no decay.
#' @param summate A logical value indicating whether to allow each impulse response function to "add" up. Default is `TRUE`.
#' @param normalize A logical value indicating whether to rescale the output so that the peak of the output is 1. Default is `FALSE`.
#' @param ... Extra arguments passed to the HRF function.
#' @family gen_hrf
#'
#' @return A \code{function} representing the blocked HRF.
#'
#' @importFrom purrr partial
#' @export
gen_hrf_blocked <- function(hrf=hrf_gaussian, width=5, precision=.1, 
                            half_life=Inf, summate=TRUE, normalize=FALSE, ...) {
  force(hrf)
  purrr::partial(convolve_block, hrf=hrf, width=width, 
                 precision=precision, half_life=half_life, 
                 summate=summate, normalize=normalize, ...)
}

#' @export
#' @aliases gen_hrf_blocked
#' @describeIn gen_hrf_blocked alias for gen_hrf_blocked
#' @return A \code{function} representing the blocked HRF.
hrf_blocked <- gen_hrf_blocked





#' Soft-threshold function
#'
#' This function applies soft-thresholding to the input values, setting values below the threshold to zero
#' and shrinking the remaining values by the threshold amount.
#'
#' @param x A numeric vector of input values
#' @param threshold A non-negative threshold value for the soft-thresholding operation
#'
#' @return A numeric vector with the soft-thresholded values
#'
#' @noRd
#' @keywords internal
soft_threshold <- function(x, threshold) {
  if (threshold < 0) {
    stop("Threshold value should be non-negative.")
  }

  sign(x) * pmax(0, abs(x) - threshold)
}



#' List all available hemodynamic response functions (HRFs)
#'
#' @description
#' Reads the internal HRF registry to list available HRF types.
#'
#' @param details Logical; if TRUE, attempt to add descriptions (basic for now).
#' @return A data frame with columns: name, type (object/generator), nbasis_default.
#' @export
list_available_hrfs <- function(details = FALSE) {
  # Get names directly from the registry
  hrf_names <- names(HRF_REGISTRY)
  
  # Determine type and default nbasis by inspecting registry entries
  hrf_info <- lapply(hrf_names, function(name) {
    entry <- HRF_REGISTRY[[name]]
    type <- if (inherits(entry, 'HRF')) "object" else if (is.function(entry)) "generator" else "unknown"
    
    nbasis_default <- NA
    if (type == "object") {
      nbasis_default <- tryCatch(nbasis(entry), error = function(e) NA)
    } else if (type == "generator") {
      fmls <- formals(entry)
      if ("nbasis" %in% names(fmls)) {
        nb_val <- fmls$nbasis
        if(is.numeric(nb_val)) nbasis_default <- nb_val 
      } 
      if(is.na(nbasis_default)) nbasis_default <- "variable"
    }
    
    # Check if this name is an alias (points to the same object/func as another primary name)
    is_alias <- FALSE
    if (type == "object") {
      primary_names <- names(HRF_REGISTRY)[sapply(HRF_REGISTRY, identical, entry)]
      is_alias <- length(primary_names) > 1 && name %in% primary_names[primary_names != name]
    } else if (type == "generator") {
      # More complex for functions, check if it points to the same generator function
      # For now, let's assume aliases only exist for objects, or mark known ones
      is_alias <- name %in% c("gam", "bs")
    }

    list(name = name, type = type, nbasis_default = as.character(nbasis_default), is_alias = is_alias) 
  })
  
  # Combine into a data frame
  hrf_df <- do.call(rbind.data.frame, c(hrf_info, list(stringsAsFactors = FALSE)))
  
  # Add basic descriptions if requested
  if (details) {
      hrf_df$description <- paste(hrf_df$name, "HRF", 
                                  ifelse(hrf_df$type == "generator", "(generator)", "(object)"),
                                  ifelse(hrf_df$is_alias, "(alias)", ""))
  }
  
  hrf_df
}

# Define Static HRF Objects -----

#' Pre-defined Hemodynamic Response Function Objects
#' 
#' A collection of pre-defined HRF objects for common fMRI analysis scenarios.
#' These objects can be used directly in model specifications or as templates
#' for creating custom HRFs.
#' 
#' @section Canonical HRFs:
#' \describe{
#'   \item{\code{HRF_SPMG1}}{SPM canonical HRF (single basis function)}
#'   \item{\code{HRF_SPMG2}}{SPM canonical HRF with temporal derivative (2 basis functions)}
#'   \item{\code{HRF_SPMG3}}{SPM canonical HRF with temporal and dispersion derivatives (3 basis functions)}
#'   \item{\code{HRF_GAMMA}}{Gamma function-based HRF}
#'   \item{\code{HRF_GAUSSIAN}}{Gaussian function-based HRF}
#' }
#' 
#' @section Flexible Basis Sets:
#' \describe{
#'   \item{\code{HRF_BSPLINE}}{B-spline basis HRF (5 basis functions)}
#'   \item{\code{HRF_FIR}}{Finite Impulse Response (FIR) basis HRF (default 12 basis functions)}
#' }
#' 
#' @section Usage:
#' All HRF objects can be:
#' \itemize{
#'   \item Called as functions with time argument: \code{HRF_SPMG1(t)}
#'   \item Used in model specifications: \code{hrf(condition, basis = HRF_SPMG1)}
#'   \item Evaluated with \code{evaluate()} method
#'   \item Combined with decorators like \code{lag_hrf()} or \code{block_hrf()}
#' }
#' 
#' @param t Numeric vector of time points (in seconds) at which to evaluate the HRF
#' @param P1,P2 Shape parameters for SPM canonical HRF (default: P1=5, P2=15)
#' @param A1 Amplitude parameter for SPM canonical HRF (default: 0.0833)
#' @param shape,rate Parameters for gamma distribution HRF (default: shape=6, rate=1)
#' @param mean,sd Parameters for Gaussian HRF (default: mean=6, sd=2)
#' 
#' @return 
#' When called as functions, return numeric vectors or matrices of HRF values.
#' When used as objects, they are HRF objects with class \code{c("HRF", "function")}.
#' 
#' @examples
#' # Evaluate HRFs at specific time points
#' times <- seq(0, 20, by = 0.5)
#' 
#' # Single basis canonical HRF
#' canonical_response <- HRF_SPMG1(times)
#' plot(times, canonical_response, type = "l", main = "SPM Canonical HRF")
#' 
#' # Multi-basis HRF with derivatives
#' multi_response <- HRF_SPMG3(times)  # Returns 3-column matrix
#' matplot(times, multi_response, type = "l", main = "SPM HRF with Derivatives")
#' 
#' # Gamma and Gaussian HRFs
#' gamma_response <- HRF_GAMMA(times)
#' gaussian_response <- HRF_GAUSSIAN(times)
#' 
#' # Compare different HRF shapes
#' plot(times, canonical_response, type = "l", col = "blue", 
#'      main = "HRF Comparison", ylab = "Response")
#' lines(times, gamma_response, col = "red")
#' lines(times, gaussian_response, col = "green")
#' legend("topright", c("SPM Canonical", "Gamma", "Gaussian"), 
#'        col = c("blue", "red", "green"), lty = 1)
#' 
#' # Use in model specification
#' \dontrun{
#' # In an event model
#' evmodel <- event_model(
#'   onsets ~ hrf(condition, basis = HRF_SPMG1),
#'   data = event_data,
#'   sampling_frame = sframe
#' )
#' 
#' # With multiple basis functions
#' evmodel2 <- event_model(
#'   onsets ~ hrf(condition, basis = HRF_SPMG3),
#'   data = event_data,
#'   sampling_frame = sframe
#' )
#' }
#' 
#' @name HRF_objects
#' @aliases HRF_SPMG1 HRF_SPMG2 HRF_SPMG3 HRF_GAMMA HRF_GAUSSIAN HRF_BSPLINE HRF_FIR
#' @family hrf
#' @seealso \code{\link{evaluate.HRF}}, \code{\link{gen_hrf}}, \code{\link{list_available_hrfs}}
NULL

#' @rdname HRF_objects
#' @export
HRF_GAMMA <- as_hrf(hrf_gamma, name="gamma", params=list(shape=6, rate=1))

#' @rdname HRF_objects
#' @export
HRF_GAUSSIAN <- as_hrf(hrf_gaussian, name="gaussian", params=list(mean=6, sd=2))

#' @rdname HRF_objects
#' @export
HRF_SPMG1 <- as_hrf(hrf_spmg1, name="SPMG1", params=list(P1=5, P2=15, A1=0.0833))

#' @rdname HRF_objects
#' @export
HRF_SPMG2 <- bind_basis(
  as_hrf(hrf_spmg1, name="SPMG1_canonical", params=list(P1=5, P2=15, A1=0.0833)),
  as_hrf(hrf_spmg1_deriv, name="SPMG1_temporal_deriv", params=list(P1=5, P2=15, A1=0.0833))
)
attr(HRF_SPMG2, "name") <- "SPMG2"

#' @rdname HRF_objects
#' @export
HRF_SPMG3 <- bind_basis(
  as_hrf(hrf_spmg1, name="SPMG1_canonical", params=list(P1=5, P2=15, A1=0.0833)),
  as_hrf(hrf_spmg1_deriv, name="SPMG1_temporal_deriv", params=list(P1=5, P2=15, A1=0.0833)),
  as_hrf(hrf_spmg1_second_deriv, name="SPMG1_dispersion_deriv", params=list(P1=5, P2=15, A1=0.0833))
)
attr(HRF_SPMG3, "name") <- "SPMG3"

# Define HRF Generators (Functions returning HRF objects) -----
hrf_bspline_generator <- function(nbasis=5, span=24) {
  # Validate inputs
  if (nbasis < 1) {
    stop("nbasis must be at least 1", call. = FALSE)
  }
  if (span <= 0) {
    stop("span must be positive", call. = FALSE)
  }
  
  # Ensure nbasis is integer
  nbasis <- as.integer(nbasis)
  
  degree <- 3 # Default cubic B-splines
  effective_nbasis <- max(1, nbasis) 
  
  f_bspline <- function(t) {
    valid_t_idx <- t >= 0 & t <= span
    if (!any(valid_t_idx)) {
      return(matrix(0, nrow = length(t), ncol = effective_nbasis))
    }
    
    res_mat <- matrix(0, nrow = length(t), ncol = effective_nbasis)
    
    bs_matrix <- tryCatch({
        splines::bs(t[valid_t_idx], df = effective_nbasis, degree = degree, 
                    Boundary.knots = c(0, span), intercept = FALSE)
    }, error = function(e) {
        warning(sprintf("splines::bs failed for effective_nbasis=%d, span=%d: %s", 
                        effective_nbasis, span, e$message), call. = FALSE)
        NULL 
    })

    if (!is.null(bs_matrix) && ncol(bs_matrix) == effective_nbasis) {
      res_mat[valid_t_idx, ] <- bs_matrix
    } else if (!is.null(bs_matrix)) {
       warning(sprintf("splines::bs returned %d columns, expected %d for effective_nbasis=%d, span=%d. Returning zeros.", 
                      ncol(bs_matrix), effective_nbasis, effective_nbasis, span), call. = FALSE)
    }
    return(res_mat)
  }

  as_hrf(
    f = f_bspline,
    name = "bspline", nbasis = as.integer(effective_nbasis), span = span,
    params = list(nbasis = effective_nbasis, degree = degree, span = span) 
  )
}

hrf_tent_generator <- function(nbasis=5, span=24) {
  as_hrf(
    f = function(t) hrf_bspline(t, span=span, N=nbasis, degree=1), # hrf_bspline from hrf-functions.R
    name="tent", nbasis=as.integer(nbasis), span=span,
    params=list(N=nbasis, degree=1, span=span)
  )
}

hrf_fourier_generator <- function(nbasis=5, span=24) {
  as_hrf(
    f = function(t) hrf_fourier(t, span=span, nbasis=nbasis), # hrf_fourier from hrf-functions.R
    name="fourier", nbasis=as.integer(nbasis), span=span,
    params=list(nbasis=nbasis, span=span)
  )
}

hrf_daguerre_generator <- function(nbasis=3, scale=4) {
  as_hrf(
    f = function(t) daguerre_basis(t, n_basis=nbasis, scale=scale), # daguerre_basis from hrf-functions.R
    name="daguerre", nbasis=as.integer(nbasis), span=24, # Default span, daguerre is time-scaled
    params=list(n_basis=nbasis, scale=scale)
  )
}

hrf_fir_generator <- function(nbasis = 12, span = 24) {
  assertthat::assert_that(
    is.numeric(nbasis) && length(nbasis) == 1 && nbasis >= 1,
    msg = "`nbasis` must be a single positive integer."
  )
  assertthat::assert_that(
    is.numeric(span) && length(span) == 1 && span > 0,
    msg = "`span` must be a single positive number."
  )
  nbasis <- as.integer(nbasis)
  bin_width <- span / nbasis

  f_fir <- function(t) {
    if (!is.numeric(t) || length(t) == 0) {
      return(matrix(0, nrow = 0, ncol = nbasis))
    }
    output_matrix <- matrix(0, nrow = length(t), ncol = nbasis)
    for (i in seq_along(t)) {
      current_t <- t[i]
      if (!is.na(current_t) && current_t >= 0 && current_t < span) {
        bin_index <- if (current_t == 0) 1 else floor(current_t / bin_width) + 1
        bin_index <- min(bin_index, nbasis) 
        output_matrix[i, bin_index] <- 1
      }
    }
    return(output_matrix)
  }

  as_hrf(
    f = f_fir,
    name = "fir",
    nbasis = nbasis,
    span = span,
    params = list(nbasis = nbasis, span = span, bin_width = bin_width)
  )
}

# Define HRF Registry -----

#' @keywords internal
HRF_REGISTRY <- list(
  spmg1    = HRF_SPMG1,
  spmg2    = HRF_SPMG2,
  spmg3    = HRF_SPMG3,
  gamma    = HRF_GAMMA,
  gaussian = HRF_GAUSSIAN,
  bspline  = hrf_bspline_generator,
  tent     = hrf_tent_generator,
  fourier  = hrf_fourier_generator,
  daguerre = hrf_daguerre_generator,
  fir      = hrf_fir_generator,
  lwu      = hrf_lwu 
)

HRF_REGISTRY$gam <- HRF_REGISTRY$gamma
HRF_REGISTRY$bs  <- HRF_REGISTRY$bspline

# getHRF function using the registry (Minimal Version) -----

#' getHRF
#'
#' Retrieves an HRF by name from the registry and applies decorators.
#'
#' @param ... Additional arguments passed to generator functions (e.g., `scale` for daguerre).
#' @return An HRF object.
#' @keywords internal
#' @noRd
getHRF <- function(name = "spmg1", # Default to spmg1
                   nbasis=5, span=24,
                   lag=0, width=0,
                   summate=TRUE, normalize=FALSE, ...) {

  key   <- match.arg(tolower(name), names(HRF_REGISTRY))
  entry <- HRF_REGISTRY[[key]]

  base <- if (inherits(entry, "HRF")) {
            entry # Use pre-defined object
      } else {
            # Call generator, passing nbasis, span, and any relevant ... args
            gen_args <- c(list(nbasis=as.integer(nbasis), span=span), list(...))
            # Only pass args the generator actually accepts
            valid_args <- gen_args[names(gen_args) %in% names(formals(entry))]
            do.call(entry, valid_args)
          }

  # Apply decorators
  if (width != 0) {
      stopifnot(width > 0)
      base <- block_hrf(base, width = width, summate = summate)
  }
  if (lag != 0) {
      base <- lag_hrf(base, lag = lag)
  }
  if (normalize) {
      base <- normalise_hrf(base)
  }

  attr(base, "name") <- key # Set name attribute to the matched registry key
  base
}

#' Evaluate an HRF Object
#'
#' This function evaluates a hemodynamic response function (HRF) object for a given set of time points (grid) and other parameters.
#' It handles both point evaluation (duration=0) and block evaluation (duration > 0).
#'
#' @param x The HRF object (inherits from `HRF` and `function`).
#' @param grid A numeric vector of time points at which to evaluate the HRF.
#' @param amplitude The scaling value for the event (default: 1).
#' @param duration The duration of the event (seconds). If > 0, the HRF is evaluated over this duration (default: 0).
#' @param precision The temporal resolution for evaluating responses when duration > 0 (default: 0.2).
#' @param summate Logical; whether the HRF response should accumulate over the duration (default: TRUE). If FALSE, the maximum response within the duration window is taken (currently only supported for single-basis HRFs).
#' @param normalize Logical; scale output so that the peak absolute value is 1 (default: FALSE). Applied *after* amplitude scaling and duration processing.
#' @param ... Additional arguments (unused).
#' @return A numeric vector or matrix of HRF values at the specified time points.
#' @export
evaluate.HRF <- function(x, grid, amplitude = 1, duration = 0,
                         precision = .2, summate = TRUE, normalize = FALSE, ...) {
  
  # Base function incorporating amplitude
  base <- function(g) amplitude * x(g)

  # Evaluate based on duration
  out <- if (duration < precision) {
      # Point evaluation
    base(grid)
  } else {
      # Block evaluation
    offs <- seq(0, duration, by = precision)
      # Evaluate HRF at shifted time points for each offset
      # Use lapply to handle potential matrix output from multi-basis HRFs
      hlist <- lapply(offs, function(o) base(grid - o))
      
      # Check if the result for the first offset is a matrix (multi-basis)
      is_multi_basis <- is.matrix(hlist[[1]])
      
      if (is_multi_basis) {
          # Combine matrices (summation is standard for multi-basis)
          if (summate) {
             Reduce("+", hlist)
    } else {
             # Taking max per-basis-column across offsets is non-standard and complex.
             # Sticking to summation for multi-basis block designs.
             warning("summate=FALSE is not typically used with multi-basis HRFs during block evaluation. Using summation.", call. = FALSE)
             Reduce("+", hlist)
          }
      } else {
          # Single basis HRF: hlist contains vectors, bind them into a matrix
        hmat <- do.call(cbind, hlist)
          if (summate) {
              rowSums(hmat)
          } else {
              # For single basis, take the max across the duration window at each grid point
              apply(hmat, 1, max, na.rm = TRUE) 
              # Alternative: find which offset gives max? apply(hmat, 1, function(vals) vals[which.max(vals)])
          }
      }
  }

  # Apply normalization if requested, handling matrix/vector case
  if (normalize) {
      if (is.matrix(out)) {
          apply(out, 2, function(col) {
              peak_val <- max(abs(col), na.rm = TRUE)
              if (!is.na(peak_val) && peak_val != 0) col / peak_val else col
          })
      } else {
          peak_val <- max(abs(out), na.rm = TRUE)
          if (!is.na(peak_val) && peak_val != 0) out / peak_val else out
      }
  } else {
    out
  }
}

# Create pre-defined HRF objects using generators -----

#' @rdname HRF_objects
#' @export
HRF_BSPLINE <- hrf_bspline_generator(nbasis=5, span=24)

#' @rdname HRF_objects  
#' @export
HRF_FIR <- hrf_fir_generator(nbasis=12, span=24)
