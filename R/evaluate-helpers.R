#' @importFrom memoise memoise
#' @keywords internal
#' @noRd
.memo_hrf <- memoise::memoise(function(hrf, span, dt) {
    times <- seq(0, span, by = dt)
    # Evaluate HRF - ensure it returns a matrix
    val <- evaluate(hrf, times)
    if (is.vector(val)) matrix(val, ncol = 1) else val
})

#' Prepare Inputs for Regressor Evaluation Engines
#' 
#' Internal helper function to perform common setup steps before calling 
#' a specific evaluation engine (fft, conv, loop, Rconv).
#' Handles filtering of events, evaluation/memoization of HRF on fine grid.
#' 
#' @param x A `Reg` object.
#' @param grid The target evaluation time grid (numeric vector).
#' @param precision The precision for internal calculations (numeric scalar).
#' @return A list containing prepared inputs:
#'   * `nb`: Number of basis functions.
#'   * `hrf_span`: The span of the HRF.
#'   * `valid_ons`: Filtered onset times relevant to the grid.
#'   * `valid_durs`: Corresponding durations.
#'   * `valid_amp`: Corresponding amplitudes.
#'   * `grid`: The original target grid.
#'   * `precision`: The precision value.
#'   * `hrf_fine_matrix`: HRF values evaluated on the fine time grid (potentially memoized).
#'   * `fine_grid`: The fine time grid itself (if needed by Rconv/loop).
#'   * `summate`: Logical summation flag from the regressor.
#'   * `hrf`: The original HRF object.
#' @keywords internal
#' @noRd
#' @importFrom stats approx median convolve
prep_reg_inputs <- function(x, grid, precision) {
  
  # Ensure grid is sorted (Correctness 1.4)
  if (is.unsorted(grid)) {
      warning("Input grid is unsorted. Sorting grid for evaluation.")
      grid <- sort(grid)
  }
    
  nb <- nbasis(x$hrf) 
  hrf_span <- x$span 
  
  # Filter events based on grid boundaries and HRF span
  onset_min_bound <- grid[1] - hrf_span
  onset_max_bound <- grid[length(grid)]
  
  # Start with potentially already filtered data from Reg constructor
  keep_indices <- which(x$onsets >= onset_min_bound & x$onsets <= onset_max_bound)
  
  # Note: Amplitude filtering already done in Reg(), no need to repeat here
  valid_ons <- x$onsets[keep_indices]
  valid_durs <- x$duration[keep_indices]
  valid_amp <- x$amplitude[keep_indices]

  if (length(valid_ons) == 0) {
    # Return minimal info needed to signal zero output
    return(list(nb = nb, grid = grid, valid_ons = numeric(0)))
  }
  
  # Prepare/Memoize finely sampled HRF (Efficiency 2.1 / Ticket D-1)
  hrf_fine_matrix <- .memo_hrf(x$hrf, hrf_span, precision)
  
  # Prepare fine grid (needed for Rconv/loop interpolation)
  # Use full range of onsets when constructing the fine grid
  # to handle unsorted event inputs without reordering events
  fine_grid_start <- min(grid[1], min(valid_ons)) - hrf_span
  fine_grid_end <- max(grid[length(grid)], max(valid_ons) + max(valid_durs)) + hrf_span
  fine_grid <- seq(fine_grid_start, fine_grid_end, by = precision)

  return(list(
    nb         = nb,
    hrf_span   = hrf_span,
    valid_ons  = valid_ons,
    valid_durs = valid_durs,
    valid_amp  = valid_amp,
    grid       = grid,
    precision  = precision,
    hrf_fine_matrix = hrf_fine_matrix,
    fine_grid  = fine_grid, 
    summate    = x$summate,
    hrf        = x$hrf
  ))
}

# Internal Evaluation Engines -----

#' FFT-based Regressor Evaluation Engine
#' @param p A list returned by prep_reg_inputs.
#' @param ... Additional arguments.
#' @keywords internal
#' @noRd
eval_fft <- function(p, ...) {
  # Call the unified C++ wrapper
  result <- evaluate_regressor_cpp(
              grid = p$grid,
              onsets = p$valid_ons,
              durations = p$valid_durs,
              amplitudes = p$valid_amp,
              hrf_matrix = p$hrf_fine_matrix,
              hrf_span = p$hrf_span,
              precision = p$precision,
              method = "fft"
            )
  result
}

#' Direct Convolution Regressor Evaluation Engine
#' @param p A list returned by prep_reg_inputs.
#' @param ... Additional arguments.
#' @keywords internal
#' @noRd
eval_conv <- function(p, ...) {
  # Call the unified C++ wrapper
  result <- evaluate_regressor_cpp(
              grid = p$grid,
              onsets = p$valid_ons,
              durations = p$valid_durs,
              amplitudes = p$valid_amp,
              hrf_matrix = p$hrf_fine_matrix,
              hrf_span = p$hrf_span,
              precision = p$precision,
              method = "conv"
            )
  result
}

#' R Convolution Regressor Evaluation Engine
#' @param p A list returned by prep_reg_inputs.
#' @param ... Additional arguments.
#' @keywords internal
#' @noRd
#' @importFrom stats convolve approx
eval_Rconv <- function(p, ...) {
  # Check conditions (moved from evaluate.Reg)
  is_regular_grid <- length(p$grid) > 1 && length(unique(round(diff(p$grid), 8))) == 1
  is_constant_duration <- length(unique(p$valid_durs)) <= 1
  
  if (!is_regular_grid || !is_constant_duration) {
    warning("Method 'Rconv' requires a regular grid and constant event durations. Falling back to 'loop' method.")
    return(eval_loop(p, ...)) # Call the loop engine directly as fallback
  }
  
  # Proceed with R convolution using stats::convolve
  delta <- numeric(length(p$fine_grid))
  onset_indices <- round((p$valid_ons - p$fine_grid[1]) / p$precision) + 1
  valid_onset_indices <- onset_indices >= 1 & onset_indices <= length(p$fine_grid)
  delta[onset_indices[valid_onset_indices]] <- p$valid_amp[valid_onset_indices]
  
  samhrf <- p$hrf_fine_matrix # Already evaluated and potentially memoized
  nb <- p$nb
  
  if (nb > 1) {
    lowres <- matrix(0, length(p$grid), nb)
    for (b in 1:nb) {
      highres_conv <- stats::convolve(delta, rev(samhrf[, b]), type = "open")
      valid_len <- length(p$fine_grid)
      highres_trimmed <- highres_conv[1:valid_len]
      interp_res <- approx(p$fine_grid, highres_trimmed, xout = p$grid, rule = 2)$y
      lowres[, b] <- interp_res
    }
    result <- lowres
  } else {
    highres_conv <- stats::convolve(delta, rev(as.vector(samhrf)), type = "open")
    valid_len <- length(p$fine_grid)
    highres_trimmed <- highres_conv[1:valid_len]
    result <- approx(p$fine_grid, highres_trimmed, xout = p$grid, rule = 2)$y
  }
  result
}

#' R Loop Regressor Evaluation Engine
#' @param p A list returned by prep_reg_inputs.
#' @param ... Additional arguments passed to evaluate.HRF.
#' @keywords internal
#' @noRd
eval_loop <- function(p, ...) {
  # Add check for p$hrf
  if (is.null(p$hrf) || !inherits(p$hrf, "HRF")) {
      stop("Error inside eval_loop: p$hrf is NULL or not an HRF object.")
  }
  
  nb <- p$nb
  hrf_span <- p$hrf_span
  grid <- p$grid
  valid_ons <- p$valid_ons
  valid_durs <- p$valid_durs
  valid_amp <- p$valid_amp
  precision <- p$precision
  summate <- p$summate
  
  dspan <- hrf_span / stats::median(diff(grid), na.rm=TRUE) # Approx span in grid units
  
  # Pre-calculate nearest grid indices for onsets (more robust than RANN for this)
  # Find the index of the grid point *just before or at* each onset
  nidx <- findInterval(valid_ons, grid)
  nidx[nidx == 0] <- 1 
  
  outmat <- matrix(0, length(grid), length(valid_ons) * nb)

  for (i in seq_along(valid_ons)) { 
    start_grid_idx <- nidx[i]
    end_grid_idx <- min(start_grid_idx + ceiling(dspan) + 5, length(grid)) 
    if (start_grid_idx > length(grid)) next 
    grid.idx <- start_grid_idx:end_grid_idx
      
    relOns <- grid[grid.idx] - valid_ons[i]
    valid_rel_idx <- which(relOns >= 0 & relOns <= hrf_span)
      
    if (length(valid_rel_idx) > 0) {
        target_indices_outmat <- grid.idx[valid_rel_idx]
        # Call evaluate S3 generic, should dispatch to evaluate.HRF
        resp <- evaluate(p$hrf, relOns[valid_rel_idx], amplitude=valid_amp[i], 
                         duration=valid_durs[i], 
                         precision=precision,
                         summate=summate, ...)
                           
        if (!is.matrix(resp) && nb > 1) {
            resp <- matrix(resp, ncol=nb)
        }
        if (!is.matrix(resp) && nb == 1) {
            resp <- matrix(resp, ncol=1)
        }

        if (nrow(resp) != length(target_indices_outmat)){
            warning("Dimension mismatch between response and target indices in loop.")
            next
        }
                          
        if (nb > 1) {
            start_col <- (i-1) * nb + 1
            end_col <- i*nb 
            outmat[target_indices_outmat, start_col:end_col] <- resp
        } else {
            outmat[target_indices_outmat, i] <- resp
        }
    }
  }
  
  # Sum contributions across onsets
  if (length(valid_ons) > 1) {
    if (nb == 1) {
      result <- matrix(rowSums(outmat), ncol=1)
    } else {
      result <- do.call(cbind, lapply(1:nb, function(basis_idx) {
        rowSums(outmat[, seq(basis_idx, by=nb, length.out=length(valid_ons)), drop = FALSE])
      }))
    }
  } else { 
    if (nb == 1) {
        result <- matrix(outmat[,1], ncol=1) 
    } else {
        result <- outmat[, 1:nb, drop=FALSE] # Use drop=FALSE
    }
  }
  result
}
