#' @keywords internal
#' @noRd
lss_fast <- function(dset, bdes, Y=NULL, use_cpp = TRUE) {
  # Data preparation
  if (is.null(Y)) {
    data_matrix <- get_data_matrix(dset)
  } else {
    data_matrix <- Y
  }

  # Check and validate Y dimensions
  if (!is.null(Y)) {
    expected_rows <- nrow(bdes$dmat_ran) # Use design matrix rows as reference
    if (nrow(Y) != expected_rows) {
      stop(sprintf("Y must have %d rows to match design matrix dimensions, but has %d rows", 
                  expected_rows, nrow(Y)))
    }
  }

  # Design matrices
  dmat_base <- as.matrix(bdes$dmat_base)
  dmat_fixed <- if (!is.null(bdes$fixed_ind)) as.matrix(bdes$dmat_fixed) else NULL
  dmat_ran <- as.matrix(bdes$dmat_ran)
  
  # Prepare nuisance matrix (baseline + fixed effects)
  if (!is.null(dmat_fixed)) {
    nuisance_matrix <- cbind(dmat_base, dmat_fixed)
  } else {
    nuisance_matrix <- dmat_base
  }
  
  # Use fmrilss package
  method <- if (use_cpp) "cpp_optimized" else "r_optimized"
  beta_matrix <- fmrilss::lss(Y = data_matrix, X = dmat_ran, Z = NULL, 
                              Nuisance = nuisance_matrix, method = method)
  
  return(beta_matrix)
}

#' @keywords internal
#' @noRd
lss_compute_r <- function(Q_dmat_ran, residual_data) {
  # This function is now deprecated - use fmrilss::lss instead
  stop("lss_compute_r is deprecated. Use fmrilss::lss instead.")
}