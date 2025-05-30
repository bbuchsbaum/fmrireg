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

  n_timepoints <- nrow(data_matrix)
  n_voxels <- ncol(data_matrix)
  
  # Design matrices
  dmat_base <- as.matrix(bdes$dmat_base)
  dmat_fixed <- if (!is.null(bdes$fixed_ind)) as.matrix(bdes$dmat_fixed) else NULL
  dmat_ran <- as.matrix(bdes$dmat_ran)
  
  # Prepare baseline and fixed design matrix
  if (!is.null(dmat_fixed)) {
    X_base_fixed <- cbind(dmat_base, dmat_fixed)
  } else {
    X_base_fixed <- dmat_base
  }
  
  if (use_cpp) {
    ##print("using cpp")
    res <- compute_residuals_cpp(X_base_fixed, data_matrix, dmat_ran)
    beta_matrix <- lss_compute_cpp(res$Q_dmat_ran, res$residual_data)
  } else {
    ##print("using r")
    # Pure R implementation
    P_base_fixed <- MASS::ginv(X_base_fixed)
    Q_base_fixed <- diag(n_timepoints) - X_base_fixed %*% P_base_fixed
    residual_data <- Q_base_fixed %*% data_matrix
    Q_dmat_ran <- Q_base_fixed %*% dmat_ran
    beta_matrix <- lss_compute_r(Q_dmat_ran, residual_data)
  }
  
  return(beta_matrix)
}

#' @keywords internal
#' @noRd
lss_compute_r <- function(Q_dmat_ran, residual_data) {
  n_timepoints <- nrow(Q_dmat_ran)
  n_events <- ncol(Q_dmat_ran)
  n_voxels <- ncol(residual_data)
  
  # Initialize beta matrix
  beta_matrix <- matrix(NA, nrow = n_events, ncol = n_voxels)
  
  # Loop over trials
  for (i in seq_len(n_events)) {
    # Current trial regressor
    c <- Q_dmat_ran[, i, drop = FALSE]
    
    # Skip if current regressor has very low variance
    c_norm <- drop(crossprod(c))
    if (c_norm < 1e-10) {
      beta_matrix[i, ] <- 0
      next
    }
    
    # If there's only one event, we can just do a simple regression
    if (n_events == 1) {
      s <- c / c_norm
      beta_matrix[i, ] <- drop(crossprod(s, residual_data))
      next
    }
    
    # Sum of other trial regressors
    b <- Q_dmat_ran[, -i, drop = FALSE]
    
    # If all other regressors are very close to zero, just use this regressor
    if (ncol(b) == 0 || all(colSums(b^2) < 1e-10)) {
      s <- c / c_norm
      beta_matrix[i, ] <- drop(crossprod(s, residual_data))
      next
    }
    
    # Sum the columns of b
    b_sum <- rowSums(b)
    b <- matrix(b_sum, ncol=1)
    
    # Check if b is essentially zero
    b_norm <- drop(crossprod(b))
    if (b_norm < 1e-10) {
      # If other regressors sum to approximately zero, use this regressor directly
      s <- c / c_norm
      beta_matrix[i, ] <- drop(crossprod(s, residual_data))
      next
    }
    
    # Compute v = c - b(b'b)^(-1)b'c
    bc <- drop(crossprod(b, c))
    v <- c - (bc/b_norm) * b
    
    # Compute cvdot = c'v
    cvdot <- drop(crossprod(c, v))
    
    # Handle numerical stability - only set to zero if truly negligible
    # compared to both c_norm and v_norm
    v_norm <- drop(crossprod(v))
    if (abs(cvdot) < 1e-8 * sqrt(c_norm * v_norm)) {
      # Try direct regression as fallback
      s <- c / c_norm
    } else {
      s <- v / cvdot
    }
    
    # Compute beta for current trial
    beta_matrix[i, ] <- drop(crossprod(s, residual_data))
  }
  
  return(beta_matrix)
}