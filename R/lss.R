#' @keywords internal
#' @noRd
lss_fast <- function(dset, bdes, use_cpp = TRUE) {
  # Data preparation
  data_matrix <- get_data_matrix(dset)
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
    print("using cpp")
    res <- compute_residuals_cpp(X_base_fixed, data_matrix, dmat_ran)
    beta_matrix <- lss_compute_cpp(res$Q_dmat_ran, res$residual_data)
  } else {
    print("using r")
    # Pure R implementation
    P_base_fixed <- MASS::ginv(X_base_fixed)
    Q_base_fixed <- diag(n_timepoints) - X_base_fixed %*% P_base_fixed
    residual_data <- Q_base_fixed %*% data_matrix
    Q_dmat_ran <- Q_base_fixed %*% dmat_ran
    beta_matrix <- lss_compute_r(Q_dmat_ran, residual_data)
  }
  
  return(list(beta_matrix = beta_matrix, estimated_hrf = NULL))
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
    
    # Sum of other trial regressors
    b <- rowSums(Q_dmat_ran[, -i, drop = FALSE])
    b <- matrix(b, ncol=1)
    
    # Compute v = c - b(b'b)^(-1)b'c
    b_norm <- drop(crossprod(b))
    if(b_norm > 1e-10) {
      bc <- drop(crossprod(b, c))
      v <- c - (bc/b_norm) * b
    } else {
      v <- c
    }
    
    # Compute cvdot = c'v
    cvdot <- drop(crossprod(c, v))
    
    # Handle numerical stability
    if (abs(cvdot) < 1e-5 * drop(crossprod(c))) {
      s <- rep(0, n_timepoints)
    } else {
      s <- v / cvdot
    }
    
    # Compute beta for current trial
    beta_matrix[i, ] <- drop(crossprod(s, residual_data))
  }
  
  return(beta_matrix)
}