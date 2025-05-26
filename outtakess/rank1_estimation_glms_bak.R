r1_glms_betas <- function(X_list, X_all, y, Z = NULL, hrf_basis, hrf_ref, maxit = 100) {
  # Input validation
  n <- length(y)
  k <- length(X_list)           
  m <- ncol(hrf_basis)          
  q <- if (!is.null(Z)) ncol(Z) else 0
  
  # Construct X and X_all
  X <- do.call(cbind, X_list)  # n x (k*m)
  if (is.null(X_all)) {
    E <- kronecker(matrix(1, k, 1), diag(m))  # (k*m) x m
    X_all <- X %*% E  # n x m
  }
  
  # Debug dimensions
  cat("Dimensions check:\n")
  cat("X:", dim(X), "\n")
  cat("X_all:", dim(X_all), "\n")
  cat("y:", length(y), "\n")
  
  # Initial values
  init_h <- rep(1, m)            
  init_beta <- rep(1, k)         
  init_z <- rep(1, k)
  init_par <- c(init_h, init_beta, init_z)
  
  # Objective function
  fn <- function(w) {
    h <- w[1:m]
    beta <- w[(m + 1):(m + k)]
    z <- w[(m + k + 1):(m + k + k)]
    
    # Compute outer products
    hbeta <- tcrossprod(h, beta)  # m x k
    hz <- tcrossprod(h, z)        # m x k
    
    # Compute X_all contribution
    X_all_hz <- X_all %*% h  # n x 1
    
    # Compute residuals
    norm <- 0
    for (j in 1:k) {
      Xi <- X_list[[j]]  # n x m
      pred_i <- Xi %*% (hbeta[,j]) + X_all_hz * z[j]
      res_i <- y - pred_i
      norm <- norm + 0.5 * sum(res_i^2)
    }
    
    # Add regularization
    norm <- norm - 0.5 * sum(h^2)
    
    # Add HRF correlation constraint
    hrf_est <- hrf_basis %*% h
    corr <- sum(hrf_est * hrf_ref)
    if (corr < 0) norm <- norm + 1e10
    
    return(norm)
  }
  
  # Gradient function
  gradient <- function(w) {
    h <- w[1:m]
    beta <- w[(m + 1):(m + k)]
    z <- w[(m + k + 1):(m + k + k)]
    
    grad <- numeric(m + 2*k)
    hbeta <- tcrossprod(h, beta)  # m x k
    hz <- tcrossprod(h, z)        # m x k
    X_all_hz <- X_all %*% h  # n x 1
    
    # Initialize gradient components
    grad_h <- matrix(0, m, 1)
    grad_beta <- matrix(0, k, 1)
    grad_z <- matrix(0, k, 1)
    
    for (j in 1:k) {
      Xi <- X_list[[j]]  # n x m
      pred_i <- Xi %*% (hbeta[,j]) + X_all_hz * z[j]
      res_i <- y - pred_i
      
      # Gradient components
      Xi_res <- crossprod(Xi, res_i)  # m x 1
      X_all_res <- crossprod(X_all, res_i)  # m x 1
      
      grad_beta[j] <- sum(Xi_res * h)
      grad_z[j] <- sum(X_all_hz * res_i)
      grad_h <- grad_h + Xi_res * beta[j] + X_all_res * z[j]
    }
    
    # Add regularization gradient for h
    grad_h <- grad_h + h
    
    grad[1:m] <- grad_h
    grad[(m + 1):(m + k)] <- grad_beta
    grad[(m + k + 1):(m + 2*k)] <- grad_z
    
    return(-grad)
  }
  
  # Optimize with bounds
  opt <- stats::optim(init_par, fn, gr = gradient, method = "L-BFGS-B",
                      control = list(maxit = maxit))
  
  # Extract and process results
  h_opt <- opt$par[1:m]
  beta_opt <- opt$par[(m + 1):(m + k)]
  
  # Estimate HRF
  hrf_est <- hrf_basis %*% h_opt
  
  # Scale results
  scale_factor <- max(abs(hrf_est))
  if (scale_factor > 1e-10) {
    h_opt <- h_opt / scale_factor
    beta_opt <- beta_opt * scale_factor
    hrf_est <- hrf_est / scale_factor
  }
  
  # Ensure positive correlation
  corr <- sum(hrf_est * hrf_ref)
  if (corr < 0) {
    h_opt <- -h_opt
    beta_opt <- -beta_opt
    hrf_est <- -hrf_est
  }
  
  attr(beta_opt, "estimated_hrf") <- as.vector(hrf_est)
  attr(beta_opt, "basis_weights") <- h_opt
  
  return(beta_opt)
}

estimate_r1_glms <- function(dset, xdat, hrf_basis, hrf_ref, maxit = 100) {
  # Convert hrf_basis and hrf_ref to correct types
  hrf_basis <- as.matrix(hrf_basis)
  hrf_ref <- as.numeric(hrf_ref)
  
  # Extract data vectors
  vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) > 0))
  nvoxels <- length(vecs)
  n <- nrow(xdat$X)
  m <- ncol(hrf_basis)
  k <- ncol(xdat$X) / m
  
  # Create X_list
  X_list <- lapply(1:k, function(i) {
    Xi <- xdat$X[, ((i - 1) * m + 1):(i * m)]
    as.matrix(Xi)
  })
  
  estimated_hrfs <- matrix(NA, nrow = nrow(hrf_basis), ncol = nvoxels)
  
  # Estimate betas for each voxel
  res <- do.call(cbind, purrr::imap(vecs, function(v, i) {
    v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
    beta_v <- r1_glms_betas(
      X_list = X_list,
      X_all = NULL,  # Not used in this version
      y = as.numeric(v0),
      Z = NULL,
      hrf_basis = hrf_basis,
      hrf_ref = hrf_ref,
      maxit = maxit
    )
    estimated_hrfs[, i] <<- attr(beta_v, "estimated_hrf")
    beta_v
  }))
  
  list(beta_matrix = res, estimated_hrf = estimated_hrfs)
}

