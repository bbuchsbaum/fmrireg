

#' @keywords internal
r1_glm_betas <- function(X, y, Z = NULL, hrf_basis, hrf_ref, maxit = 100) {
  # Dimensions
  n <- length(y)
  m <- ncol(hrf_basis)           # number of basis functions (3)
  k <- ncol(X) / m               # number of events (70)
  q <- if (!is.null(Z)) ncol(Z) else 0
  
  # Print dimensions for debugging
  message("Dimensions check:")
  message(sprintf("X: %d x %d", nrow(X), ncol(X)))
  message(sprintf("k (events): %d", k))
  message(sprintf("m (basis functions): %d", m))
  
  # Initial values
  init_beta <- rep(1, k)         # one per event (70)
  init_h <- rep(1, m)            # one per basis function (3)
  init_omega <- if (q > 0) rep(0, q) else numeric(0)
  z0 <- c(init_beta, init_h, init_omega)
  
  # Define bounds
  if (m == 3) {
    lower <- c(rep(-Inf, k), 1, -1, -1)
    upper <- c(rep(Inf, k), 1, 1, 1)
  } else {
    lower <- c(rep(-Inf, k), rep(-Inf, m))
    upper <- c(rep(Inf, k), rep(Inf, m))
  }
  if (q > 0) {
    lower <- c(lower, rep(-Inf, q))
    upper <- c(upper, rep(Inf, q))
  }
  
  # Objective function
  fn <- function(z) {
    beta <- z[1:k]               # event betas (70)
    h <- z[(k+1):(k+m)]         # basis weights (3)
    omega <- if (q > 0) z[(k+m+1):(k+m+q)] else numeric(0)
    
    # Reshape X to handle basis functions properly
    X_reshaped <- matrix(0, nrow=n, ncol=k)
    for(j in 1:m) {
      idx <- seq(j, ncol(X), by=m)
      X_reshaped <- X_reshaped + h[j] * X[,idx]
    }
    
    # Compute predicted response
    pred <- X_reshaped %*% beta
    if (q > 0) pred <- pred + Z %*% omega
    
    # Compute residual sum of squares
    res <- y - pred
    fval <- 0.5 * sum(res^2)
    return(fval)
  }
  
  # Gradient function
  gr <- function(z) {
    beta <- z[1:k]               # event betas (70)
    h <- z[(k+1):(k+m)]         # basis weights (3)
    omega <- if (q > 0) z[(k+m+1):(k+m+q)] else numeric(0)
    
    # Reshape X and compute residuals
    X_reshaped <- matrix(0, nrow=n, ncol=k)
    for(j in 1:m) {
      idx <- seq(j, ncol(X), by=m)
      X_reshaped <- X_reshaped + h[j] * X[,idx]
    }
    
    pred <- X_reshaped %*% beta
    if (q > 0) pred <- pred + Z %*% omega
    res <- y - pred
    
    # Gradient for betas
    grad_beta <- -t(X_reshaped) %*% res
    
    # Gradient for basis weights
    grad_h <- numeric(m)
    for(j in 1:m) {
      idx <- seq(j, ncol(X), by=m)
      grad_h[j] <- -sum(res * (X[,idx] %*% beta))
    }
    
    # Gradient for nuisance
    grad_omega <- if (q > 0) -t(Z) %*% res else numeric(0)
    
    return(c(grad_beta, grad_h, grad_omega))
  }
  
  # Run optimization
  res <- optim(z0, fn, gr, method = "L-BFGS-B", 
               lower = lower, upper = upper,
               control = list(maxit = maxit))
  
  # Extract results
  beta_opt <- res$par[1:k]
  h_opt <- res$par[(k+1):(k+m)]
  
  # Generate estimated HRF
  estimated_hrf <- hrf_basis %*% h_opt
  
  # Normalize and ensure positive correlation
  scale_factor <- max(abs(estimated_hrf))
  if (scale_factor < 1e-10) scale_factor <- 1
  h_opt <- h_opt / scale_factor
  beta_opt <- beta_opt * scale_factor
  
  # Ensure positive correlation with reference
  corr <- sum(estimated_hrf * hrf_ref)
  if (corr < 0) {
    h_opt <- -h_opt
    beta_opt <- -beta_opt
    estimated_hrf <- -estimated_hrf
  }
  
  # Attach attributes
  attr(beta_opt, "estimated_hrf") <- estimated_hrf
  attr(beta_opt, "basis_weights") <- h_opt
  
  return(beta_opt)
}
#' @keywords internal
IaXb <- function(X, a, b) {
  # Computes (I ⊗ a^T) %*% t(X) %*% b
  # X: [n x (d*k)]
  # a: [k x 1]
  # b: [n x 1]
  
  # Compute X^T * b
  Xt_b <- crossprod(X, b)  # [(d*k) x 1]

  # Reshape Xt_b to a matrix of shape [d, k]
  B <- matrix(Xt_b, nrow = d, ncol = length(a))
  
  # Multiply by a
  res <- B %*% a  # [d x 1]
  return(res)
}

#' @keywords internal
aIXb <- function(X, a, b) {
  # Computes (a^T ⊗ I) %*% t(X) %*% b
  # X: [n x (d*k)]
  # a: [d x 1]
  # b: [n x 1]
  
  # Compute X^T * b
  Xt_b <- crossprod(X, b)  # [(d*k) x 1]
  
  # Reshape Xt_b to a matrix of shape [k, d]
  B <- matrix(Xt_b, nrow = length(a), ncol = ncol(X) / length(a))
  B <- t(B)  # Now B is [k x d]
  
  # Multiply by a
  res <- B %*% a  # [k x 1]
  return(res)
}

#' @keywords internal
estimate_r1 <- function(dset, xdat, hrf_basis, hrf_ref, maxit = 100) {
  vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) > 0))
  xdat$X <- as.matrix(xdat$X)
  nvoxels <- length(vecs)
  
  estimated_hrfs <- matrix(NA, nrow = nrow(hrf_basis), ncol = nvoxels)
  message("Estimating betas using Rank-1 GLM method...")
  
  res <- do.call(cbind, purrr::imap(vecs, function(v, i) {
    v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
    beta_v <- r1_glm_betas(
      X = as.matrix(xdat$X),
      y = as.numeric(v0),
      Z = NULL,
      hrf_basis = as.matrix(hrf_basis),
      hrf_ref = as.numeric(hrf_ref),
      maxit = maxit
    )
    print(i)
    print(attr(beta_v, "estimated_hrf"))
    estimated_hrfs[, i] <<- as.vector(attr(beta_v, "estimated_hrf"))
    beta_v
  }))
  
  list(beta_matrix = as.matrix(res), estimated_hrf = estimated_hrfs)
}



#' Estimate activation coefficients and HRF using Rank-1 GLM
#'
#' This function implements the Rank-1 constrained GLM (R1-GLM) for joint 
#' estimation of activation coefficients (betas) and the HRF for a single voxel.
#' It uses a quasi-Newton approach (L-BFGS-B) for optimization. A Rank-1 
#' constraint enforces a shared HRF across all conditions, while allowing each 
#' condition to have its own amplitude.
#'
#' This implementation includes:
#'   - Optional QR decomposition for efficiency when n > p.
#'   - Optional box constraints for HRF parameters when m=3.
#'   - An optional sign-flip step to ensure positive correlation with a reference HRF.
#'   - A small penalty to prevent the HRF from collapsing to a near-zero vector.
#'
#' @param X Numeric matrix of size \eqn{n \times (k*m)}, where \eqn{n} is the number 
#'   of time points, \eqn{k} is the number of conditions, and \eqn{m} is the number 
#'   of basis functions. Typically constructed by convolving condition onsets 
#'   with HRF basis functions.
#' @param y Numeric vector of length \eqn{n} containing the BOLD signal at a single voxel.
#' @param Z Optional numeric matrix of nuisance regressors, of size \eqn{n \times q}. Default is NULL.
#' @param hrf_basis Numeric matrix of size \eqn{T \times m} with HRF basis functions 
#'   (one column per basis function).
#' @param hrf_ref Numeric vector of length \eqn{T} representing a reference HRF.
#' @param maxit Integer, maximum number of iterations for the optimizer. Default is 100.
#' @param flip_sign Logical, whether to enforce positive correlation with \code{hrf_ref}. 
#'   If TRUE and the estimated HRF is anti-correlated, it will be sign-flipped together 
#'   with the betas. If FALSE, negative deflections are allowed naturally. Default is FALSE.
#' @param use_box_constraints Logical, whether to use box constraints on HRF parameters 
#'   if \code{m=3}. If TRUE, HRF parameters are constrained to \cdoe{(-1,1)} with \code{h[1]=1}. 
#'   This strongly constrains the HRF shape and scale. Default is FALSE.
#'
#' @return A list with elements:
#'   \item{beta}{Numeric vector of length \eqn{k}, the estimated activation coefficients.}
#'   \item{h}{Numeric vector of length \eqn{m}, the estimated HRF coefficients in the chosen basis.}
#'   \item{omega}{Numeric vector of length \eqn{q} with nuisance parameters (if any).}
#'   \item{converged}{Logical indicating whether the optimizer converged.}
#'   \item{value}{Final value of the objective function.}
#'
#' @examples
#' # Mock data example:
#' set.seed(42)
#' n <- 200; k <- 10; m <- 3; q <- 5
#' X <- matrix(rnorm(n*k*m), n, k*m)
#' y <- rnorm(n)
#' Z <- matrix(rnorm(n*q), n, q)
#' hrf_basis <- matrix(rnorm(32*m), 32, m)
#' hrf_ref <- dgamma(seq(0, 31, length.out=32), shape=6, rate=1)
#'
#' res <- r1_glm_betas(X, y, Z, hrf_basis, hrf_ref, maxit=50, flip_sign=FALSE, use_box_constraints=FALSE)
#' str(res)
#'
#' @references
#' Pedregosa F, Eickenberg M, Ciuciu P, Thirion B, Gramfort A (2015).
#' Data-driven HRF estimation for encoding and decoding models. NeuroImage 104:209-220.
#'
#' @export
r1_glm_betas2 <- function(X, y, Z = NULL, hrf_basis, hrf_ref, maxit = 100,
                         flip_sign = FALSE, use_box_constraints = FALSE) {
  # Dimensions
  n <- length(y)
  m <- ncol(hrf_basis)       # number of basis functions
  k <- ncol(X) / m           # number of conditions
  q <- if (!is.null(Z)) ncol(Z) else 0
  
  if (round(k) != k) {
    stop("Number of conditions k is not an integer. Check dimensions of X and hrf_basis.")
  }
  
  # Handle nuisance regressors
  if (is.null(Z)) {
    Z <- matrix(0, n, 0)
    q <- 0
  }
  
  # Thin QR decomposition for efficiency if possible
  use_QR <- FALSE
  if (nrow(X) > ncol(X)) {
    qrX <- qr(X)
    if (qrX$rank == ncol(X)) {
      Q <- qr.Q(qrX)
      R <- qr.R(qrX)
      X_tilde <- R
      y_tilde <- t(Q) %*% y
      Z_tilde <- if (q > 0) t(Q) %*% Z else matrix(0, ncol(X), 0)
      use_QR <- TRUE
    } else {
      X_tilde <- X
      y_tilde <- y
      Z_tilde <- Z
    }
  } else {
    X_tilde <- X
    y_tilde <- y
    Z_tilde <- Z
  }
  
  # Initialization:
  # Start beta and h at zero for simplicity
  init_beta <- rep(0, k)      
  init_h <- rep(0, m)         
  init_omega <- if (q > 0) rep(0, q) else numeric(0)
  z0 <- c(init_beta, init_h, init_omega)
  
  # Set box constraints if requested and m=3
  if (use_box_constraints && m == 3) {
    # Force h[1]=1, and others in [-1,1]
    lower <- c(rep(-Inf, k), 1, -1, -1)
    upper <- c(rep( Inf, k), 1,  1,  1)
  } else {
    lower <- rep(-Inf, length(z0))
    upper <- rep( Inf, length(z0))
  }
  if (q > 0) {
    lower <- c(lower, rep(-Inf, q))
    upper <- c(upper, rep( Inf, q))
  }
  
  # Objective and Gradient Function
  # z = [beta(1:k), h(1:m), omega(1:q)]
  obj_grad <- function(z) {
    beta <- z[1:k]
    h <- z[(k+1):(k+m)]
    omega <- if (q > 0) z[(k+m+1):(k+m+q)] else numeric(0)
    
    # Compute predicted response and residual
    hb <- as.vector(outer(h, beta)) # Flattened by columns: vec(h beta^T)
    if (use_QR) {
      pred <- X_tilde %*% hb
      if (q > 0) pred <- pred + Z_tilde %*% omega
      res <- y_tilde - pred
    } else {
      pred <- X %*% hb
      if (q > 0) pred <- pred + Z %*% omega
      res <- y - pred
    }
    
    # Objective: 0.5 * ||res||^2 + small penalty
    f <- 0.5 * sum(res^2)
    # Prevent h from collapsing to zero vector:
    f <- f - 1e-6 * (h[1]^2)
    
    # Gradient computations
    X_full <- if (use_QR) X_tilde else X
    
    # grad w.r.t. omega: -Z^T res
    grad_omega <- if (q > 0) (-t(if (use_QR) Z_tilde else Z) %*% res) else numeric(0)
    
    # grad w.r.t. h and beta:
    # We have k blocks of size m in X_full.
    # grad_h = -Σ_j beta_j X_j^T res (summed over j), with penalty on h[1]
    # grad_beta_j = -(X_j h)^T res
    grad_h <- rep(0, m)
    grad_beta <- rep(0, k)
    
    for (j in seq_len(k)) {
      idx_start <- (j-1)*m + 1
      idx_end <- j*m
      Xj <- X_full[, idx_start:idx_end, drop=FALSE]
      XjTr <- t(Xj) %*% res           # m x 1
      grad_h <- grad_h - beta[j] * XjTr
      
      Xj_h <- Xj %*% h                # n x 1
      grad_beta[j] <- - as.numeric(t(Xj_h) %*% res)
    }
    
    # Add penalty gradient on h
    grad_h[1] <- grad_h[1] - 2e-6 * h[1]
    
    grad <- c(grad_beta, grad_h, grad_omega)
    list(value = f, gradient = grad)
  }
  
  # Optimize using L-BFGS-B
  res <- optim(
    par = z0,
    fn = function(z) obj_grad(z)$value,
    gr = function(z) obj_grad(z)$gradient,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = maxit, factr = 1e7)
  )
  
  z_opt <- res$par
  beta_opt <- z_opt[1:k]
  h_opt <- z_opt[(k+1):(k+m)]
  omega_opt <- if (q > 0) z_opt[(k+m+1):(k+m+q)] else numeric(0)
  
  # If not using box constraints with m=3, we should fix scale ambiguity.
  if (!(use_box_constraints && m == 3)) {
    Bh <- hrf_basis %*% h_opt
    scale_h <- max(abs(Bh))
    if (scale_h < 1e-12) scale_h <- 1.0
    h_opt <- h_opt / scale_h
    beta_opt <- beta_opt * scale_h
    Bh <- Bh / scale_h
  } else {
    # If using box constraints for m=3, scale is fixed by h[1]=1
    Bh <- hrf_basis %*% h_opt
  }
  
  # Optionally flip sign to match reference HRF direction
  if (flip_sign) {
    corr <- sum(Bh * hrf_ref)
    if (corr < 0) {
      h_opt <- -h_opt
      beta_opt <- -beta_opt
      Bh <- -Bh
    }
  }
  
  list(
    beta = beta_opt,
    h = h_opt,
    omega = omega_opt,
    converged = (res$convergence == 0),
    value = res$value
  )
}

#' Estimate activation coefficients and HRF using Rank-1 GLM with Separate Designs (R1-GLMS)
#'
#' This function implements the Rank-1 constrained GLM with Separate Designs (R1-GLMS) 
#' for joint estimation of activation coefficients (betas) and the HRF given a set of 
#' design matrices and BOLD response. It enforces a shared HRF across all conditions, 
#' yet allows the activation coefficients to differ. It uses a quasi-Newton approach 
#' (L-BFGS-B) for optimization.
#'
#' This version includes optional parameters to:
#'  - Apply a sign-flip step to enforce a positive correlation between the 
#'    estimated HRF and a given reference HRF.
#'  - Apply simple box constraints on HRF parameters when using exactly three 
#'    basis functions, mirroring the approach used in the provided code snippet.
#'
#' @param X_list A list of numeric matrices, each of size \eqn{n \times (k*m)}, where \eqn{n} is the number 
#'   of time-points, \eqn{k} is the number of conditions, and \eqn{m} is the number 
#'   of basis functions.  Each matrix in the list should be a separate design matrix 
#'   for a given condition. Usually constructed by convolving condition onsets with basis functions.
#' @param y Numeric vector of length \eqn{n} containing the BOLD signal at a single voxel.
#' @param Z Optional numeric matrix of nuisance regressors of size \eqn{n \times q}.
#' @param hrf_basis Numeric matrix \eqn{T \times m} with HRF basis functions 
#'   (one column per basis function).
#' @param hrf_ref Numeric vector of length \eqn{T} representing a reference HRF.
#' @param maxit Integer, maximum number of iterations for the optimizer. Default is 100.
#' @param flip_sign Logical, whether to enforce a positive correlation between the 
#'   estimated HRF and the reference HRF. Default is TRUE.
#' @param use_box_constraints Logical, whether to use box constraints on the HRF parameters 
#'   if \code{m = 3}. The constraints used are similar to those provided in the other snippet. 
#'   By default FALSE. If TRUE and \code{m=3}, HRF parameters are forced into \code{(-1,1)} and 
#'   scaled so that the first parameter is constrained to 1. This is a strong constraint.
#'
#' @return A list with elements:
#'   \item{beta}{Numeric vector of length \eqn{k} representing the estimated activation coefficients.}
#'   \item{h}{Numeric vector of length \eqn{m} representing the estimated HRF coefficients in the chosen basis.}
#'   \item{omega}{Numeric vector of length \eqn{q} with nuisance parameters (if any).}
#'   \item{r}{Numeric vector of length \eqn{k} with nuisance parameters used for the model with separate designs (if any)}
#'   \item{converged}{Logical indicating whether the optimizer converged.}
#'   \item{value}{Final value of the objective function.}
#'
#' @examples
#' # Example (mock data):
#' set.seed(42)
#' n <- 200; k <- 10; m <- 3; q <- 5
#' X_list <- lapply(1:k, function(i) matrix(rnorm(n*m*2), n, m*2))
#' y <- rnorm(n)
#' Z <- matrix(rnorm(n*q), n, q)
#' hrf_basis <- matrix(rnorm(32*m), 32, m) # 32 time points HRF model
#' hrf_ref <- dgamma(seq(0, 31, length.out=32), shape=6, rate=1) # a gamma HRF shape
#'
#' res <- r1_glm_separate_betas(
#'   X_list, y, Z, hrf_basis, hrf_ref, maxit=50, flip_sign=FALSE, use_box_constraints=TRUE
#' )
#' str(res)
#'
#' @references
#' Pedregosa F, Eickenberg M, Ciuciu P, Thirion B, Gramfort A (2015) 
#' Data-driven HRF estimation for encoding and decoding models. NeuroImage 104:209-220.
#'
#' @export
r1_glm_separate_betas2 <- function(X_list, y, Z = NULL, hrf_basis, hrf_ref, maxit = 100,
                                  flip_sign = TRUE, use_box_constraints = FALSE) {
  # Dimensions
  n <- length(y)                  # Number of timepoints
  m <- ncol(hrf_basis)          # number of basis functions
  k <- length(X_list)             # number of conditions
  q <- if (!is.null(Z)) ncol(Z) else 0 # number of nuisance regressors
  
  # Handle Z if NULL
  if (is.null(Z)) {
    Z <- matrix(0, n, 0)
    q <- 0
  }
  
  # Check dimensions
  for (X in X_list){
    if (nrow(X) != n) {
      stop("Number of rows in design matrix and y must match.")
    }
  }
  
  # Thin QR decomposition for efficiency if beneficial
  # For simplicity, only perform QR decomposition if all matrices are full rank
  use_QR <- FALSE
  qr_list <- lapply(X_list, function(X) if (nrow(X) > ncol(X)) qr(X) else NULL)
  if (all(sapply(qr_list, function(qrX) !is.null(qrX) && qrX$rank == ncol(qrX)))) {
    Q_list <- lapply(qr_list, qr.Q)
    R_list <- lapply(qr_list, qr.R)
    X_tilde_list <- R_list
    y_tilde <- t(Reduce('+', Q_list) / k) %*% y # Average QR over designs
    if (q > 0) {
      Z_tilde <- t(Reduce('+', Q_list) / k) %*% Z
    } else {
      Z_tilde <- matrix(0, ncol(X_list[[1]]), 0)
    }
    use_QR <- TRUE
  } else {
    X_tilde_list <- X_list
    y_tilde <- y
    Z_tilde <- Z
  }
  
  # Initial values
  init_beta <- rep(0, k)          # Initial beta
  init_h <- rep(1 / m, m)          # Initial h
  init_r <- rep(0, k)            # Initial r (for separate designs)
  init_omega <- if (q > 0) rep(0, q) else numeric(0)  # Initial omega (nuisance regressors)
  z0 <- c(init_h, init_beta, init_r, init_omega)     # Combine initial values into a single vector for optimization
  
  
  # Set box constraints if requested
  # If use_box_constraints=TRUE and m=3, we replicate the idea of forcing h in [-1,1],
  # and also fix h[1]=1 (enforcing scale?). This is a strong constraint.
  # Otherwise, no constraints.
  if (use_box_constraints && m == 3) {
    lower <- c(1, -1, -1, rep(-Inf, k), rep(-Inf, k))
    upper <- c(1,  1,  1, rep( Inf, k), rep( Inf, k))
  } else {
    lower <- rep(-Inf, length(z0) )
    upper <- rep( Inf, length(z0) )
  }
  if (q > 0) {
    lower <- c(lower, rep(-Inf, q))
    upper <- c(upper, rep( Inf, q))
  }
  
  # Objective and gradient function
  obj_grad <- function(z) {
    # Extract parameters from optimization vector
    h <- z[1:m]
    beta <- z[(m + 1):(m + k)]
    r <- z[(m + k + 1):(m + 2*k)] # nuisance regressor for separate designs
    omega <- if (q > 0) z[(m + 2*k + 1):(m + 2*k + q)] else numeric(0)
    
    f <- 0 # Initialize objective function value
    
    grad_h <- rep(0, m)        # Initialize gradient for h
    grad_beta <- rep(0, k)    # Initialize gradient for beta
    grad_r <- rep(0, k)        # Initialize gradient for r
    
    # Loop over each condition design matrix
    for (j in seq_len(k)) {
      
      if (use_QR) {
        Xj <- X_tilde_list[[j]]
        # Construct predicted response for this condition
        pred <-  (Xj %*% h)*beta[j] + (Xj %*% h) * r[j]
        if (q>0) pred <- pred + Z_tilde %*% omega
        res <- y_tilde - pred
        # Gradient Computation
        XjTr <- t(Xj) %*% res
        grad_h <- grad_h - (XjTr * (beta[j] + r[j]))
        grad_beta[j] <- - as.numeric(t(Xj %*% h) %*% res)
        grad_r[j] <- - as.numeric(t(Xj %*% h) %*% res)
        
      } else {
        Xj <- X_list[[j]]
        # Construct predicted response for this condition
        pred <-  (Xj %*% h)*beta[j] + (Xj %*% h) * r[j]
        if (q>0) pred <- pred + Z %*% omega
        res <- y - pred
        # Gradient Computation
        XjTr <- t(Xj) %*% res
        grad_h <- grad_h - (XjTr * (beta[j] + r[j]))
        grad_beta[j] <- - as.numeric(t(Xj %*% h) %*% res)
        grad_r[j] <- - as.numeric(t(Xj %*% h) %*% res)
      }
      
      f <- f + 0.5 * sum(res^2)
    }
    
    # small negative penalty to avoid h collapsing
    f <- f - 1e-6 * (h[1]^2)
    
    # Gradient wrt omega
    grad_omega <- if (q > 0) -t(if (use_QR) Z_tilde else Z) %*% res else numeric(0)
    
    # Add penalty grad wrt h
    grad_h[1] <- grad_h[1] - 2e-6 * h[1]
    
    grad <- c(grad_h, grad_beta, grad_r, grad_omega)   # Combine gradients into a single vector
    list(value = f, gradient = grad)                 # Return objective function value and gradient
  }
  
  # Run optimization
  res <- optim(
    par = z0,
    fn = function(z) { obj_grad(z)$value },
    gr = function(z) { obj_grad(z)$gradient },
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = maxit, factr = 1e7)
  )
  
  z_opt <- res$par                                       # Extract optimized parameters
  h_opt <- z_opt[1:m]                                  # HRF coefficients
  beta_opt <- z_opt[(m + 1):(m + k)]                    # Activation coefficients
  r_opt <- z_opt[(m + k + 1):(m + 2 * k)]                  # Nuisance parameter for separate design
  omega_opt <- if (q > 0) z_opt[(m + 2 * k + 1):(m + 2 * k + q)] else numeric(0) # Nuisance regressors
  
  # Normalization step (only if not using strict box constraints that fix h):
  # If box constraints are not used or do not fix scale, we can still normalize
  # so that ||B h||∞ = 1. This fixes scale ambiguity.
  # Note: If we used the constraints with m=3, h[1]=1, this may already fix scale.
  if (!(use_box_constraints && m == 3)) {
    Bh <- hrf_basis %*% h_opt
    scale_h <- max(abs(Bh))
    if (scale_h < 1e-12) scale_h <- 1.0
    h_opt <- h_opt / scale_h
    beta_opt <- beta_opt * scale_h
    r_opt <- r_opt * scale_h
    Bh <- Bh / scale_h
  } else {
    # If we used constraints forcing h into [-1,1], scale is somewhat already fixed.
    # Do nothing here or could do additional checks if desired.
    Bh <- hrf_basis %*% h_opt
  }
  
  # If flip_sign is TRUE, ensure positive correlation with reference HRF:
  # If user wants to allow negative deflections, set flip_sign=FALSE.
  if (flip_sign) {
    corr <- sum(Bh * hrf_ref)
    if (corr < 0) {
      h_opt <- -h_opt
      beta_opt <- -beta_opt
      r_opt <- -r_opt
      Bh <- -Bh
    }
  }
  
  list(
    h = h_opt,             # Estimated HRF coefficients
    beta = beta_opt,       # Estimated activation coefficients
    r = r_opt,
    omega = omega_opt,       # Nuisance regressor weights
    converged = (res$convergence == 0), # Optimization convergence status
    value = res$value         # Final objective function value
  )
}







