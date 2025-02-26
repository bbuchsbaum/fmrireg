#' Rank-1 GLM Solver
#'
#' Jointly estimates activation coefficients (betas) and an HRF shape (in a given basis)
#' under a Rank-1 constraint. Allows optional nuisance regressors. 
#' Can optionally apply box constraints if \code{m=3} and sign-flip to ensure 
#' positive correlation with a reference HRF.
#'
#' @param X A numeric matrix of size \eqn{n \times (k*m)}, with n timepoints,
#'   k conditions, and m basis functions. Typically formed by convolving each
#'   event onset with \code{m} basis functions, then horizontally stacking them for k conditions.
#' @param y A length-n numeric vector of fMRI data (single voxel).
#' @param Z Optional numeric matrix of nuisance regressors, \eqn{n \times q}. If NULL, no nuisance.
#' @param hrf_basis A \eqn{T \times m} matrix of basis functions for the HRF.
#' @param hrf_ref A length-T numeric vector of some reference HRF shape, for an optional sign-flip check.
#' @param maxit Max number of L-BFGS-B iterations. Default 100.
#' @param flip_sign Logical. If \code{TRUE}, we ensure the final HRF is positively correlated with \code{hrf_ref}.
#' @param use_box_constraints Logical. If \code{TRUE} and \code{m=3}, we forcibly set \code{h[1]=1} and \code{h[2:3] \in [-1,1]}.
#'
#' @return A list with:
#' \describe{
#'   \item{beta}{Numeric vector of length k (event amplitudes).}
#'   \item{h}{Numeric vector of length m (basis weights).}
#'   \item{omega}{Numeric vector of length q for nuisance, or numeric(0) if none.}
#'   \item{converged}{Logical, TRUE if the L-BFGS-B optimizer converged.}
#'   \item{value}{The final objective (residual sum of squares / 2).}
#' }
#'
#' @examples
#' # Minimal usage example
#' set.seed(42)
#' n <- 200; k <- 10; m <- 3
#' X <- matrix(rnorm(n*k*m), n, k*m)
#' y <- rnorm(n)
#' hrf_basis <- matrix(rnorm(32*m), 32, m)
#' hrf_ref   <- dgamma(seq(0, 31, length.out=32), shape=6, rate=1)
#'
#' fit <- r1_glm_betas(X, y, NULL, hrf_basis, hrf_ref, maxit=50)
#' str(fit)
#' @export
r1_glm_betas <- function(X, y, Z = NULL, 
                         hrf_basis, 
                         hrf_ref,
                         maxit = 100,
                         flip_sign = FALSE,
                         use_box_constraints = FALSE) 
{
  n <- length(y)
  m <- ncol(hrf_basis)   # number of basis functions
  k <- ncol(X) / m       # number of conditions
  if (round(k) != k) {
    stop("X does not align with 'm'; X should have 'k*m' columns for integer k.")
  }
  q <- if (!is.null(Z)) ncol(Z) else 0
  
  # If Z is absent, treat it as empty matrix
  if (is.null(Z)) {
    Z <- matrix(0, n, 0)
    q <- 0
  }
  
  # Possibly do a thin-QR if n >> (k*m)
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
  
  # Initial guesses
  init_beta <- rep(0, k)      
  init_h    <- rep(0, m)         
  init_omega <- if (q > 0) rep(0, q) else numeric(0)
  z0 <- c(init_beta, init_h, init_omega)
  
  # Box constraints if m=3
  if (use_box_constraints && m == 3) {
    # Force h[1]=1, rest in [-1,1]
    lower <- c(rep(-Inf, k), 1, -1, -1)
    upper <- c(rep( Inf, k), 1,  1,  1)
  } else {
    lower <- rep(-Inf, length(z0))
    upper <- rep( Inf, length(z0))
  }
  
  if (q > 0) {
    # Extend for nuisance
    lower <- c(lower, rep(-Inf, q))
    upper <- c(upper, rep( Inf, q))
  }
  
  # The objective + gradient
  obj_grad <- function(z) {
    beta <- z[1:k]
    h    <- z[(k+1):(k+m)]
    omega <- if (q>0) z[(k+m+1):(k+m+q)] else numeric(0)
    
    # Pred = X * (beta ⊗ h) + Z * omega
    # We store them as "hb = outer(h, beta)" (dim [m x k]) => Flatten => length m*k
    # Then X %*% hb = predicted contribution from rank-1
    if (use_QR) {
      hb <- as.vector(outer(h, beta))
      pred <- X_tilde %*% hb
      if (q>0) pred <- pred + Z_tilde %*% omega
      res <- y_tilde - pred
    } else {
      hb <- as.vector(outer(h, beta))
      pred <- X %*% hb
      if (q>0) pred <- pred + Z %*% omega
      res <- y - pred
    }
    
    f <- 0.5*sum(res^2)
    # Slight penalty to deter h from being zero
    f <- f - 1e-6*h[1]^2
    
    # gradient wrt. (beta, h, omega)
    # (The direct summation approach is used here for clarity.)
    
    # grad wrt omega
    grad_omega <- if (q>0) -t(if(use_QR) Z_tilde else Z) %*% res else numeric(0)
    
    # We'll re-construct X if not using QR
    X_full <- if (use_QR) X_tilde else X
    
    grad_beta <- rep(0, k)
    grad_h    <- rep(0, m)
    
    # Partition X into k blocks of size m
    for (j in seq_len(k)) {
      idx_start <- (j-1)*m + 1
      idx_end   <- j*m
      Xj <- X_full[, idx_start:idx_end, drop=FALSE]
      # residual-based partial derivative
      # grad wrt beta[j] = -(Xj*h)^T * res
      Xj_h <- Xj %*% h
      grad_beta[j] <- - as.numeric(t(Xj_h) %*% res)
      
      # grad wrt h
      # For each h[i], partial derivative is sum(...) of -res * beta[j]*Xj
      # simpler as: -beta[j] * t(Xj) %*% res
      tmp <- t(Xj) %*% res
      grad_h <- grad_h - beta[j]*tmp
    }
    # small penalty on h[1]
    grad_h[1] <- grad_h[1] - 2e-6*h[1]
    
    grad <- c(grad_beta, grad_h, grad_omega)
    list(value=f, gradient=grad)
  }
  
  # Do the L-BFGS-B
  fit <- optim(
    par     = z0,
    fn      = function(z) obj_grad(z)$value,
    gr      = function(z) obj_grad(z)$gradient,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = list(maxit = maxit, factr = 1e7)
  )
  
  z_opt <- fit$par
  beta_opt <- z_opt[1:k]
  h_opt    <- z_opt[(k+1):(k+m)]
  omega_opt <- if (q>0) z_opt[(k+m+1):(k+m+q)] else numeric(0)
  
  # If not using box constraints, fix scale: max(|hrf|) = 1 => scale betas
  Bh <- as.vector(hrf_basis %*% h_opt)
  if (!(use_box_constraints && (m==3))) {
    sc <- max(abs(Bh))
    if (sc < 1e-12) sc <- 1
    h_opt    <- h_opt / sc
    beta_opt <- beta_opt * sc
    Bh       <- Bh / sc
  }
  
  # Optionally ensure HRF correlation is positive
  if (flip_sign) {
    if (sum(Bh*hrf_ref) < 0) {
      h_opt    <- -h_opt
      beta_opt <- -beta_opt
      Bh       <- -Bh
    }
  }
  
  list(
    beta      = beta_opt,
    h         = h_opt,
    omega     = omega_opt,
    converged = (fit$convergence == 0),
    value     = fit$value
  )
}



#' keywords internal
estimate_r1 <- function(dset, xdat, hrf_basis, hrf_ref, maxit = 100,
                       flip_sign = TRUE, use_box_constraints = FALSE) {
  # Ensure correct data types and validate inputs
  if (is.null(hrf_basis)) {
    stop("hrf_basis must not be NULL for the r1 method")
  }
  if (is.null(hrf_ref)) {
    stop("hrf_ref must not be NULL for the r1 method")
  }
  
  hrf_basis <- as.matrix(hrf_basis)
  hrf_ref   <- as.numeric(hrf_ref)
  
  # Pull out voxel data from the fmri_dataset
  vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) > 0))
  nvoxels <- length(vecs)
  
  if (nvoxels == 0) {
    stop("No voxels found in dataset mask")
  }
  
  # Show message about estimation (from original function)
  message("Estimating betas with rank-1 across ", nvoxels, " voxels ...")
  
  # Dimensions
  n <- nrow(xdat$X)        # number of time points
  m <- ncol(hrf_basis)     # number of basis functions
  
  if (m <= 0) {
    stop("hrf_basis must have at least one column")
  }
  
  # The stacked X has shape n x (k*m). So k is (ncol(xdat$X) / m).
  k <- ncol(xdat$X) / m
  if (round(k) != k) {
    stop("xdat$X does not match 'm'; ensure (k*m) columns for some integer k.")
  }
  
  if (k <= 0 || is.na(k)) {
    stop("Invalid k value: ", k, " (must be positive)")
  }
  
  # For the betas: each sub-design -> one amplitude => total k amplitudes
  # We'll store them in a k x nvox matrix
  beta_res <- matrix(NA, k, nvoxels)
  
  # Prepare to store the estimated HRF for each voxel
  estimated_hrfs <- matrix(NA, nrow = nrow(hrf_basis), ncol = nvoxels)
  
  for (i in seq_len(nvoxels)) {
    v <- vecs[[i]]
    # remove baseline
    v0 <- resid(lsfit(xdat$Base, v, intercept=FALSE))
    
    fit <- r1_glm_betas(X=as.matrix(xdat$X),
                        y=as.numeric(v0),
                        Z=NULL,
                        hrf_basis=hrf_basis,
                        hrf_ref=hrf_ref,
                        maxit=maxit,
                        flip_sign=flip_sign,
                        use_box_constraints=use_box_constraints)
    
    beta_res[,i] <- fit$beta
    # Reconstruct the actual HRF shape
    hshape <- as.vector(hrf_basis %*% fit$h)
    estimated_hrfs[,i] <- hshape
  }
  
  # Return the big matrix of betas and hrfs
  list(
    beta_matrix = beta_res,
    estimated_hrf = estimated_hrfs
  )
}

#' Estimate betas using Rank-1 GLM with separate designs (Mumford-style approach)
#'
#' This function estimates single-trial (or single-condition) betas in a rank-1
#' framework by splitting the design matrix \code{xdat$X} into a list of separate
#' sub-designs (\code{X_list}), one for each condition (or event). It then calls
#' \code{r1_glms_betas()} to jointly estimate the HRF shape and the per-condition
#' amplitudes. 
#'
#' @param dset An \code{fmri_dataset} object containing voxel data and mask
#' @param xdat A list typically returned by \code{get_X()}, containing:
#'    \itemize{
#'      \item \code{X}: The \emph{stacked} design matrix of dimension \eqn{n \times (k*m)}
#'      \item \code{Base}: The baseline design matrix, dimension \eqn{n \times p}
#'    }
#' @param hrf_basis A \eqn{T \times m} matrix of HRF basis functions
#' @param hrf_ref A length-\eqn{T} vector representing a reference HRF shape
#' @param maxit Maximum number of L-BFGS-B iterations (default 100)
#'
#' @return A list containing:
#'   \item{beta_matrix}{A \eqn{k \times nvox} matrix of estimated condition amplitudes}
#'   \item{estimated_hrf}{A \eqn{T \times nvox} matrix of the per-voxel estimated HRF shape}
#'
#' @details
#' This is the \emph{"Rank-1 GLM with Separate Designs"} approach (sometimes
#' called the Mumford method). We partition \code{xdat$X} into \eqn{k} sub-designs,
#' each \eqn{n \times m}. For each voxel, we remove the baseline using
#' \code{\link{lsfit}} and call \code{\link{r1_glms_betas}} to solve for the
#' event-wise betas under a rank-1 HRF constraint.
#'
#' @keywords internal
estimate_r1_glms <- function(dset, xdat, hrf_basis, hrf_ref, maxit = 100) {
  # Ensure correct data types
  hrf_basis <- as.matrix(hrf_basis)
  hrf_ref   <- as.numeric(hrf_ref)
  
  # Pull out voxel data from the fmri_dataset
  vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) > 0))
  nvoxels <- length(vecs)
  
  # Dimensions
  n <- nrow(xdat$X)        # number of time points
  m <- ncol(hrf_basis)     # number of basis functions
  # The stacked X has shape n x (k*m). So k is (ncol(xdat$X) / m).
  k <- ncol(xdat$X) / m
  if (round(k) != k) {
    stop("xdat$X does not match 'm'; ensure (k*m) columns for some integer k.")
  }
  
  # Partition xdat$X into a list of sub-designs: each sub-design is n x m
  X_list <- lapply(seq_len(k), function(i) {
    idx_start <- (i - 1)*m + 1
    idx_end   <- i*m
    Xi <- xdat$X[, idx_start:idx_end, drop=FALSE]
    as.matrix(Xi)
  })
  
  # Prepare to store the estimated HRF for each voxel
  estimated_hrfs <- matrix(NA, nrow = nrow(hrf_basis), ncol = nvoxels)
  
  # For the betas: each sub-design -> one amplitude => total k amplitudes
  # We'll store them in a k x nvox matrix
  beta_res <- matrix(NA, k, nvoxels)
  
  # Loop over voxels
  res <- purrr::imap_dfc(vecs, function(v, i) {
    # Remove baseline from the voxel time series
    v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
    
    # Call the "separate designs" rank-1 solver
    # r1_glms_betas is your function that handles sub-designs for each condition
    beta_v <- r1_glms_betas(
      X_list    = X_list,
      X_all     = NULL,    # Not used in this version
      y         = as.numeric(v0),
      Z         = NULL,    # Could pass if you want nuisance
      hrf_basis = hrf_basis,
      hrf_ref   = hrf_ref,
      maxit     = maxit
    )
    
    # The solver presumably attaches "estimated_hrf" as an attribute
    # or returns it as well. Let's assume it's in attr(..., "estimated_hrf")
    hrf_v <- attr(beta_v, "estimated_hrf")
    if (!is.null(hrf_v)) {
      estimated_hrfs[, i] <<- hrf_v
    }
    
    # We'll store the betas in the results object
    beta_v
  })
  
  # 'res' is a tibble with columns for each voxel. We want a numeric matrix k x nvox.
  # We can convert it easily:
  betas_mat <- as.matrix(res)
  # That is dimension [n, ???], so let's confirm that each column is length k:
  # Actually, each column might be length k, so betas_mat is k x nvox if the tibble is rowwise
  # or the opposite. Let's do t() if needed:
  
  # If 'res' was constructed via imap_dfc, each column is the entire vector "beta_v".
  # So each column is length k. That means betas_mat is k x nvox => perfect.
  
  list(
    beta_matrix   = betas_mat,
    estimated_hrf = estimated_hrfs
  )
}





