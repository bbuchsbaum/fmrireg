#' @noRd
#' @keywords internal
ridge_betas <- function(X, Y, penalty_factor=rep(1:ncol(X)), lambda=.01) {
  with_package("glmnet")
  fit <- glmnet::glmnet(X, Y, penalty.factor=penalty_factor, alpha=0,lambda=lambda)
  coef(fit)[,1,drop=FALSE]
}


#' @noRd
#' @keywords internal
pls_betas <- function(X, Y, ncomp=3) {
  with_package("pls")
  dx <- list(X=as.matrix(X), Y=Y)
  fit <- pls::plsr(Y ~ X, data=dx, ncomp=ncomp, method="simpls", scale=TRUE, maxit=500)
  coef(fit, ncomp=ncomp)[,,1]
}


#' @keywords internal
#' @noRd
pls_global_betas <- function(X, Y, ncomp=3) {
  with_package("pls")
  dx <- list(X=as.matrix(X), Y=Y)
  fit <- pls::plsr(Y ~ X, data=dx, ncomp=ncomp, method="widekernelpls", scale=TRUE, maxit=500)
  coef(fit, ncomp=ncomp)[,,1]
}


#' @keywords internal
#' @noRd
ols_betas <- function(X, Y) {
  fit <- lm.fit(as.matrix(X),Y)
  coef(fit)
}


#' @keywords internal
#' @noRd
slm_betas <- function(X, Y) {
  with_package("care")
  slm.1 <- care::slm(X, Y, verbose=FALSE)
  b2 <- coef(slm.1)[,-(1:2)]
  b1 <- coef(slm.1)[,1]
  b1 + b2
}


#' @noRd
#' @keywords internal
mixed_betas <- function(X, Y, ran_ind, fixed_ind) {
  with_package("rrBLUP")  
  fit <- rrBLUP::mixed.solve(Y, Z=X[,ran_ind], X=X[,c(fixed_ind)], bounds=c(c(1e-05, .2)))
  c(fit$u, fit$b)
}

mixed_betas_cpp <- function(X, Y, ran_ind, fixed_ind) {
  fit <- mixed_solve_internal(as.matrix(Y), Z=X[,ran_ind,drop=FALSE], X=X[,c(fixed_ind),drop=FALSE],  bounds=c(c(1e-05, .2)))
  c(fit$u, fit$b)
}



#' Estimate betas using the Rank-1 GLM methods (R1 and R1-GLMS with optimizations)
#'
#' This function estimates betas (regression coefficients) and the hemodynamic response function (HRF)
#' simultaneously using the Rank-1 GLM (`r1`) and Rank-1 GLM with Mumford's separate beta estimation (`r1_glms`) methods.
#' It includes optimizations to improve computational efficiency.
#'
#' @param x An object of class `fmri_dataset` representing the fMRI dataset.
#' @param fixed A formula specifying the fixed regressors that model constant effects (i.e., non-varying over trials).
#' @param ran A formula specifying the random (trialwise) regressors that model single trial effects.
#' @param block A formula specifying the block factor.
#' @param method The regression method for estimating trialwise betas; use `"r1"` for the Rank-1 GLM method or `"r1_glms"` for the Rank-1 GLM with Mumford's approach.
#' @param basemod A `baseline_model` instance to regress out of data before beta estimation (default: NULL).
#' @param hrf_basis A matrix of basis functions for the HRF (default: NULL).
#' @param hrf_ref A reference HRF vector for initializing and constraining the HRF estimation (default: NULL).
#' @param maxit Maximum number of iterations for the optimization (default: 100).
#' @param ... Additional arguments passed to the estimation method.
#'
#' @return A list of class "fmri_betas" containing the following components:
#'   * betas_fixed: NeuroVec object representing the fixed effect betas.
#'   * betas_ran: NeuroVec object representing the random effect betas.
#'   * design_ran: Design matrix for random effects.
#'   * design_fixed: Design matrix for fixed effects.
#'   * design_base: Design matrix for baseline model.
#'   * basemod: Baseline model object.
#'   * fixed_model: Fixed effect model object.
#'   * ran_model: Random effect model object.
#'   * estimated_hrf: The estimated HRF vector.
#'
#' @details
#' The `r1` method uses the Rank-1 GLM approach to jointly estimate the HRF and activation coefficients.
#' The `r1_glms` method implements the Mumford approach by estimating each beta individually to reduce correlations,
#' treating all events as coming from one condition.
#' This implementation includes optimizations to improve computational efficiency:
#' - Precomputing the total sum of all trial regressors to avoid redundant computations.
#' - Precomputing the QR decomposition of design matrices to speed up linear algebra operations.
#'
#' @references
#' Pedregosa, F., et al. (2015). GLM with Rank-1 constraint (R1-GLM): a fast, spatially adaptive model for single trial fMRI data.
#' \emph{NeuroImage}, 104, 271–285.
#'
#' Mumford, J. A., Turner, B. O., Ashby, F. G., & Poldrack, R. A. (2012). Deconvolving BOLD activation in event-related designs for multivoxel pattern classification analyses. \emph{NeuroImage}, 59(3), 2636–2643.
#'
#' @seealso \code{\link{fmri_dataset}}, \code{\link{baseline_model}}, \code{\link{event_model}}
#'
#' @examples
#' \dontrun{
#' facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
#' facedes$frun <- factor(facedes$run)
#' scans <- paste0("rscan0", 1:6, ".nii")
#'
#' dset <- fmri_dataset(scans=scans, mask="mask.nii", TR=1.5, run_length=rep(436,6), event_table=facedes)
#' fixed = onset ~ hrf(run)
#' ran = onset ~ trialwise()
#' block = ~ run
#'
#' betas <- estimate_betas(dset, fixed=fixed, ran=ran, block=block, method="r1_glms")
#' }
#' @export
estimate_betas.fmri_dataset <- function(x, fixed = NULL, ran, block,
                                        method = c("mixed", "mixed_cpp", "lss", "lss_naive", "lss_cpp",
                                                   "r1", "pls", "pls_searchlight", 
                                                   "pls_global", "ols", "lowrank_hrf"),
                                        basemod = NULL,
                                        hrf_basis = NULL,
                                        hrf_ref = NULL,
                                        maxit = 100,
                                        ...) {
  method <- match.arg(method)
  dset <- x
  bvec <- get_data(dset)
  mask <- get_mask(dset)
  
  if (method == "r1_glms" || method == "r1") {
    if (is.null(hrf_basis)) {
      stop("Please provide 'hrf_basis', a matrix of HRF basis functions.")
    }
    if (is.null(hrf_ref)) {
      stop("Please provide 'hrf_ref', a reference HRF vector.")
    }
  }
  
  bmod <- if (is.null(basemod)) {
    baseline_model("constant", sframe = dset$sampling_frame)
  } else {
    basemod
  }
  
  bdes <- gen_beta_design(fixed, ran, block, bmod, dset, method = method)
  betas <- run_estimate_betas(bdes, dset, method, hrf_basis = hrf_basis, hrf_ref = hrf_ref, maxit = maxit, ...)
  
  # Check dimensions before indexing
  message(sprintf("beta_matrix dimensions: %d x %d", 
                  nrow(betas$beta_matrix), ncol(betas$beta_matrix)))
  message(sprintf("ran_ind length: %d", length(bdes$ran_ind)))
  
  
  nbetas <- nrow(betas$beta_matrix)
  ospace_ran <- neuroim2::add_dim(neuroim2::space(mask), length(bdes$ran_ind))
  
  if (!is.null(bdes$fixed_ind)) {
    ospace_fixed <- neuroim2::add_dim(neuroim2::space(mask), length(bdes$fixed_ind))
    fixed <- neuroim2::NeuroVec(as.matrix(betas$beta_matrix[bdes$fixed_ind, , drop = FALSE]), ospace_fixed, mask = mask)
  } else {
    fixed <- NULL
  }
  
  ran <- neuroim2::NeuroVec(as.matrix(betas$beta_matrix[bdes$ran_ind, , drop = FALSE]), ospace_ran, mask = mask)
  
  ret <- list(betas_fixed = fixed,
              betas_ran = ran,
              design_ran = bdes$dmat_ran,
              design_fixed = bdes$dmat_fixed,
              design_base = bdes$dmat_base,
              basemod = basemod,
              fixed_model = bdes$emod_fixed,
              ran_model = bdes$emod_ran,
              estimated_hrf = betas$estimated_hrf)
  
  class(ret) <- "fmri_betas"
  ret
}


#' @keywords internal
run_estimate_betas <- function(bdes, dset, method, hrf_basis = NULL, hrf_ref = NULL, maxit = 100, ncomp=4, ...) {
  get_X <- function() {
    X <- if (is.null(bdes$fixed)) bdes$dmat_ran else cbind(bdes$dmat_ran, bdes$dmat_fixed)
    Base <- as.matrix(bdes$dmat_base)
    X[is.na(X)] <- 0
    list(Base = Base, X = X)
  }
  
  if (method == "mixed") {
    vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) > 0))
    xdat <- get_X()
    res <- do.call(cbind, furrr::future_map(vecs, function(v) {
      v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
      mixed_betas(xdat$X, v0, ran_ind = seq_len(ncol(bdes$dmat_ran)), 
                 fixed_ind = if(!is.null(bdes$dmat_fixed)) {
                   (ncol(bdes$dmat_ran) + 1):(ncol(bdes$dmat_ran) + ncol(bdes$dmat_fixed))
                 } else {
                   NULL
                 })
    }))
    list(beta_matrix=as.matrix(res), estimated_hrf=NULL)
  } else if (method == "mixed_cpp") {
    vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) > 0))
    xdat <- get_X()
    res <- do.call(cbind, furrr::future_map(vecs, function(v) {
      v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
      mixed_betas_cpp(as.matrix(xdat$X), v0, ran_ind = seq_len(ncol(bdes$dmat_ran)), 
                  fixed_ind = if(!is.null(bdes$dmat_fixed)) {
                    (ncol(bdes$dmat_ran) + 1):(ncol(bdes$dmat_ran) + ncol(bdes$dmat_fixed))
                  } else {
                    NULL
                  })
    }))
    
    list(beta_matrix=as.matrix(res), estimated_hrf=NULL)
    
  
  }  else if (method == "lss_naive") {
    lss_naive(dset, bdes)
  } else if (method == "r1_glms") {
    xdat <- get_X()
    estimate_r1_glms(dset, xdat, hrf_basis, hrf_ref, maxit)
  } else if (method == "r1") {
    xdat <- get_X()
    estimate_r1(dset, xdat, hrf_basis, hrf_ref, maxit)
  } else if (method == "lss") {
    lss_fast(dset, bdes, use_cpp = FALSE)
  } else if (method == "lss_cpp") {
    lss_fast(dset, bdes, use_cpp = TRUE)
  } else if (method == "pls") {
    vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) >0))
    xdat <- get_X()
   
    res <-
      do.call(cbind, furrr::future_map(vecs, function(v) {
        v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
        pls_betas(xdat$X, v0, ncomp = ncomp)
      }))
    list(beta_matrix=as.matrix(res), estimated_hrf=NULL)
  } else if (method == "pls_global") {
    vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) >0))
    xdat <- get_X()
    Y <- do.call(cbind, lapply(vecs, function(v) v))
    
    if (ncomp < log(ncol(Y))) {
      warning("'ncomp' for pls_global method is less than log(nvoxels), consider increasing.")
    }
    
    Y0 <- resid(lsfit(xdat$Base, Y, intercept = FALSE))
    list(beta_matrix=pls_global_betas(xdat$X, Y0, ncomp=ncomp), estimated_hrf=NULL)
  } else if (method == "ols") {
    vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) >0))
    xdat <- get_X()
    Y <- do.call(cbind, lapply(vecs, function(v) v))
    Y0 <- resid(lsfit(xdat$Base, Y, intercept = FALSE))
    list(beta_matrix=ols_betas(xdat$X, Y0), estimated_hrf=NULL)
  } else if (method == "lowrank_hrf") {
    # 1. Find best HRFs with clustering
    hrf_results <- find_best_hrf(dset, onset_var = bdes$emod_ran$model_spec$onset_var, 
                                cluster_series = TRUE)
    
    # 2. For each unique HRF cluster:
    unique_hrfs <- hrf_results$cluster_representatives
    voxel_clusters <- split(seq_len(ncol(get_data(dset))), hrf_results$reassigned_hrfs)
    
    # 3. Create design matrices and estimate betas for each cluster
    beta_list <- lapply(names(unique_hrfs), function(hrf_id) {
      # Get voxels for this cluster
      voxel_indices <- voxel_clusters[[hrf_id]]
      
      # Create event model with this HRF
      emod_ran_hrf <- event_model(
        bdes$emod_ran$model_spec$formula,
        data = dset$event_table,
        block = bdes$emod_ran$model_spec$block,
        sampling_frame = dset$sampling_frame,
        hrf = gen_empirical_hrf(rsam, L[, unique_hrfs[hrf_id]])
      )
      
      # Get design matrix for this HRF
      X_ran <- design_matrix(emod_ran_hrf)
      X <- if (is.null(bdes$dmat_fixed)) {
        cbind(X_ran, bdes$dmat_base)
      } else {
        cbind(X_ran, bdes$dmat_fixed, bdes$dmat_base)
      }
      
      # Use LSS for this cluster's voxels
      lss_fast(dset, bdes, X = X, voxel_indices = voxel_indices)
    })
    
    # 4. Combine results
    beta_matrix <- matrix(0, nrow = nrow(beta_list[[1]]), ncol = ncol(get_data(dset)))
    for (i in seq_along(beta_list)) {
      voxel_indices <- voxel_clusters[[i]]
      beta_matrix[, voxel_indices] <- beta_list[[i]]
    }
    
    list(beta_matrix = beta_matrix, 
         estimated_hrf = hrf_results$L_best_hrfs)
  } else {
    stop("Invalid method. Supported methods are 'r1_glms', 'r1', 'mixed', 'pls', 'pls_global', and 'ols'")
  }
}


#' @keywords internal
gen_beta_design <- function(fixed = NULL, ran, block, bmod, dset, method = NULL) {
  # Get the base design matrices
  if (!is.null(fixed)) {
    emod_fixed <- event_model(fixed, data = dset$event_table, block = block, sampling_frame = dset$sampling_frame)
    dmat_fixed <- design_matrix(emod_fixed)
  } else {
    emod_fixed <- NULL
    dmat_fixed <- NULL
  }
  
  emod_ran <- event_model(ran, data = dset$event_table, block = block, sampling_frame = dset$sampling_frame)
  dmat_ran <- design_matrix(emod_ran)
  dmat_base <- design_matrix(bmod)
  
  # For R1 methods, we need to track both expanded and collapsed indices
  is_r1_method <- !is.null(method) && method %in% c("r1", "r1_glms")
  if (is_r1_method) {
    # Get the number of basis functions from the design matrix
    n_basis <- ncol(dmat_ran) / nrow(emod_ran$model_spec$event_table)
    n_events <- ncol(dmat_ran)/n_basis
  
    # Create both expanded and collapsed indices
    ran_ind <- 1:n_events  # Collapsed indices for R1 methods
    ran_ind_expanded <- 1:(n_events * n_basis)  # Expanded indices for other methods
  } else {
    ran_ind <- 1:ncol(dmat_ran)
    ran_ind_expanded <- ran_ind
  }
  
  # Combine design matrices
  dmat_all <- if (is.null(fixed)) {
    cbind(dmat_ran, dmat_base)
  } else {
    cbind(dmat_ran, dmat_fixed, dmat_base)
  }
  
  # Calculate fixed and base indices
  start_fixed <- ncol(dmat_ran) + 1
  if (is.null(fixed)) {
    start_base <- start_fixed
    fixed_ind <- NULL
  } else {
    start_base <- start_fixed + ncol(dmat_fixed)
    fixed_ind <- start_fixed:(start_base - 1)
  }
  base_ind <- start_base:ncol(dmat_all)
  
  # Return list with both expanded and collapsed indices
  list(
    bmod = bmod,
    emod_fixed = emod_fixed,
    emod_ran = emod_ran,
    dmat_fixed = dmat_fixed,
    dmat_ran = dmat_ran,
    dmat_base = dmat_base,
    ran_ind = ran_ind,
    ran_ind_expanded = ran_ind_expanded,
    fixed_ind = fixed_ind,
    base_ind = base_ind,
    n_basis = if(is_r1_method) n_basis else NULL
  )
}

#' Estimate betas for a matrix dataset
#'
#' This function estimates betas (regression coefficients) for fixed and random effects
#' in a matrix dataset using various methods.
#'
#' @param x An object of class `matrix_dataset` representing the matrix dataset
#' @param fixed A formula specifying the fixed regressors that model constant effects (i.e., non-varying over trials)
#' @param ran A formula specifying the random (trialwise) regressors that model single trial effects
#' @param block A formula specifying the block factor
#' @param method The regression method for estimating trialwise betas; one of "mixed", "pls", "pls_global", or "ols" (default: "mixed")
#' @param basemod A `baseline_model` instance to regress out of data before beta estimation (default: NULL)
#' @param ncomp Number of PLS components for the "pls" and "pls_global" methods (default: 4)
#' @param lambda Lambda parameter (not currently used; default: 0.01)
#' @param ... Additional arguments passed to the estimation method
#' 
#' @family estimate_betas
#'
#' @return A list of class "fmri_betas" containing the following components:
#'   * betas_fixed: Matrix representing the fixed effect betas
#'   * betas_ran: Matrix representing the random effect betas
#'   * design_ran: Design matrix for random effects
#'   * design_fixed: Design matrix for fixed effects
#'   * design_base: Design matrix for baseline model
#'
#' @seealso \code{\link{matrix_dataset}}, \code{\link{baseline_model}}
#'
#' @export
estimate_betas.matrix_dataset <- function(x,fixed=NULL, ran, block,  
                                        method=c("r1_glms", "r1", "lss", "lss_naive", "mixed", "pls", "pls_global", "ols"), 
                                        basemod=NULL,
                                        hrf_basis = NULL,
                                        hrf_ref = NULL,
                                        ncomp=4, lambda=.01,...) {
  
  method <- match.arg(method)
  dset <- x
  mask <- get_mask(dset)
 
  bmod <- if (is.null(basemod)) {
    baseline_model("constant", sframe=dset$sampling_frame)
  } else {
    basemod
  }

  bdes <- gen_beta_design(fixed, ran, block, bmod, dset)
  betas <- run_estimate_betas(bdes, dset, method, ncomp=ncomp)
  
  ran <- as.matrix(betas[bdes$ran_ind,,drop=FALSE])
  fixed <- as.matrix(betas[bdes$fixed_ind,,drop=FALSE])
  
  ret <- list(betas_fixed=fixed,
              betas_ran=ran,
              #design=bdes$dmat_all,
              design_ran=bdes$dmat_ran,
              design_fixed=bdes$dmat_fixed,
              design_base=bdes$dmat_base)
  
  class(ret) <-  c("fmri_betas")
  ret
  
}

#' Estimate betas for a latent dataset
#'
#' This function estimates betas (regression coefficients) for fixed and random effects
#' in a matrix dataset using various methods.
#'
#' @param x An object of class `matrix_dataset` representing the matrix dataset
#' @param fixed A formula specifying the fixed regressors that model constant effects (i.e., non-varying over trials)
#' @param ran A formula specifying the random (trialwise) regressors that model single trial effects
#' @param block A formula specifying the block factor
#' @param method The regression method for estimating trialwise betas; one of "mixed", "pls", "pls_global", or "ols" (default: "mixed")
#' @param basemod A `baseline_model` instance to regress out of data before beta estimation (default: NULL)
#' @param ncomp Number of PLS components for the "pls" and "pls_global" methods (default: 4)
#' @param lambda Lambda parameter (not currently used; default: 0.01)
#' @param prewhiten currently experimental, default to \code{FALSE}.
#' @param ... Additional arguments passed to the estimation method
#'
#' @return A list of class "fmri_betas" containing the following components:
#'   * betas_fixed: Matrix representing the fixed effect betas
#'   * betas_ran: Matrix representing the random effect betas
#'   * design_ran: Design matrix for random effects
#'   * design_fixed: Design matrix for fixed effects
#'   * design_base: Design matrix for baseline model
#'
#' @seealso \code{\link{matrix_dataset}}, \code{\link{baseline_model}}
#'
#' @family estimate_betas
#'
#' @export
#' @rdname estimate_betas
estimate_betas.latent_dataset <- function(x, fixed=NULL, ran, block, 
                                          method=c("mixed", "pls", "pls_global", "ols"), 
                                          basemod=NULL, ncomp=4, lambda=.01, prewhiten=FALSE,...) {
  
  method <- match.arg(method)
  dset <- x
  mask <- get_mask(dset)
  
  bmod <- if (is.null(basemod)) {
    baseline_model("constant", sframe=dset$sampling_frame)
  } else {
    basemod
  }
  
  bdes <- gen_beta_design(fixed, ran, block, bmod, dset)
  
  if (prewhiten) {
    wmat <- auto_whiten(dset@basis, fixed)
    ## hack
    ## swap in whitened matrix
    dset@basis <- wmat
    ###
  }
  
  betas <- run_estimate_betas(bdes, dset, method, ncomp=ncomp)
  
  ran <- as.matrix(betas[bdes$ran_ind,,drop=FALSE])
  fixed <- as.matrix(betas[bdes$fixed_ind,,drop=FALSE])
  
  ret <- list(betas_fixed=fixed,
              betas_ran=ran,
              #design=bdes$dmat_all,
              design_ran=bdes$dmat_ran,
              design_fixed=bdes$dmat_fixed,
              design_base=bdes$dmat_base,
              prewhiten=prewhiten)
  
  class(ret) <-  c("fmri_latent_betas", "fmri_betas")
  ret
}


#' Estimate hemodynamic response function (HRF) using Generalized Additive Models (GAMs)
#'
#' This function estimates the HRF using GAMs from the `mgcv` package.
#' The HRF can be estimated with or without fixed effects.
#'
#' @param form A formula specifying the event model for the conditions of interest
#' @param fixed A formula specifying the fixed regressors that model constant effects (i.e., non-varying over trials); default is NULL
#' @param block A formula specifying the block factor
#' @param dataset An object representing the fMRI dataset
#' @param bs Basis function for the smooth term in the GAM; one of "tp" (default), "ts", "cr", or "ps"
#' @param rsam A sequence of time points at which the HRF is estimated (default: seq(0, 20, by = 1))
#' @param basemod A `baseline_model` instance to regress out of data before HRF estimation (default: NULL)
#' @param k the dimension of the basis, default is 8
#' @param fx indicates whether the term is a fixed d.f. regression spline (TRUE) or a penalized regression spline (FALSE); default is TRUE.
#' 
#' @return A matrix with the estimated HRF values for each voxel
#'
#' @importFrom neuroim2 vectors
#' @importFrom furrr future_map
#' @seealso \code{\link{baseline_model}}, \code{\link{event_model}}, \code{\link{design_matrix}}
#' @examples 
#' 
#' # To be added
#'
#' @export 
#' @autoglobal
estimate_hrf <- function(form, fixed=NULL, block, dataset, 
                           bs=c("tp", "ts", "cr", "ps"), 
                           rsam=seq(0,20,by=1),
                           basemod=NULL,
                           k=8,
                           fx=TRUE) {
  with_package("mgcv")
  dset <- dataset
  bvec <- get_data(dset)
  mask <- get_mask(dset)
  
  onset_var <- lazyeval::f_lhs(form)
  dvars <- lazyeval::f_rhs(form)
  
  bmod <- if (is.null(basemod)) {
    baseline_model("constant", sframe=dset$sampling_frame)
  } else {
    basemod
  }
  
  if (!is.null(fixed)) {
    emod_fixed <- event_model(fixed, data=dset$event_table, block=block, sampling_frame=dset$sampling_frame)
    X_fixed <- as.matrix(design_matrix(emod_fixed))
    has_fixed=TRUE
  } else {
    has_fixed=FALSE
  }
  
  emat_cond <- event_model(form, data=dset$event_table, block=block, 
                           sampling_frame=dset$sampling_frame)
  

  X_base <- as.matrix(design_matrix(bmod))
  X_cond <- as.matrix(design_matrix(emat_cond))
  #browser()
  
  res <- do.call(cbind, furrr::future_map(neuroim2::vectors(bvec, subset=which(mask>0)), function(v) {
    gam.1 <- if (has_fixed) {
      mgcv::gam(v ~ mgcv::s(X_cond, bs=bs, fx=TRUE, k=8) + X_fixed + X_base)
    } else {
      mgcv::gam(v ~ mgcv::s(X_cond, bs=bs, fx=TRUE, k=8) + X_base)
    }
    
    time <- rsam
    xb <- matrix(colMeans(X_base), length(time),ncol(X_base), byrow=TRUE)
    ##xf <- matrix(colMeans(X_fixed), length(time),ncol(X_fixed), byrow=TRUE)
    predict(gam.1, list(X_cond=time, X_base=xb, X_fixed=xf))
  }))
  
  res
  
}
#' 
#' #' @keywords internal
#' #' @noRd
#' setup_hrf_library <- function(hrflib = NULL, rsam = seq(0, 24, by = 1), ncomp = 5) {
#'   if (is.null(hrflib)) {
#'     params <- expand.grid(
#'       h1 = seq(1, 3, by = 0.33),
#'       h2 = seq(3, 7, by = 0.33),
#'       h3 = seq(6, 8, by = 1),
#'       h4 = seq(5, 9, by = 1),
#'       f1 = seq(0, 0.2, by = 0.1),
#'       f2 = seq(0, 0.2, by = 0.1)
#'     )
#'     hrflib <- gen_hrf_library(hrf_half_cosine, params)
#'   }
#'   
#'   L <- hrflib(rsam)
#'   pca_L <- multivarious::pca(L, preproc = multivarious::pass())
#'   Lk <- pca_L$u[, 1:ncomp, drop = FALSE]
#'   
#'   basis_set <- lapply(1:ncol(Lk), function(i) {
#'     gen_empirical_hrf(rsam, Lk[, i])
#'   })
#'   
#'   list(L = L, Lk = Lk, basis_set = basis_set)
#' }
#' 
#' #' @keywords internal
#' #' @noRd
#' compute_best_hrfs <- function(X_proj, L_proj, L, cluster_series = TRUE, cluster_quantile = .5) {
#'   basis_euc <- proxy::dist(L_proj, X_proj, method = "euclidean")
#'   best_match_indices <- apply(as.matrix(basis_euc), 2, which.min)
#'   
#'   if (!cluster_series) {
#'     return(list(
#'       best_hrf_indices = best_match_indices,
#'       L_best_hrfs = L[, best_match_indices, drop = FALSE],
#'       L_proj = L_proj,
#'       X_proj = X_proj
#'     ))
#'   }
#'   
#'   unique_best_hrfs <- unique(best_match_indices)
#'   L_best_hrfs <- L[, unique_best_hrfs, drop = FALSE]
#'   hrf_dist <- dist(t(L_best_hrfs))
#'   hclust_result <- hclust(hrf_dist, method = "average")
#'   clusters <- cutree(hclust_result, h = quantile(hrf_dist, cluster_quantile))
#'   
#'   cluster_medoid <- function(cluster_indices) {
#'     sub_dist <- as.matrix(hrf_dist)[cluster_indices, cluster_indices]
#'     medoid_index <- cluster_indices[which.min(rowSums(sub_dist))]
#'     return(unique_best_hrfs[medoid_index])
#'   }
#'   
#'   unique_clusters <- unique(clusters)
#'   cluster_representatives <- sapply(unique_clusters, function(cluster) {
#'     cluster_indices <- which(clusters == cluster)
#'     cluster_medoid(cluster_indices)
#'   })
#'   
#'   reassigned_hrfs <- sapply(best_match_indices, function(idx) {
#'     cluster <- clusters[which(unique_best_hrfs == idx)]
#'     cluster_representatives[as.character(cluster)]
#'   })
#'   
#'   list(
#'     best_hrf_indices = unique_best_hrfs,
#'     clusters = clusters,
#'     cluster_representatives = cluster_representatives,
#'     reassigned_hrfs = reassigned_hrfs,
#'     L_best_hrfs = L_best_hrfs,
#'     L_proj = L_proj,
#'     X_proj = X_proj
#'   )
#' }
#' 
#' 
#' 
#' #' @export
#' find_best_hrf <- function(dataset, onset_var, hrflib = NULL, rsam = seq(0, 24, by = 1),
#'                          ncomp = 5, nsplits = 3, block = NULL, basemod = NULL, 
#'                          cluster_series = TRUE, 
#'                          cluster_quantile = .5) {
#'   UseMethod("find_best_hrf")
#' }
#' 
#' #' @export
#' find_best_hrf.matrix_dataset <- function(dataset, onset_var, hrflib = NULL, rsam = seq(0, 24, by = 1),
#'                                        ncomp = 5, nsplits = 3, block = NULL, basemod = NULL, 
#'                                        cluster_series = TRUE, 
#'                                        cluster_quantile = .5) {
#'   
#'   # Setup HRF library and basis functions
#'   hrf_setup <- setup_hrf_library(hrflib, rsam, ncomp)
#'   hrf_list <- do.call(gen_hrf_set, hrf_setup$basis_set)
#'   
#'   # Setup models
#'   if (is.null(basemod)) {
#'     bmod <- baseline_model("constant", sframe = dataset$sampling_frame)
#'   } else {
#'     bmod <- basemod
#'   }
#'   
#'   # Create split variable - FIXED THIS LINE
#'   split_ids <- rep(1:nsplits, length.out = nrow(dataset$event_table))
#'   fac_var <- factor(paste0("split_", split_ids))
#'   dataset$event_table <- dataset$event_table %>% mutate(.splitvar = fac_var)
#'   
#'   # Setup event model
#'   form <- as.formula(paste0(onset_var, " ~ ", "hrf(.splitvar, basis = hrf_list)"))
#'   emod <- event_model(
#'     x = form,
#'     formula = form,
#'     data = dataset$event_table,
#'     block = block,
#'     sampling_frame = dataset$sampling_frame
#'   )
#'   
#'   # Get data and compute residuals
#'   X <- get_data_matrix(dataset)
#'   X_base <- as.matrix(design_matrix(bmod))
#'   X_resid <- resid(lsfit(X_base, X, intercept = FALSE))
#'   
#'   # Dimensionality reduction
#'   pca_X <- multivarious::pca(X_resid, preproc = multivarious::center())
#'   ncomp_X <- min(ncomp, ncol(pca_X$u))
#'   Mk <- pca_X$u[, 1:ncomp_X, drop = FALSE]
#'   
#'   # Compute projections
#'   Rk <- design_matrix(emod)
#'   C <- t(Rk) %*% Mk
#'   svd_C <- svd(C)
#'   
#'   X_proj <- t(X_resid) %*% Mk %*% svd_C$v
#'   L_proj <- t(hrf_setup$L) %*% hrf_setup$Lk %*% svd_C$u
#'   
#'   # Compute best HRFs and return results
#'   compute_best_hrfs(X_proj, L_proj, hrf_setup$L, cluster_series, cluster_quantile)
#' }
#' 
#' #' @export
#' find_best_hrf.fmri_mem_dataset <- function(dataset, onset_var, hrflib = NULL, rsam = seq(0, 24, by = 1),
#'                                          ncomp = 5, nsplits = 3, block = NULL, basemod = NULL, 
#'                                          cluster_series = TRUE, 
#'                                          cluster_quantile = .5) {
#'   # Use the same implementation as matrix_dataset but with get_data_matrix
#'   find_best_hrf.matrix_dataset(
#'     dataset = dataset,
#'     onset_var = onset_var,
#'     hrflib = hrflib,
#'     rsam = rsam,
#'     ncomp = ncomp,
#'     nsplits = nsplits,
#'     block = block,
#'     basemod = basemod,
#'     cluster_series = cluster_series,
#'     cluster_quantile = cluster_quantile
#'   )
#' }
#' 
#' 
#' #' @keywords internal
#' #' @noRd
#' smooth_hrf_assignments <- function(reassigned_hrfs, L_best_hrfs, voxel_coords, 
#'                                    lambda = 0.1, n_iterations = 5, k = 6) {
#'   # This function will take the initial HRF assignments and encourage spatial smoothness.
#'   # Steps:
#'   # 1. Build a neighbor structure using voxel_coords.
#'   # 2. Iteratively reassign HRFs to minimize:
#'   #    Cost = distance(voxel's HRF, voxel's data) + lambda * sum(differences to neighbors)
#'   #
#'   # For simplicity, we'll consider "differences to neighbors" as the average difference
#'   # in HRF shape compared to neighbors. HRF shape difference can be approximated by 
#'   # Euclidean distance in the original HRF space (L_best_hrfs) for the assigned HRFs.
#'   
#'   n_voxels <- length(reassigned_hrfs)
#'   if (n_voxels == 0) return(reassigned_hrfs)
#'   
#'   # Build a neighbor index for each voxel
#'   # We'll use a quick k-nearest neighbors approach:
#'   # Note: For large datasets, consider a more efficient approach. Here we do a simple R-based solution.
#'   dists <- as.matrix(dist(voxel_coords))
#'   # For each voxel, find the k nearest neighbors (excluding itself)
#'   neighbor_indices <- apply(dists, 1, function(row) {
#'     # order by distance and take the first k+1 (including self)
#'     # then drop self and take k neighbors
#'     neigh <- order(row)
#'     neigh <- neigh[neigh != which.min(row)] # remove self
#'     head(neigh, k)
#'   })
#'   # neighbor_indices is now a matrix with k columns, each column is for one voxel
#'   
#'   # Precompute HRF shapes for assigned indices
#'   # Reassigned_hrfs are indices into L_best_hrfs
#'   # L_best_hrfs is a matrix [time, hrf_index], reassigned_hrfs indicates which column each voxel gets
#'   # Actually, L_best_hrfs contains only the unique best HRFs. We need to handle differences:
#'   # We'll assume reassigned_hrfs indexes directly into columns of L_best_hrfs.
#'   # If not, ensure that these align. The code from the original snippet does align them.
#'   
#'   # Function to get HRF shape given an index in L_best_hrfs
#'   get_hrf_shape <- function(idx) L_best_hrfs[, idx]
#'   
#'   # On each iteration, update assignments
#'   for (iter in seq_len(n_iterations)) {
#'     # We'll compute a simple smoothness cost:
#'     # For voxel v with HRF h:
#'     # Cost = 0 (data-fitting already done)
#'     #       + lambda * average distance to neighbors' HRFs
#'     # We'll reassign h to the closest HRF among the currently chosen set of HRFs (if we have multiple HRFs?).
#'     # In this simplified scenario, we only have one chosen HRF per voxel from the cluster representatives.
#'     # To truly allow refinement, we might consider allowing a small set of candidate HRFs.
#'     #
#'     # For now, let's just encourage consistency. We won't re-choose from the entire library (that would be costly).
#'     # Instead, we'll nudge voxels to switch to an HRF that reduces the smoothness cost if there's an alternative HRF
#'     # among neighbors that might be better.
#'     #
#'     # A simple heuristic:
#'     # For each voxel, look at neighbors' HRFs, pick the HRF that is most frequent among neighbors
#'     # and if that HRF is different and reduces average distance cost, switch to it.
#'     
#'     new_assignments <- reassigned_hrfs
#'     for (v in seq_len(n_voxels)) {
#'       v_hrf_idx <- reassigned_hrfs[v]
#'       neigh_idx <- neighbor_indices[, v]
#'       neigh_hrfs <- reassigned_hrfs[neigh_idx]
#'       # Find the most common neighbor HRF
#'       hrf_freq <- sort(table(neigh_hrfs), decreasing = TRUE)
#'       # Consider top candidate from neighbors
#'       candidate_hrf <- as.numeric(names(hrf_freq)[1])
#'       
#'       if (candidate_hrf != v_hrf_idx) {
#'         # Check if switching reduces cost
#'         # Current cost: distance of v_hrf to neighbors
#'         v_shape <- get_hrf_shape(v_hrf_idx)
#'         neigh_shapes <- L_best_hrfs[, neigh_hrfs, drop = FALSE]
#'         current_cost <- mean(apply(neigh_shapes, 2, function(x) sqrt(sum((v_shape - x)^2))))
#'         
#'         # Candidate cost
#'         cand_shape <- get_hrf_shape(candidate_hrf)
#'         candidate_cost <- mean(apply(neigh_shapes, 2, function(x) sqrt(sum((cand_shape - x)^2))))
#'         
#'         # If candidate improves by at least lambda (some threshold),
#'         # or simply if candidate_cost < current_cost (since lambda is already controlling how strongly we weigh this)
#'         if (candidate_cost < current_cost) {
#'           new_assignments[v] <- candidate_hrf
#'         }
#'       }
#'     }
#'     reassigned_hrfs <- new_assignments
#'   }
#'   
#'   reassigned_hrfs
#' }
#' 
#' 
#' #' #' @export
#' #' find_best_hrf <- function(dataset, onset_var, hrflib = NULL, rsam = seq(0, 24, by = 1),
#' #'                           ncomp = 5, nsplits = 3, block = NULL, basemod = NULL, 
#' #'                           cluster_series = TRUE, 
#' #'                           cluster_quantile = .5,
#' #'                           dist_method = "euclidean",
#' #'                           spatial_smoothing = FALSE,
#' #'                           lambda = 0.1,
#' #'                           n_iterations = 5,
#' #'                           k_neighbors = 6) {
#' #'   UseMethod("find_best_hrf")
#' #' }
#' 
#' #' @export
#' find_best_hrf.matrix_dataset <- function(dataset, onset_var, hrflib = NULL, rsam = seq(0, 24, by = 1),
#'                                          ncomp = 5, nsplits = 3, block = NULL, basemod = NULL, 
#'                                          cluster_series = TRUE, 
#'                                          cluster_quantile = .5,
#'                                          dist_method = "euclidean",
#'                                          spatial_smoothing = FALSE,
#'                                          lambda = 0.1,
#'                                          n_iterations = 5,
#'                                          k_neighbors = 6) {
#'   
#'   # Setup HRF library and basis functions
#'   hrf_setup <- setup_hrf_library(hrflib, rsam, ncomp)
#'   hrf_list <- do.call(gen_hrf_set, hrf_setup$basis_set)
#'   
#'   # Setup models
#'   if (is.null(basemod)) {
#'     bmod <- baseline_model("constant", sframe = dataset$sampling_frame)
#'   } else {
#'     bmod <- basemod
#'   }
#'   
#'   # Create split variable
#'   split_ids <- rep(1:nsplits, length.out = nrow(dataset$event_table))
#'   fac_var <- factor(paste0("split_", split_ids))
#'   dataset$event_table <- dataset$event_table %>% mutate(.splitvar = fac_var)
#'   
#'   # Setup event model
#'   form <- as.formula(paste0(onset_var, " ~ ", "hrf(.splitvar, basis = hrf_list)"))
#'   emod <- event_model(
#'     x = form,
#'     formula = form,
#'     data = dataset$event_table,
#'     block = block,
#'     sampling_frame = dataset$sampling_frame
#'   )
#'   
#'   # Get data and compute residuals
#'   X <- get_data_matrix(dataset)
#'   X_base <- as.matrix(design_matrix(bmod))
#'   X_resid <- resid(lsfit(X_base, X, intercept = FALSE))
#'   
#'   # Dimensionality reduction
#'   pca_X <- multivarious::pca(X_resid, preproc = multivarious::center())
#'   ncomp_X <- min(ncomp, ncol(pca_X$u))
#'   Mk <- pca_X$u[, 1:ncomp_X, drop = FALSE]
#'   
#'   # Compute projections
#'   Rk <- design_matrix(emod)
#'   C <- t(Rk) %*% Mk
#'   svd_C <- svd(C)
#'   
#'   X_proj <- t(X_resid) %*% Mk %*% svd_C$v
#'   L_proj <- t(hrf_setup$L) %*% hrf_setup$Lk %*% svd_C$u
#'   
#'   # Compute best HRFs using chosen distance measure
#'   res <- compute_best_hrfs(X_proj, L_proj, hrf_setup$L, cluster_series, cluster_quantile, dist_method = dist_method)
#'   
#'   # If spatial smoothing is requested and voxel coordinates are available
#'   if (spatial_smoothing && !is.null(dataset$voxel_coords)) {
#'     # We have res$reassigned_hrfs or res$best_hrf_indices
#'     # Use res$reassigned_hrfs if available; if not, default to best_hrf_indices.
#'     # We'll also need L_best_hrfs from res. The "res" list has L_best_hrfs and final assignments.
#'     # We'll smooth over the reassigned_hrfs.
#'     final_assignments <- if (!is.null(res$reassigned_hrfs)) res$reassigned_hrfs else res$best_hrf_indices
#'     
#'     final_assignments_smoothed <- smooth_hrf_assignments(
#'       reassigned_hrfs = final_assignments,
#'       L_best_hrfs = res$L_best_hrfs,
#'       voxel_coords = dataset$voxel_coords,
#'       lambda = lambda,
#'       n_iterations = n_iterations,
#'       k = k_neighbors
#'     )
#'     
#'     # Update results with smoothed assignments
#'     res$reassigned_hrfs_smoothed <- final_assignments_smoothed
#'   }
#'   
#'   res
#' }
#' 
#' #' @export
#' find_best_hrf.fmri_mem_dataset <- function(dataset, onset_var, hrflib = NULL, rsam = seq(0, 24, by = 1),
#'                                            ncomp = 5, nsplits = 3, block = NULL, basemod = NULL, 
#'                                            cluster_series = TRUE, 
#'                                            cluster_quantile = .5,
#'                                            dist_method = "euclidean",
#'                                            spatial_smoothing = FALSE,
#'                                            lambda = 0.1,
#'                                            n_iterations = 5,
#'                                            k_neighbors = 6) {
#'   find_best_hrf.matrix_dataset(
#'     dataset = dataset,
#'     onset_var = onset_var,
#'     hrflib = hrflib,
#'     rsam = rsam,
#'     ncomp = ncomp,
#'     nsplits = nsplits,
#'     block = block,
#'     basemod = basemod,
#'     cluster_series = cluster_series,
#'     cluster_quantile = cluster_quantile,
#'     dist_method = dist_method,
#'     spatial_smoothing = spatial_smoothing,
#'     lambda = lambda,
#'     n_iterations = n_iterations,
#'     k_neighbors = k_neighbors
#'   )
#' }



