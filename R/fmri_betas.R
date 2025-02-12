#' @noRd
#' @keywords internal
ridge_betas <- function(X, Y, penalty_factor=rep(1,ncol(X)), lambda=.01) {
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
                                                   "r1", "pls",  "pls_global", "ols", "lowrank_hrf"),
                                        basemod = NULL,
                                        hrf_basis = NULL,
                                        hrf_ref = NULL,
                                        maxit = 100,
                                        ...) {
  method <- match.arg(method)
  dset <- x
  bvec <- get_data(dset)
  mask <- get_mask(dset)
  
  if (method == "r1") {
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
  betas <- run_estimate_betas(bdes, dset, method, hrf_basis = hrf_basis, hrf_ref = hrf_ref, block=block, maxit = maxit, ...)
  
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
run_estimate_betas <- function(bdes, dset, method, hrf_basis = NULL, hrf_ref = NULL, block, maxit = 100, ncomp=4, ...) {
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
    rsam <- seq(0, 24, by=.25)
    # 1. Find best HRFs with clustering
    hrf_results <- find_best_hrf(dset, onset_var = rlang::f_lhs(bdes$emod_ran$model_spec$formula), 
                                 block=block,
                                 rsam=rsam,
                                cluster_series = TRUE)

    
    # 2. For each unique HRF cluster:
    unique_hrfs <- hrf_results$cluster_representatives
    names(unique_hrfs) <- as.character(seq_along(unique_hrfs))
    voxel_clusters <- split(seq_along(hrf_results$reassigned_hrfs), hrf_results$reassigned_hrfs)
    names(voxel_clusters) <- as.character(seq_along(voxel_clusters))
    L <- hrf_results$L_unique

    Xdat <- get_data_matrix(dset)
    
    #voxel_clusters <- hrf_results$clusters
    # 3. Create design matrices and estimate betas for each cluster
    beta_list <- lapply(names(unique_hrfs), function(hrf_id) {
      # Get voxels for this cluster
      voxel_indices <- voxel_clusters[[as.character(hrf_id)]]

      
      emphrf <- gen_empirical_hrf(rsam, L[, as.integer(hrf_id)])
      new_form <- inject_basis(bdes$emod_ran$model_spec$formula, emphrf)
      # Create event model with this HRF
      emod_ran_hrf <- event_model(
        new_form,
        data = dset$event_table,
        block = bdes$emod_ran$block,
        sampling_frame = dset$sampling_frame,
    
      )
      
      # Get design matrix for this HRF
      X_ran <- design_matrix(emod_ran_hrf)
      X <- if (is.null(bdes$dmat_fixed)) {
        cbind(X_ran, bdes$dmat_base)
      } else {
        cbind(X_ran, bdes$dmat_fixed, bdes$dmat_base)
      }
     
      # Use LSS for this cluster's voxels
      lss_fast(dset, bdes, Y = Xdat[,voxel_indices])
    })


    # 4. Combine results
    beta_matrix <- matrix(0, nrow = nrow(beta_list[[1]]), ncol = length(unlist(voxel_clusters)))
    for (i in seq_along(beta_list)) {
      voxel_indices <- voxel_clusters[[i]]
      beta_matrix[, voxel_indices] <- beta_list[[i]]
    }
    
    list(beta_matrix = beta_matrix, estimated_hrf = NULL)
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



#' @noRd 
#' @keywords internal
inject_basis <- function(oldform, new_basis, fun_names = c("hrf", "trialwise")) {
  stopifnot(is.formula(oldform))
  
  # A recursive helper that descends through an expression
  # and injects `basis=new_basis` into calls named in fun_names.
  recfun <- function(expr) {
    if (!is_call(expr)) {
      return(expr)  # If it’s not a call, return as is
    }
    thisfun <- call_name(expr)
    
    # If this call is one of the functions we want to modify:
    if (thisfun %in% fun_names) {
      # 1) Recursively transform sub-expressions
      expr_args <- as.list(expr)
      for (i in seq_along(expr_args)[-1]) {
        expr_args[[i]] <- recfun(expr_args[[i]])
      }
      # 2) Rebuild the call, then override/add `basis = new_basis`
      call_rebuilt <- as.call(expr_args)
      call_modified <- call_modify(call_rebuilt, basis = new_basis)
      return(call_modified)
    } else {
      # Not hrf() or trialwise(), so keep walking
      expr_args <- as.list(expr)
      for (i in seq_along(expr_args)[-1]) {
        expr_args[[i]] <- recfun(expr_args[[i]])
      }
      return(as.call(expr_args))
    }
  }
  
  # Extract old LHS, RHS, and environment
  lhs     <- f_lhs(oldform)
  rhs_old <- f_rhs(oldform)
  f_env   <- f_env(oldform)
  
  # Recursively transform the RHS
  rhs_new <- recfun(rhs_old)
  
  # Build the new formula with the same environment
  newform <- new_formula(lhs = lhs, rhs = rhs_new, env = f_env)
  newform
}




