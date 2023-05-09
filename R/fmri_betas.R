

#' @importFrom glmnet glmnet
#' @noRd
#' @keywords internal
ridge_betas <- function(X, Y, penalty_factor=rep(1:ncol(X)), lambda=.01) {
  fit <- glmnet::glmnet(X, Y, penalty.factor=penalty_factor, alpha=0,lambda=lambda)
  coef(fit)[,1,drop=FALSE]
}

#' @importFrom pls plsr 
#' @noRd
pls_betas <- function(X, Y, ncomp=3) {
  dx <- data.frame(X=X, Y=Y)
  fit <- pls::plsr(Y ~ X, data=dx, ncomp=ncomp, method="simpls", scale=TRUE, maxit=500)
  coef(fit, ncomp=ncomp)[,,1]
}


#' @keywords internal
#' @noRd
pls_global_betas <- function(X, Y, ncomp=3) {
  dx <- data.frame(X=X, Y=Y)
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
  slm.1 <- care::slm(X, Y, verbose=FALSE)
  b2 <- coef(slm.1)[,-(1:2)]
  b1 <- coef(slm.1)[,1]
  b1 + b2
}


#' @noRd
#' @keywords internal
mixed_betas <- function(X, Y, ran_ind, fixed_ind) {
  fit <- rrBLUP::mixed.solve(Y, Z=X[,ran_ind], X=X[,c(fixed_ind)], bounds=c(c(1e-05, .2)))
  c(fit$u, fit$b)
}


#' @keywords internal
gen_beta_design <- function(fixed=NULL, ran, block, bmod, dset) {
  
  if (!is.null(fixed)) {
    emod_fixed <- event_model(fixed, data=dset$event_table, block=block, sampling_frame=dset$sampling_frame)
    dmat_fixed <- design_matrix(emod_fixed)
  } else {
    emod_fixed=NULL
    dmat_fixed=NULL
  }
  
  emod_ran <- event_model(ran, data=dset$event_table, block=block, sampling_frame=dset$sampling_frame)
  dmat_ran <- design_matrix(emod_ran)
  
  dmat_base <- design_matrix(bmod)
  #dmat_all <- cbind(scale(dmat_ran), scale(dmat_fixed), scale(dmat_base))
  
  
  dmat_all <- if (is.null(fixed)) {
    cbind(dmat_ran, dmat_base)
  } else {
    cbind(dmat_ran, dmat_fixed, dmat_base)
  }
  
  
  start_ran <- 1
  start_fixed <- ncol(dmat_ran)+1
  
  if (is.null(fixed)) {
    start_base <- start_fixed 
    fixed_ind <- NULL
  } else {
 
    start_base <- start_fixed + ncol(dmat_fixed)
    fixed_ind <- start_fixed:(start_base-1)
  }
  
  ran_ind <- 1:ncol(dmat_ran)
  base_ind <- start_base:(ncol(dmat_all))
  
  
  list(bmod=bmod,
       emod_fixed=emod_fixed,
       emod_ran=emod_ran,
       dmat_fixed=dmat_fixed,
       dmat_ran=dmat_ran,
       dmat_base=dmat_base,
       ran_ind=ran_ind,
       fixed_ind=fixed_ind,
       base_ind=base_ind)
}
  
#' @importFrom furrr future_map
run_estimate_betas <- function(bdes, dset, method, ncomp=3, niter=8, radius=8) {
  
  get_X <- function(scale=FALSE) {
    X <- if (is.null(bdes$fixed)) bdes$dmat_ran else cbind(bdes$dmat_ran, bdes$dmat_fixed)
    #X  <- cbind(bdes$dmat_ran, bdes$dmat_fixed)
    if (scale) {
      X <- scale(X)
    }
    
    Base <- as.matrix(bdes$dmat_base)
    X[is.na(X)] <- 0
    list(Base=Base, X=X)
  }

  betas <- with(bdes, {
    if (method == "slm") {
      vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) >0))
      ## TODO FIXME
      xdat <- get_X()
      
      res <- do.call(cbind, furrr::future_map(vecs, function(v) { 
        v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
        slm_betas(xdat$X, v0) 
      }))
      as.matrix(res)
      
      
    } else if (method == "mixed") {
      vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) > 0))
      xdat <- get_X()
      res <-
        do.call(cbind, furrr::future_map(vecs, function(v) {
          v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
          mixed_betas(xdat$X, v0, ran_ind = ran_ind, fixed_ind = fixed_ind)
        }))
      as.matrix(res)
    } else if (method == "pls") {
      vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) >0))
      xdat <- get_X(scale=TRUE)
      
      res <-
        do.call(cbind, furrr::future_map(vecs, function(v) {
          v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
          pls_betas(xdat$X, v0, ncomp = ncomp)
        }))
      as.matrix(res)
    } else if (method == "pls_global") {
      vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) >0))
      xdat <- get_X(scale=TRUE)
      Y <- do.call(cbind, lapply(vecs, function(v) v))
      
      if (ncomp < log(ncol(Y))) {
        warning("'ncomp' for pls_global method is less than log(nvoxels), consider increasing.")
      }
      
      Y0 <- resid(lsfit(xdat$Base, Y, intercept = FALSE))
      pls_global_betas(xdat$X, Y0, ncomp=ncomp)
    } else if (method == "ols") {
     
      vecs <- neuroim2::vectors(get_data(dset), subset = which(get_mask(dset) >0))
      xdat <- get_X()
      Y <- do.call(cbind, lapply(vecs, function(v) v))
      
      Y0 <- resid(lsfit(xdat$Base, Y, intercept = FALSE))
      ols_betas(xdat$X, Y0)
      
    } else {
     
      xdat <- get_X(scale=TRUE)
      
      bvec <- get_data(dset)
      mask <- get_mask(dset)
      res <- Reduce("+", furrr::future_map(1:niter, function(iter) {
        slight <- neuroim2::random_searchlight(mask, radius = radius)
        mset <- Reduce("+", furrr::future_map(slight, function(s) {
          cds <- neuroim2::coords(s)
          Y <- neuroim2::series(bvec, cds)
          Y0 <- resid(lsfit(xdat$Base, Y, intercept = FALSE))
          
          dx <- list(Y = Y0, X = xdat$X)
          res <-
            pls::plsr(
              Y ~ X,
              data = dx,
              ncomp = ncomp,
              validation = "none",
              method = "simpls",
              maxit=500
            )
          idx <- neuroim2::grid_to_index(mask, cds)
          B <- coef(res)[, , 1, drop = FALSE]
          m <- Matrix::sparseMatrix(
            i = rep(1:nrow(B), length(idx)),
            j = rep(idx, each = nrow(B)),
            x = as.numeric(B),
            dims = c(nrow(B), prod(dim(mask)))
          )
          #com <- colMeans(cds)
          #D <- sqrt(rowSums(sweep(cds, 2, com, "-")^2))
          m
          
        }))
      })) / niter
      
      as.matrix(res)[, mask != 0]
    }
  })
  
  
}

#' 
#' #' @export
#' estimate_betas.fmri_mem_dataset <- function(x,fixed=NULL, ran, block,  
#'                                          method=c("mixed", "pls", "pls_searchlight", "pls_global", "ols"), 
#'                                          basemod=NULL, 
#'                                          radius=8, niter=8, ncomp=4, lambda=.01,...) {
#'   
#'   estimate_betas.fmri_dataset(x,fixed,ran,block, method, basemod, radius, niter,ncomp, lambda,...)
#' }

#' Estimate betas for an fMRI dataset
#'
#' This function estimates betas (regression coefficients) for fixed and random effects
#' in an fMRI dataset using various methods.
#'
#' @param x An object of class `fmri_dataset` representing the fMRI dataset
#' @param fixed A formula specifying the fixed regressors that model constant effects (i.e., non-varying over trials)
#' @param ran A formula specifying the random (trialwise) regressors that model single trial effects
#' @param block A formula specifying the block factor
#' @param method The regression method for estimating trialwise betas; one of "mixed", "pls", "pls_searchlight", "pls_global", or "ols" (default: "mixed")
#' @param basemod A `baseline_model` instance to regress out of data before beta estimation (default: NULL)
#' @param radius The radius in mm for the `pls_searchlight` approach (default: 8)
#' @param niter Number of searchlight iterations for the "pls_searchlight" method (default: 8)
#' @param ncomp Number of PLS components for the "pls", "pls_searchlight", and "pls_global" methods (default: 4)
#' @param lambda Lambda parameter (not currently used; default: 0.01)
#' @param ... Additional arguments passed to the estimation method
#'
#' @return A list of class "fmri_betas" containing the following components:
#'   * betas_fixed: NeuroVec object representing the fixed effect betas
#'   * betas_ran: NeuroVec object representing the random effect betas
#'   * design_ran: Design matrix for random effects
#'   * design_fixed: Design matrix for fixed effects
#'   * design_base: Design matrix for baseline model
#'   * basemod: Baseline model object
#'   * fixed_model: Fixed effect model object
#'   * ran_model: Random effect model object
#'
#' @seealso \code{\link{fmri_dataset}}, \code{\link{baseline_model}}
#'
#' @importFrom care slm
#' @export
#'
#' @details The `method` argument allows for several beta estimation approaches:
#'   * "mixed": Uses a linear mixed-effects modeling of trialwise random effects as implemented in the `rrBLUP` package.
#'   * "pls": Uses separate partial least squares for each voxel to estimate trialwise betas.
#'   * "pls_searchlight": Estimates PLS solutions over searchlight windows and averages the beta estimates.
#'   * "pls_global": Estimates a single multiresponse PLS solution, where the `Y` matrix is the full data matrix.
#'   * "ols": Ordinary least squares estimate of betas – no regularization.
#'
#' @examples 
#' facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
#' facedes$frun <- factor(facedes$run)
#' scans <- paste0("rscan0", 1:6, ".nii")
#'
#' dset <- fmri_dataset(scans=scans,mask="mask.nii",TR=1.5,run_length=rep(436,6),event_table=facedes)
#' fixed = onset ~ hrf(run)
#' ran = onset ~ trialwise()
#' block = ~ run
#'
#' betas <- estimate_betas(dset, fixed=fixed, ran=ran, block=block 
estimate_betas.fmri_dataset <- function(x,fixed=NULL, ran, block,  
                           method=c("mixed", "pls", "pls_searchlight", "pls_global", "ols"), 
                           basemod=NULL, 
                           radius=8, niter=8, ncomp=4, lambda=.01,...) {
  
  method <- match.arg(method)
  dset <- x
  bvec <- get_data(dset)
  mask <- get_mask(dset)
  
  if (method == "mixed" && is.null(fixed)) {
    stop("method 'mixed' requires a fixed effects term.")
  }
  
  bmod <- if (is.null(basemod)) {
    baseline_model("constant", sframe=dset$sampling_frame)
  } else {
    basemod
  }
  
  bdes <- gen_beta_design(fixed, ran, block, bmod, dset)
  betas <- run_estimate_betas(bdes, dset, method, ncomp=ncomp, niter=niter, radius=radius)
  nbetas <- nrow(betas)
  ospace_ran <- neuroim2::add_dim(neuroim2::space(mask), length(bdes$ran_ind))
  if (!is.null(bdes$fixed_ind)) {
    ospace_fixed <- neuroim2::add_dim(neuroim2::space(mask), length(bdes$fixed_ind))
    fixed <- neuroim2::NeuroVec(as.matrix(betas[bdes$fixed_ind,,drop=FALSE]), ospace_fixed, mask=mask)
  } else {
    fixed <- NULL
  }
  
  ran <- neuroim2::NeuroVec(as.matrix(betas[bdes$ran_ind,,drop=FALSE]), ospace_ran, mask=mask)
  
  #baseline <- neuroim2::NeuroVec(as.matrix(betas[base_ind,,drop=FALSE]), ospace, mask=mask)
  
  ret <- list(betas_fixed=fixed,
       betas_ran=ran,
       design_ran=bdes$dmat_ran,
       design_fixed=bdes$dmat_fixed,
       design_base=bdes$dmat_base,
       basemod=basemod,
       fixed_model=bdes$emod_fixed,
       ran_model=bdes$emod_ran)
  
  class(ret) <- "fmri_betas"
  ret
  
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
                                        method=c("mixed", "pls", "pls_global", "ols"), 
                                        basemod=NULL,
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
                                          basemod=NULL, ncomp=4, lambda=.01, prewhiten=TRUE,...) {
  
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
#'
#' @return A matrix with the estimated HRF values for each voxel
#'
#' @importFrom mgcv gam s
#' @importFrom neuroim2 vectors
#' @importFrom furrr future_map
#' @seealso \code{\link{baseline_model}}, \code{\link{event_model}}, \code{\link{design_matrix}}
#' @examples 
#' 
#' # To be added
#'
#' @export 
estimate_hrf <- function(form, fixed=NULL, block, dataset, 
                           bs=c("tp", "ts", "cr", "ps"), 
                           rsam=seq(0,20,by=1),
                           basemod=NULL) {
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
      gam(v ~ s(X_cond, bs=bs) + X_fixed + X_base)
    } else {
      gam(v ~ s(X_cond, bs=bs) + X_base)
    }
    
    time <- rsam
    xb <- matrix(colMeans(X_base), length(time),ncol(X_base), byrow=TRUE)
    predict(gam.1, list(X_cond=time, X_base=xb))
  }))
  
  res
  
}

#' Pretty print method for fmri_betas objects
#'
#' @param x An object of class "fmri_betas"
#' @param ... Additional arguments passed to the print function
#'
#' @export
print.fmri_betas <- function(x, ...) {
  cat("fmri_betas object:\n\n")
  cat("Fixed effect betas:\n")
  print(x$betas_fixed, ...)
  cat("\nRandom effect betas:\n")
  print(x$betas_ran, ...)
  
  cat("\nRandom effects design matrix:\n")
  print(x$design_ran, ...)
  cat("\nFixed effects design matrix:\n")
  print(x$design_fixed, ...)
  cat("\nBaseline design matrix:\n")
  print(x$design_base, ...)
  
  cat("\nBaseline model:\n")
  print(x$basemod, ...)
  cat("\nFixed effect model:\n")
  print(x$fixed_model, ...)
  cat("\nRandom effect model:\n")
  print(x$ran_model, ...)
}

  
  

  
  
  
                             

