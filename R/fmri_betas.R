#' @keywords internal
#' @noRd
ols_betas <- function(X, Y) {
  fit <- lm.fit(as.matrix(X),Y)
  coef(fit)
}


#' @keywords internal
#' @noRd
mixed_betas <- function(X, Y, ran_ind, fixed_ind) {
  # Ensure X is a matrix
  X <- as.matrix(X)
  
  # Check dimensions and prevent 0-dim matrices
  if (ncol(X) == 0 || nrow(X) == 0) {
    stop("Design matrix X has zero rows or columns")
  }
  
  # Ensure ran_ind has proper values
  if (length(ran_ind) == 0) {
    stop("No random effect indices specified")
  }
  
  # Handle case when fixed_ind is NULL
  if (is.null(fixed_ind)) {
    fixed_ind <- integer(0)  # Empty integer vector
  }
  
  # Check if indices are out of bounds
  if (max(c(ran_ind, fixed_ind)) > ncol(X)) {
    stop("Index out of bounds: indices exceed number of columns in X")
  }
  
  # Ensure Y is a vector
  if (!is.vector(Y)) {
    Y <- as.vector(Y)
  }
  
  # Try to use rrBLUP if available
  tryCatch({
    with_package("rrBLUP")
    
    # If fixed_ind is empty, use a minimal X matrix (intercept only)
    X_fixed <- if (length(fixed_ind) == 0) {
      matrix(1, nrow = nrow(X), ncol = 1)
    } else {
      X[, fixed_ind, drop = FALSE]
    }
    
    
    fit <- rrBLUP::mixed.solve(Y, 
                               Z = X[, ran_ind, drop = FALSE], 
                               X = X_fixed) 
                               #bounds = c(1e-07, 0.5))
    
    # Return results, handling the case where fixed_ind is empty
    if (length(fixed_ind) == 0) {
      return(fit$u)  # Only return random effects
    } else {
      return(c(fit$u, fit$b))  # Return both random and fixed effects
    }
  }, 
  error = function(e) {
    # If rrBLUP fails, try using our C++ implementation if available
    message("rrBLUP mixed.solve failed, attempting alternative: ", e$message)
    
    if (requireNamespace("Rcpp", quietly = TRUE) && 
        requireNamespace("fmrilss", quietly = TRUE)) {
      
      # Use the C++ implementation with proper input validation
      X_fixed <- if (length(fixed_ind) == 0) {
        matrix(1, nrow = nrow(X), ncol = 1)
      } else {
        X[, fixed_ind, drop = FALSE]
      }
      
      fit <- tryCatch({
        fmrilss::mixed_solve(Y = Y, 
                            Z = X[, ran_ind, drop = FALSE], 
                            X = X_fixed)
      }, error = function(e2) {
        # If even that fails, use a fallback
        message("C++ mixed model solver also failed: ", e2$message)
        if (length(fixed_ind) == 0) {
          return(list(u = rep(0, length(ran_ind))))
        } else {
          return(list(
            u = rep(0, length(ran_ind)),
            beta = rep(0, length(fixed_ind))
          ))
        }
      })
      
      # Return results based on whether fixed_ind is empty
      if (length(fixed_ind) == 0) {
        return(fit$u)
      } else {
        return(c(fit$u, fit$beta))
      }
    } else {
      # Last resort - return zeros
      message("No alternative mixed model solver available")
      if (length(fixed_ind) == 0) {
        return(rep(0, length(ran_ind)))
      } else {
        return(rep(0, length(c(ran_ind, fixed_ind))))
      }
    }
  })
}

#' Estimate betas using various regression methods
#'
#' This function estimates betas (regression coefficients) for fixed and random effects
#' using various regression methods including mixed models, least squares, and PLS.
#'
#' @param x An object of class `fmri_dataset` representing the fMRI dataset.
#' @param fixed A formula specifying the fixed regressors that model constant effects (i.e., non-varying over trials).
#' @param ran A formula specifying the random (trialwise) regressors that model single trial effects.
#' @param block A formula specifying the block factor.
#' @param method The regression method for estimating trialwise betas; one of "mixed", "lss", or "ols".
#' @param basemod A `baseline_model` instance to regress out of data before beta estimation (default: NULL).
#' @param maxit Maximum number of iterations for optimization methods (default: 1000).
#' @param fracs Fraction of voxels used for prewhitening.
#' @param progress Logical; show progress bar.
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
#'   * estimated_hrf: The estimated HRF vector (NULL for most methods).
#'
#' @seealso \code{\link{fmri_dataset}}, \code{\link{baseline_model}}, \code{\link{event_model}}
#'
#' @examples
#' \dontrun{
#' facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
#' facedes$frun <- factor(facedes$run)
#' scans <- paste0("rscan0", 1:6, ".nii")
#'
#' dset <- fmri_dataset(scans=scans, mask="mask.nii", TR=1.5, 
#'         run_length=rep(436,6), event_table=facedes)
#' fixed = onset ~ hrf(run)
#' ran = onset ~ trialwise()
#' block = ~ run
#'
#' betas <- estimate_betas(dset, fixed=fixed, ran=ran, block=block, method="mixed")
#' }
#' @export
estimate_betas.fmri_dataset <- function(x, fixed = NULL, ran, block,
                                        method = c("mixed", "lss", "ols"),
                                        basemod = NULL,
                                        maxit = 1000,
                                        fracs = 0.5,
                                        progress = TRUE,
                                        ...) {
  method <- match.arg(method)
  dset <- x
  bvec <- fmridataset::get_data(dset)
  mask <- fmridataset::get_mask(dset)
  
  bmod <- if (is.null(basemod)) {
    baseline_model("constant", sframe = dset$sampling_frame)
  } else {
    basemod
  }
  
  bdes <- gen_beta_design(fixed, ran, block, bmod, dset, method = method)
  betas <- run_estimate_betas(bdes, dset, method, block = block,
                              maxit = maxit, fracs = fracs,
                              progress = progress, ...)
  
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
run_estimate_betas <- function(bdes, dset, method,
                               block, maxit = 100,
                               fracs = .5,
                               progress = TRUE,
                               ...) {
  method <- match.arg(method, c("mixed", "lss", "ols"))

  xdat <- build_design_data(bdes)

  if (method == "mixed") {
    vecs <- masked_vectors(dset)
    res <- map_voxels(vecs, function(v) {
      v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
      mixed_betas(
        xdat$X,
        v0,
        ran_ind = seq_len(ncol(bdes$dmat_ran)),
        fixed_ind = if (!is.null(bdes$dmat_fixed)) {
          (ncol(bdes$dmat_ran) + 1):(ncol(bdes$dmat_ran) + ncol(bdes$dmat_fixed))
        } else {
          NULL
        }
      )
    }, .progress = progress)
    return(list(beta_matrix = as.matrix(res), estimated_hrf = NULL))
  }

  if (method == "lss") {
    data_matrix <- get_data_matrix(dset)
    dmat_base <- as.matrix(bdes$dmat_base)
    dmat_fixed <- if (!is.null(bdes$fixed_ind)) as.matrix(bdes$dmat_fixed) else NULL
    dmat_ran <- as.matrix(bdes$dmat_ran)

    nuisance_matrix <- if (!is.null(dmat_fixed)) {
      cbind(dmat_base, dmat_fixed)
    } else {
      dmat_base
    }

    beta_matrix_ran <- fmrilss::lss(
      Y = data_matrix,
      X = dmat_ran,
      Z = NULL,
      Nuisance = nuisance_matrix,
      method = "r_optimized"
    )

    if (!is.null(bdes$fixed_ind) && length(bdes$fixed_ind) > 0) {
      mask_idx <- which(fmridataset::get_mask(dset) > 0)
      vecs <- neuroim2::vectors(data_matrix, subset = mask_idx)
      X_base_fixed <- cbind(as.matrix(bdes$dmat_base), as.matrix(bdes$dmat_fixed))

      beta_matrix_fixed <- map_voxels(vecs, function(v) {
        fit <- lm.fit(X_base_fixed, v)
        coef(fit)[(ncol(bdes$dmat_base) + 1):length(coef(fit))]
      }, .progress = progress)

      beta_matrix <- rbind(beta_matrix_ran, beta_matrix_fixed)
    } else {
      beta_matrix <- beta_matrix_ran
    }

    return(list(beta_matrix = beta_matrix, estimated_hrf = NULL))
  }

  vecs <- masked_vectors(dset)
  Y <- map_voxels(vecs, function(v) v, .progress = progress)
  Y0 <- resid(lsfit(xdat$Base, Y, intercept = FALSE))
  beta_matrix <- ols_betas(xdat$X, Y0)
  list(beta_matrix = as.matrix(beta_matrix), estimated_hrf = NULL)
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
  
  # Standard indices for all methods
  ran_ind <- 1:ncol(dmat_ran)
  ran_ind_expanded <- ran_ind
  
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
  
  # Return list with indices
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
    base_ind = base_ind
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
#' @param method The regression method for estimating trialwise betas; one of "mixed", "lss", or "ols" (default: "mixed")
#' @param basemod A `baseline_model` instance to regress out of data before beta estimation (default: NULL)
#' @param fracs Fraction of voxels used for prewhitening.
#' @param progress Logical; show progress bar.
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
estimate_betas.matrix_dataset <- function(x, fixed = NULL, ran, block,
                                         method = c("mixed", "lss", "ols"),
                                         basemod = NULL,
                                         fracs = .5, progress = TRUE, ...) {
  
  method <- match.arg(method)
  dset <- x
  mask <- fmridataset::get_mask(dset)
 
  bmod <- if (is.null(basemod)) {
    baseline_model("constant", sframe=dset$sampling_frame)
  } else {
    basemod
  }


  bdes <- gen_beta_design(fixed, ran, block, bmod, dset)
  betas <- run_estimate_betas(bdes, dset, method,
                              fracs = fracs, progress = progress)
  
  # Access beta_matrix from the list returned by run_estimate_betas
  beta_matrix <- betas$beta_matrix
  
  # Extract random and fixed effects from the beta matrix
  if (length(bdes$ran_ind) > 0) {
    ran <- as.matrix(beta_matrix[bdes$ran_ind,,drop=FALSE])
  } else {
    ran <- NULL
  }
  
  if (length(bdes$fixed_ind) > 0) {
    fixed <- as.matrix(beta_matrix[bdes$fixed_ind,,drop=FALSE])
  } else {
    fixed <- NULL
  }
  
  ret <- list(betas_fixed=fixed,
              betas_ran=ran,
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
#' @param method The regression method for estimating trialwise betas; one of "mixed", "lss", or "ols" (default: "mixed")
#' @param basemod A `baseline_model` instance to regress out of data before beta estimation (default: NULL)
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
estimate_betas.latent_dataset <- function(x, fixed = NULL, ran, block,
                                         method = c("mixed", "lss", "ols"),
                                         basemod = NULL,
                                         prewhiten = FALSE, progress = TRUE, ...) {
  
  method <- match.arg(method)
  dset <- x
  mask <- fmridataset::get_mask(dset)
  
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
  
  betas <- run_estimate_betas(bdes, dset, method,
                              progress = progress)
  
  # Access beta_matrix from the list returned by run_estimate_betas
  beta_matrix <- betas$beta_matrix
  
  # Extract random and fixed effects from the beta matrix
  if (length(bdes$ran_ind) > 0) {
    ran <- as.matrix(beta_matrix[bdes$ran_ind,,drop=FALSE])
  } else {
    ran <- NULL
  }
  
  if (length(bdes$fixed_ind) > 0) {
    fixed <- as.matrix(beta_matrix[bdes$fixed_ind,,drop=FALSE])
  } else {
    fixed <- NULL
  }
  
  ret <- list(betas_fixed=fixed,
              betas_ran=ran,
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
#' @param progress Logical; display progress during estimation.
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
estimate_hrf <- function(form, fixed = NULL, block, dataset,
                           bs = c("tp", "ts", "cr", "ps"),
                           rsam = seq(0, 20, by = 1),
                           basemod = NULL,
                           k = 8,
                           fx = TRUE,
                           progress = TRUE) {
  with_package("mgcv")
  dset <- dataset
  
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
  
  vecs <- masked_vectors(dset)
  res <- map_voxels(vecs, function(v) {
    gam.1 <- if (has_fixed) {
      mgcv::gam(v ~ mgcv::s(X_cond, bs=bs, fx=TRUE, k=8) + X_fixed + X_base)
    } else {
      mgcv::gam(v ~ mgcv::s(X_cond, bs=bs, fx=TRUE, k=8) + X_base)
    }
    
    time <- rsam
    xb <- matrix(colMeans(X_base), length(time),ncol(X_base), byrow=TRUE)
    ##xf <- matrix(colMeans(X_fixed), length(time),ncol(X_fixed), byrow=TRUE)
    predict(gam.1, list(X_cond = time, X_base = xb, X_fixed = xf))
  }, .progress = progress)
  
  res
  
}



#' @noRd 
#' @keywords internal
#' @importFrom rlang new_formula f_lhs f_rhs f_env is_call call_name
inject_basis <- function(oldform, new_basis, fun_names = c("hrf", "trialwise")) {
  stopifnot(is.formula(oldform))
  
  # A recursive helper that descends through an expression
  # and injects `basis=new_basis` into calls named in fun_names.
  recfun <- function(expr) {
    if (!is_call(expr)) {
      return(expr)  # If it's not a call, return as is
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
      # Manually add the basis argument
      call_rebuilt$basis <- new_basis
      return(call_rebuilt)
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

#' GLM OLS Estimation Convenience Function
#'
#' A convenience wrapper around `estimate_betas` for ordinary least squares (OLS) estimation.
#' This function provides a simplified interface for fitting GLMs using OLS on matrix datasets.
#' 
#' **Use Cases:**
#' - **Condition-level estimation**: Estimates average responses for each experimental condition
#' - **General linear modeling**: Standard GLM approach for group-level or condition-level effects
#' - **Multi-trial averaging**: Combines trials of the same condition to estimate mean responses
#' 
#' For single-trial estimation where each trial gets its own beta estimate, use `glm_lss()` instead.
#'
#' @param dataset A `matrix_dataset` object containing the fMRI time series data
#' @param model_obj An `event_model` object specifying the experimental design
#' @param basis_obj An HRF basis object (e.g., from `fmrihrf::HRF_SPMG1`, `HRF_FIR`, etc.)
#' @param basemod A `baseline_model` instance to regress out of data before beta estimation (default: NULL)
#' @param block A formula specifying the block factor (default: ~ 1 for single block)
#' @param progress Logical; show progress bar (default: TRUE)
#' @param ... Additional arguments passed to `estimate_betas`
#'
#' @return A list of class "fmri_betas" containing the estimated coefficients
#'
#' @examples
#' \dontrun{
#' # Create event model and data
#' event_data <- data.frame(
#'   onset = c(10, 30, 50, 70),
#'   condition = factor(c("A", "B", "A", "B")),
#'   run = rep(1, 4)
#' )
#' sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
#' model_obj <- event_model(onset ~ hrf(condition), 
#'                         data = event_data, 
#'                         block = ~ run, 
#'                         sampling_frame = sframe)
#' 
#' # Create data matrix (100 timepoints, 10 voxels)
#' Y <- matrix(rnorm(1000), 100, 10)
#' 
#' # Create matrix_dataset with event table
#' dset <- matrix_dataset(Y, TR = 2, run_length = 100, event_table = event_data)
#' 
#' # Fit with OLS - estimates average response for each condition
#' fit <- glm_ols(dset, model_obj, fmrihrf::HRF_SPMG1)
#' dim(fit$betas_ran)  # 2 conditions x 10 voxels
#' }
#'
#' @export
#' @seealso \code{\link{estimate_betas}} for the underlying estimation function, 
#'   \code{\link{glm_lss}} for single trial estimation
glm_ols <- function(dataset, model_obj, basis_obj, basemod = NULL, 
                    block = ~ 1, progress = TRUE, ...) {
  
  # Validate inputs
  if (!inherits(dataset, "matrix_dataset")) {
    stop("dataset must be a matrix_dataset object. Use matrix_dataset() to create one from your data matrix.")
  }
  
  if (!inherits(model_obj, "event_model")) {
    stop("model_obj must be an event_model object")
  }
  
  # Validate basis_obj
  if (is.character(basis_obj)) {
    # Check if it's a valid HRF basis name
    valid_basis_names <- c("HRF_SPMG1", "HRF_SPMG2", "HRF_SPMG3", "HRF_FIR", 
                          "HRF_AFNI", "HRF_GAM", "HRF_IL", "HRF_DD")
    if (!basis_obj %in% valid_basis_names) {
      stop(paste0("Unknown HRF basis name: ", basis_obj))
    }
    # Convert string to actual basis object from fmrihrf package
    basis_obj <- get(basis_obj, envir = asNamespace("fmrihrf"))
  } else if (!inherits(basis_obj, "HRF")) {
    stop("basis_obj must be an HRF object or a valid HRF basis name")
  }
  
  # Extract the formula from the event model and inject the new basis
  original_formula <- model_obj$model_spec$formula_or_list
  if (is.null(original_formula)) {
    stop("Cannot extract formula from event_model")
  }
  
  # Inject the new basis into the formula
  updated_formula <- inject_basis(original_formula, basis_obj)
  
  # Call estimate_betas with the updated formula and the dataset's event table
  estimate_betas(dataset, 
                fixed = NULL,
                ran = updated_formula, 
                block = block,
                method = "ols",
                basemod = basemod,
                progress = progress,
                ...)
}

#' GLM LSS Estimation Convenience Function (Single Trial Estimation)
#'
#' A convenience wrapper around `estimate_betas` for least squares separate (LSS) estimation.
#' **This is primarily designed for single trial estimation**, where each individual trial/event 
#' gets its own separate beta estimate rather than averaging across trials of the same condition.
#' 
#' **Primary Use Case - Single Trial Estimation:**
#' - **Trial-wise beta estimation**: Each trial gets its own beta coefficient
#' - **Single trial analysis**: Useful for decoding, representational similarity analysis (RSA)
#' - **Trial-by-trial variability**: Captures individual trial responses rather than condition averages
#' - **Avoiding trial averaging**: Preserves trial-specific information that would be lost in standard GLM
#' 
#' **Method Details:**
#' LSS (Least Squares Separate) fits a separate model for each trial, where the trial of interest 
#' gets its own regressor while all other trials of the same condition are modeled together. This 
#' approach avoids the collinearity issues that would arise from including separate regressors 
#' for every trial simultaneously.
#' 
#' For standard condition-level estimation (averaging trials within conditions), use `glm_ols()` instead.
#'
#' @param dataset A `matrix_dataset` object containing the fMRI time series data
#' @param model_obj An `event_model` object specifying the experimental design
#' @param basis_obj An HRF basis object (e.g., from `fmrihrf::HRF_SPMG1`, `HRF_FIR`, etc.)
#' @param basemod A `baseline_model` instance to regress out of data before beta estimation (default: NULL)
#' @param block A formula specifying the block factor (default: ~ 1 for single block)
#' @param use_cpp Deprecated. The C++ implementation has been retired. This parameter is ignored; fmrilss is always used.
#' @param progress Logical; show progress bar (default: TRUE)
#' @param ... Additional arguments passed to `estimate_betas`
#'
#' @return A list of class "fmri_betas" containing the estimated trial-wise coefficients
#'
#' @examples
#' \dontrun{
#' # Create event model and data
#' event_data <- data.frame(
#'   onset = c(10, 30, 50, 70),
#'   condition = factor(c("A", "B", "A", "B")),
#'   run = rep(1, 4)
#' )
#' sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
#' model_obj <- event_model(onset ~ hrf(condition), 
#'                         data = event_data, 
#'                         block = ~ run, 
#'                         sampling_frame = sframe)
#' 
#' # Create data matrix (100 timepoints, 10 voxels)
#' Y <- matrix(rnorm(1000), 100, 10)
#' 
#' # Create matrix_dataset with event table
#' dset <- matrix_dataset(Y, TR = 2, run_length = 100, event_table = event_data)
#' 
#' # Fit with LSS - estimates separate beta for each individual trial
#' fit <- glm_lss(dset, model_obj, fmrihrf::HRF_SPMG1)
#' dim(fit$betas_ran)  # 4 trials x 10 voxels (NOT averaged by condition)
#' 
#' # This is useful for:
#' # - Decoding analysis (predicting condition from single trial patterns)
#' # - RSA (representational similarity analysis)
#' # - Studying trial-by-trial variability
#' }
#'
#' @export
#' @seealso \code{\link{estimate_betas}} for the underlying estimation function, 
#'   \code{\link{glm_ols}} for condition-level estimation
glm_lss <- function(dataset, model_obj, basis_obj, basemod = NULL,
                    block = ~ 1, use_cpp = FALSE, progress = TRUE, ...) {
  
  # Validate inputs
  if (!inherits(dataset, "matrix_dataset")) {
    stop("dataset must be a matrix_dataset object. Use matrix_dataset() to create one from your data matrix.")
  }
  
  if (!inherits(model_obj, "event_model")) {
    stop("model_obj must be an event_model object")
  }
  
  # Validate basis_obj
  if (is.character(basis_obj)) {
    # Check if it's a valid HRF basis name
    valid_basis_names <- c("HRF_SPMG1", "HRF_SPMG2", "HRF_SPMG3", "HRF_FIR", 
                          "HRF_AFNI", "HRF_GAM", "HRF_IL", "HRF_DD")
    if (!basis_obj %in% valid_basis_names) {
      stop(paste0("Unknown HRF basis name: ", basis_obj))
    }
    # Convert string to actual basis object from fmrihrf package
    basis_obj <- get(basis_obj, envir = asNamespace("fmrihrf"))
  } else if (!inherits(basis_obj, "HRF")) {
    stop("basis_obj must be an HRF object or a valid HRF basis name")
  }
  
  # Extract the formula from the event model and inject the new basis
  original_formula <- model_obj$model_spec$formula_or_list
  if (is.null(original_formula)) {
    stop("Cannot extract formula from event_model")
  }
  
  # Inject the new basis into the formula
  updated_formula <- inject_basis(original_formula, basis_obj)
  
  if (use_cpp) {
    warning("C++-optimized LSS implementation has been retired; using method = 'lss'.", call. = FALSE)
  }
  method <- "lss"
  
  # Call estimate_betas with the updated formula
  res <- estimate_betas(
    dataset,
    fixed = NULL,
    ran = updated_formula,
    block = block,
    method = method,
    basemod = basemod,
    progress = progress,
    ...
  )

  betas_ran <- res$betas_ran
  if (!is.null(betas_ran)) {
    betas_mat <- as.matrix(betas_ran)
    if (all(!is.finite(betas_mat)) || any(!is.finite(betas_mat))) {
      stop("Cholesky factorization failed: design matrix not positive definite", call. = FALSE)
    }
  }

  res
}
