#' @keywords internal
#' @noRd
ols_betas <- function(X, Y) {
  fit <- lm.fit(as.matrix(X),Y)
  coef(fit)
}


#' @keywords internal
#' @noRd
mixed_betas <- function(X, Y, ran_ind, fixed_ind, solver = NULL) {
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
  
  all_indices <- c(ran_ind, fixed_ind)
  if (anyNA(all_indices) || any(all_indices < 1L) ||
      any(all_indices != as.integer(all_indices))) {
    stop("Random and fixed effect indices must be positive integers", call. = FALSE)
  }

  # Check if indices are out of bounds
  if (max(all_indices) > ncol(X)) {
    stop("Index out of bounds: indices exceed number of columns in X")
  }
  
  # Ensure Y is a vector
  if (!is.vector(Y)) {
    Y <- as.vector(Y)
  }

  if (length(Y) != nrow(X)) {
    stop("Y must have the same number of rows as X", call. = FALSE)
  }
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("X and Y must contain only finite values", call. = FALSE)
  }

  failure_value <- function() {
    rep(NA_real_, length(ran_ind) + length(fixed_ind))
  }

  if (!is.null(solver) && !is.function(solver)) {
    stop("solver must be NULL or a function", call. = FALSE)
  }
  
  # Try C++ mixed model implementation
  tryCatch({
    
    if (!is.null(solver) ||
        (requireNamespace("Rcpp", quietly = TRUE) &&
         requireNamespace("fmrilss", quietly = TRUE))) {
      
      # Use the C++ implementation with proper input validation
      X_fixed <- if (length(fixed_ind) == 0) {
        matrix(1, nrow = nrow(X), ncol = 1)
      } else {
        X[, fixed_ind, drop = FALSE]
      }
      
      solve_fun <- if (is.null(solver)) fmrilss::mixed_solve else solver
      fit <- tryCatch({
        solve_fun(Y = Y,
                  Z = X[, ran_ind, drop = FALSE],
                  X = X_fixed)
      }, error = function(e2) {
        warning(
          "Mixed model solver failed: ", conditionMessage(e2),
          "; returning NA coefficients for this response.",
          call. = FALSE
        )
        NULL
      })

      if (is.null(fit)) return(failure_value())
      
      # Return results based on whether fixed_ind is empty
      if (length(fixed_ind) == 0) {
        return(fit$u)
      } else {
        return(c(fit$u, fit$beta))
      }
    } else {
      warning(
        "No mixed model solver is available; returning NA coefficients for this response.",
        call. = FALSE
      )
      return(failure_value())
    }
  })
}

#' Estimate betas using various regression methods
#'
#' This function estimates betas (regression coefficients) for fixed and random effects
#' using various regression methods including mixed models, least squares, and PLS.
#'
#' @param x An `fmri_frame` (see [matrix_frame()], [neurovec_frame()],
#'   [nifti_frame()], and [latent_frame()]). When the frame's feature space is a
#'   `volume_space`, the fixed and random betas are returned as `NeuroVec`
#'   objects on that grid; otherwise (index or basis spaces) they are returned
#'   as coefficient-by-feature matrices.
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
#'   * betas_fixed: fixed effect betas (NeuroVec for volumetric frames, matrix otherwise).
#'   * betas_ran: random effect betas (NeuroVec for volumetric frames, matrix otherwise).
#'   * design_ran: Design matrix for random effects.
#'   * design_fixed: Design matrix for fixed effects.
#'   * design_base: Design matrix for baseline model.
#'   * basemod: Baseline model object.
#'   * fixed_model: Fixed effect model object.
#'   * ran_model: Random effect model object.
#'   * estimated_hrf: The estimated HRF vector (NULL for most methods).
#'
#' @seealso \code{\link{matrix_frame}}, \code{\link{baseline_model}}, \code{\link{event_model}}
#'
#' @examples
#' \dontrun{
#' facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
#' facedes$frun <- factor(facedes$run)
#' scans <- paste0("rscan0", 1:6, ".nii")
#'
#' dset <- nifti_frame(scans=scans, mask="mask.nii", TR=1.5,
#'         run_length=rep(436,6), event_table=facedes)
#' fixed = onset ~ hrf(run)
#' ran = onset ~ trialwise()
#' block = ~ run
#'
#' betas <- estimate_betas(dset, fixed=fixed, ran=ran, block=block, method="mixed")
#' }
#' @family estimate_betas
#' @rdname estimate_betas
#' @export
estimate_betas.fmri_frame <- function(x, fixed = NULL, ran, block,
                                      method = c("mixed", "lss", "ols"),
                                      basemod = NULL,
                                      maxit = 1000,
                                      fracs = 0.5,
                                      progress = TRUE,
                                      ...) {
  method <- match.arg(method)
  dset <- x
  volumetric <- inherits(.dset_space(dset), "volume_space")
  spatial <- if (volumetric) .fmri_dataset_mask_space(dset, "beta image reconstruction") else NULL

  bmod <- if (is.null(basemod)) {
    baseline_model("constant", sframe = .dset_sampling_frame(dset))
  } else {
    basemod
  }

  bdes <- gen_beta_design(fixed, ran, block, bmod, dset, method = method)
  betas <- run_estimate_betas(bdes, dset, method, block = block,
                              maxit = maxit, fracs = fracs,
                              progress = progress, ...)
  beta_matrix <- betas$beta_matrix

  if (volumetric) {
    mask <- spatial$mask
    ospace_ran <- neuroim2::add_dim(spatial$space, length(bdes$ran_ind))
    if (!is.null(bdes$fixed_ind)) {
      ospace_fixed <- neuroim2::add_dim(spatial$space, length(bdes$fixed_ind))
      fixed <- neuroim2::NeuroVec(as.matrix(beta_matrix[bdes$fixed_ind, , drop = FALSE]), ospace_fixed, mask = mask)
    } else {
      fixed <- NULL
    }
    ran <- neuroim2::NeuroVec(as.matrix(beta_matrix[bdes$ran_ind, , drop = FALSE]), ospace_ran, mask = mask)
  } else {
    ran <- if (length(bdes$ran_ind) > 0) {
      as.matrix(beta_matrix[bdes$ran_ind, , drop = FALSE])
    } else {
      NULL
    }
    fixed <- if (length(bdes$fixed_ind) > 0) {
      as.matrix(beta_matrix[bdes$fixed_ind, , drop = FALSE])
    } else {
      NULL
    }
  }

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
    data_matrix <- .dset_data_matrix(dset)
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
      vecs <- neuroim2::vectors(data_matrix, subset = seq_len(ncol(data_matrix)))
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
    emod_fixed <- event_model(fixed, data = .dset_event_table(dset), block = block, sampling_frame = .dset_sampling_frame(dset))
    dmat_fixed <- design_matrix(emod_fixed)
  } else {
    emod_fixed <- NULL
    dmat_fixed <- NULL
  }
  
  emod_ran <- event_model(ran, data = .dset_event_table(dset), block = block, sampling_frame = .dset_sampling_frame(dset))
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

#' @noRd 
#' @keywords internal
#' @importFrom rlang new_formula f_lhs f_rhs f_env is_call call_name
inject_basis <- function(oldform, new_basis, fun_names = c("hrf", "trialwise", "feature")) {
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
      # Not hrf() / trialwise() / feature(), so keep walking
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
#' @param dataset An `fmri_frame` containing the fMRI time series data (see [matrix_frame()])
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
#' # Create an fmri_frame with event table
#' dset <- matrix_frame(Y, TR = 2, run_length = 100, event_table = event_data)
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
  if (!.is_fmri_frame(dataset)) {
    stop("dataset must be an fmri_frame. Use matrix_frame() to create one from your data matrix.")
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
#' @param dataset An `fmri_frame` containing the fMRI time series data (see [matrix_frame()])
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
#' # Create an fmri_frame with event table
#' dset <- matrix_frame(Y, TR = 2, run_length = 100, event_table = event_data)
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
  if (!.is_fmri_frame(dataset)) {
    stop("dataset must be an fmri_frame. Use matrix_frame() to create one from your data matrix.")
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
