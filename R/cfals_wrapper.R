#' Fit Rank-1 HRF Using CF-ALS
#'
#' High level wrapper implementing CF-ALS HRF estimation
#' for data stored in standard `fmrireg` objects.
#'
#' @param fmri_data_obj An `fmri_dataset` or numeric matrix of BOLD data
#'   (time points x voxels). If a dataset, sampling information is
#'   taken from the object.
#' @param event_model An `event_model` describing the stimuli to use
#'   for HRF estimation.
#' @param hrf_basis An `HRF` basis object used for the convolution
#'   design matrices.
#' @param confound_obj Optional matrix of confound regressors with the
#'   same number of rows as the data matrix.
#' @param lam_beta Ridge penalty for the beta update step.
#' @param lam_h Ridge penalty for the h update step.
#' @param fullXtX Logical; if TRUE include cross condition terms in the
#'   h update.
#' @param max_alt Number of alternating updates after initialisation.
#' @return An object of class `fmrireg_cfals_fit` containing the
#'   estimated HRF coefficients and amplitudes.
#' @export
fmrireg_hrf_cfals <- function(fmri_data_obj,
                              event_model,
                              hrf_basis,
                              confound_obj = NULL,
                              lam_beta = 10,
                              lam_h = 1,
                              fullXtX = FALSE,
                              max_alt = 1) {

  if (inherits(fmri_data_obj, "fmri_dataset")) {
    Y <- get_data_matrix(fmri_data_obj)
    sframe <- fmri_data_obj$sampling_frame
  } else if (is.matrix(fmri_data_obj)) {
    Y <- fmri_data_obj
    sframe <- attr(fmri_data_obj, "sampling_frame")
    if (is.null(sframe)) {
      stop("Matrix input must have a 'sampling_frame' attribute")
    }
  } else {
    stop("'fmri_data_obj' must be an 'fmri_dataset' or matrix")
  }

  if (!inherits(event_model, "event_model")) {
    stop("'event_model' must be an 'event_model' object")
  }

  if (!inherits(hrf_basis, "HRF")) {
    stop("'hrf_basis' must be an object of class 'HRF'")
  }

  # build regressors for each condition using the provided basis
  reg_lists <- lapply(event_model$terms, regressors.event_term,
                      hrf = hrf_basis,
                      sampling_frame = sframe,
                      summate = FALSE,
                      drop.empty = TRUE)
  regs <- unlist(reg_lists, recursive = FALSE)
  cond_names <- names(regs)
  sample_times <- samples(sframe, global = TRUE)
  X_list <- lapply(regs, function(r) {
    evaluate(r, sample_times, precision = sframe$precision)
  })
  names(X_list) <- cond_names

  proj <- project_confounds(X_list, Y, confound_obj)
  Xp <- proj$X_list
  Yp <- proj$Y

  fit <- cf_als_engine(Xp, Yp,
                       lambda_b = lam_beta,
                       lambda_h = lam_h,
                       fullXtX_flag = fullXtX,
                       max_alt = max_alt)

  # reconstruct HRF shapes on the sampling grid
  Phi <- reconstruction_matrix(hrf_basis, sframe)
  recon_hrf <- Phi %*% fit$h

  # predicted BOLD and residuals in the projected space
  n <- nrow(Yp)
  v <- ncol(Yp)
  pred_p <- Reduce(`+`, Map(function(Zc, bc) {
    Zc %*% (fit$h * bc)
  }, Xp, asplit(fit$beta, 1)))
  resids <- Yp - pred_p

  # simple R^2 per voxel computed on projected data
  SST <- colSums((Yp - matrix(colMeans(Yp), n, v, TRUE))^2)
  SSE <- colSums(resids^2)
  r2 <- 1 - SSE / SST

  out <- list(h = fit$h,
              beta = fit$beta,
              reconstructed_hrfs = recon_hrf,
              residuals = resids,
              gof_per_voxel = r2,
              hrf_basis_used = hrf_basis,
              lambda_used = c(beta = lam_beta, h = lam_h),
              design_info = list(d = nbasis(hrf_basis),
                                 k = length(X_list),
                                 n = nrow(Y),
                                 v = ncol(Y),
                                 fullXtX = fullXtX))
  class(out) <- c("fmrireg_cfals_fit", "list")
  out
}

