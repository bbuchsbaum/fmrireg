#' Estimate Rank-1 HRF Using LS+SVD/CF-ALS Methods
#'
#' High level wrapper that prepares design matrices and dispatches to
#' the desired estimation engine.  This function supports the fast
#' \emph{LS+SVD} initialisation, the one-step refinement
#' \emph{LS+SVD+1ALS}, or the full alternating scheme implemented in
#' `cf_als_engine`.
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
#' @param method Estimation method. One of "ls_svd_only",
#'   "ls_svd_1als" (default) or "cf_als".
#' @param lambda_init Ridge penalty for the initial GLM solve used by
#'   `ls_svd` based methods.
#' @param lambda_b Ridge penalty for the beta update step.
#' @param lambda_h Ridge penalty for the h update step.
#' @param fullXtX Logical; if `TRUE` include cross condition terms in
#'   the h update (where supported).
#' @param max_alt Number of alternating updates after initialisation
#'   when `method = "cf_als"`.
#' @return An object of class `fmrireg_cfals_fit` containing the
#'   estimated HRF coefficients and amplitudes.
#' @details
#' The `method` argument selects between the closed-form
#' \code{"ls_svd_only"}, the default \code{"ls_svd_1als"} which adds one
#' ALS refinement step, or the iterative \code{"cf_als"} engine.  The
#' ridge penalties \code{lambda_init}, \code{lambda_b} and
#' \code{lambda_h} control regularisation of the initial solve, the
#' beta-update and the h-update respectively.  Setting
#' \code{fullXtX = TRUE} includes cross-condition terms in the h-update
#' (when supported by the chosen engine).  R\eqn{^2} is computed on the
#' data after confound projection.
#'
#' @examples
#' sframe <- sampling_frame(blocklens = 40, TR = 1)
#' ev_df <- data.frame(onset = c(5, 15, 25), block = 1)
#' emod <- event_model(onset ~ hrf(constant), data = ev_df,
#'                     block = ~ block, sampling_frame = sframe)
#' Y <- matrix(rnorm(40 * 2), 40, 2)
#' fit <- fmrireg_cfals(Y, emod, HRF_SPMG1)
#' print(fit)
#' @export
fmrireg_cfals <- function(fmri_data_obj,
                         event_model,
                         hrf_basis,
                         confound_obj = NULL,
                         method = c("ls_svd_1als", "ls_svd_only", "cf_als"),
                         lambda_init = 1,
                         lambda_b = 10,
                         lambda_h = 1,
                         fullXtX = FALSE,
                         max_alt = 1) {

  method <- match.arg(method)

  if (inherits(fmri_data_obj, "fmri_dataset")) {
    Y <- get_data_matrix(fmri_data_obj)
  } else if (is.matrix(fmri_data_obj)) {
    Y <- fmri_data_obj
  } else {
    stop("'fmri_data_obj' must be an 'fmri_dataset' or matrix")
  }

  sframe <- if (inherits(fmri_data_obj, "fmri_dataset"))
    fmri_data_obj$sampling_frame else attr(fmri_data_obj, "sampling_frame")

  if (is.null(sframe)) {
    stop("Sampling information could not be determined from input")
  }

  design <- create_fmri_design(event_model, hrf_basis)
  X_list <- design$X_list
  cond_names <- names(X_list)

  proj <- project_confounds(Y, X_list, confound_obj)
  Xp <- proj$X_list
  Yp <- proj$Y

  fit <- switch(method,
    ls_svd_only = ls_svd_engine(Xp, Yp,
                                lambda_init = lambda_init,
                                h_ref_shape_norm = design$h_ref_shape_norm),
    ls_svd_1als = ls_svd_1als_engine(Xp, Yp,
                                     lambda_init = lambda_init,
                                     lambda_b = lambda_b,
                                     lambda_h = lambda_h,
                                     fullXtX_flag = fullXtX,
                                     h_ref_shape_norm = design$h_ref_shape_norm),
    cf_als = cf_als_engine(Xp, Yp,
                           lambda_b = lambda_b,
                           lambda_h = lambda_h,
                           fullXtX_flag = fullXtX,
                           h_ref_shape_norm = design$h_ref_shape_norm,
                           max_alt = max_alt)
  )

  rownames(fit$beta) <- cond_names

  Phi <- design$Phi
  recon_hrf <- Phi %*% fit$h

  n <- nrow(Yp)
  v <- ncol(Yp)
  pred_p <- Reduce(`+`, Map(function(Xc, bc) {
    Xc %*% (fit$h * matrix(bc, nrow = nrow(fit$h), ncol = v, byrow = TRUE))
  }, Xp, asplit(fit$beta, 1)))
  resids <- Yp - pred_p

  SST <- colSums((Yp - matrix(colMeans(Yp), n, v, TRUE))^2)
  SSE <- colSums(resids^2)
  r2 <- 1 - SSE / SST

  out <- fmrireg_cfals_fit(h_coeffs = fit$h,
                           beta_amps = fit$beta,
                           method = method,
                           lambdas = c(init = lambda_init,
                                       beta = lambda_b,
                                       h = lambda_h),
                           call = match.call(),
                           hrf_basis = hrf_basis,
                           design_info = list(d = design$d,
                                              k = length(X_list),
                                              n = n,
                                              v = v,
                                              fullXtX = fullXtX),
                           residuals = resids,
                           recon_hrf = recon_hrf,
                           gof = r2)
  out
}

#' Fit Rank-1 HRF Using CF-ALS
#'
#' Convenience wrapper for the original CF-ALS implementation.  This
#' function simply calls [fmrireg_cfals()] with `method = "cf_als"` and
#' retains the historical argument names.
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
#' @details R\eqn{^2} is computed on the data after confound projection.
#' @export

fmrireg_hrf_cfals <- function(fmri_data_obj,
                              event_model,
                              hrf_basis,
                              confound_obj = NULL,
                              lam_beta = 10,
                              lam_h = 1,
                              fullXtX = FALSE,
                              max_alt = 1) {
  fmrireg_cfals(fmri_data_obj,
                event_model,
                hrf_basis,
                confound_obj = confound_obj,
                method = "cf_als",
                lambda_b = lam_beta,
                lambda_h = lam_h,
                fullXtX = fullXtX,
                max_alt = max_alt)
}
