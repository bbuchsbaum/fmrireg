#' Robust Iterative Fitting Engine
#'
#' Performs Iteratively Reweighted Least Squares (IRLS) for robust regression
#' using either Huber or Tukey bisquare weights. This function takes an initial
#' GLM context (typically from OLS on original or whitened data) and iteratively
#' refines the fit using robust weights.
#'
#' @param initial_glm_ctx A \code{glm_context} object containing initial fit results
#' @param cfg_robust_options List of robust fitting options from \code{fmri_lm_control$robust}
#' @param X_orig_for_resid The design matrix corresponding to Y in initial_glm_ctx
#'   before any robust weighting (needed for residual calculation)
#' @param sigma_fixed Optional fixed sigma value (for global scale estimation)
#' @param obs_weights Optional prior observation weights already applied to the
#'   rows of `X`/`Y` (e.g. volume weights). Used only to count the
#'   observations that take part in the fit: a row with zero prior weight is
#'   not an observation.
#' @param weights_fixed Optional robust weights estimated beforehand (one per
#'   observation). The IRLS loop is skipped and the result is the WLS fit with
#'   these weights -- the same final solve the loop would end with had it
#'   produced them. Used to fit voxel chunks with weights estimated once over
#'   all voxels (`.integrated_noise_estimate()`); the MAD scale fields are
#'   then `NULL`.
#'
#' @return A list containing:
#' \describe{
#'   \item{betas_robust}{Final robust coefficient estimates: exactly the
#'     weighted least-squares fit with `robust_weights_final`.}
#'   \item{XtWXi_final}{(X'WX)^-1 for those weights}
#'   \item{sigma2}{Residual variance of that WLS fit, per response column:
#'     \eqn{\sum_i w_i r_i^2 / (n_{w>0} - \mathrm{rank})}, as
#'     `stats::lm(weights = w)` reports it. This is the scale used for
#'     inference by every robust path.}
#'   \item{rss}{Weighted residual sum of squares \eqn{\sum_i w_i r_i^2}.}
#'   \item{sigma_robust_scale_final}{The MAD scale that produced the final
#'     weights (diagnostic; not used for standard errors).}
#'   \item{robust_weights_final}{The weights used for the final solve.}
#'   \item{dfres}{Residual df of the WLS fit: observations with non-zero
#'     weight minus the rank of the weighted design.}
#' }
#'
#' @details
#' The IRLS algorithm:
#' 1. Calculate residuals using X_orig_for_resid and current betas
#' 2. Estimate robust scale (MAD) unless sigma_fixed is provided
#' 3. Calculate robust weights based on scaled residuals
#' 4. Create weighted GLM context and solve
#' 5. Iterate until convergence or max_iter
#'
#' The returned weights, betas, \eqn{(X'WX)^{-1}} and scale all belong to the
#' same (last) weighted solve, so the fit is reproducible as a WLS fit.
#'
#' @keywords internal
#' @importFrom matrixStats rowMedians
#' @noRd
robust_iterative_fitter <- function(initial_glm_ctx,
                                    cfg_robust_options,
                                    X_orig_for_resid,
                                    sigma_fixed = NULL,
                                    obs_weights = NULL,
                                    weights_fixed = NULL) {
  
  # Validate inputs
  stopifnot(is.glm_context(initial_glm_ctx))
  stopifnot(is.list(cfg_robust_options))
  stopifnot(!is.null(X_orig_for_resid))
  
  # Extract options
  psi_type <- cfg_robust_options$type
  if (isFALSE(psi_type)) {
    stop("robust_iterative_fitter called with robust type = FALSE")
  }
  
  k_huber <- cfg_robust_options$k_huber
  c_tukey <- cfg_robust_options$c_tukey
  max_iter <- max(1L, as.integer(cfg_robust_options$max_iter %||% 2L))
  tol <- as.numeric(cfg_robust_options$tol %||% 1e-4)
  min_iter <- as.integer(cfg_robust_options$min_iter %||% 1L)
  if (!is.finite(tol) || tol < 0) tol <- 1e-4
  if (!is.finite(min_iter) || min_iter < 1L) min_iter <- 1L
  scope_opt <- cfg_robust_options$scale_scope
  if (is.null(scope_opt)) scope_opt <- "run"
  scale_scope <- tolower(as.character(scope_opt))
  if (identical(scale_scope, "local")) scale_scope <- "voxel"
  if (!(scale_scope %in% c("run", "global", "voxel"))) {
    scale_scope <- "run"
  }
  
  # Ensure matrices
  if (!is.matrix(X_orig_for_resid)) X_orig_for_resid <- as.matrix(X_orig_for_resid)
  Y <- initial_glm_ctx$Y
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  
  # Check for NAs
  if (anyNA(Y) || anyNA(X_orig_for_resid)) {
    stop("NA values detected in 'X_orig_for_resid' or 'Y' for robust_iterative_fitter")
  }
  
  if (!is.null(weights_fixed)) {
    weights_fixed <- as.numeric(weights_fixed)
    if (length(weights_fixed) != nrow(Y) || anyNA(weights_fixed)) {
      stop("weights_fixed must hold one finite weight per observation")
    }
  } else {
    # Initial fit from context (could be OLS or whitened OLS)
    betas <- solve_glm_core(initial_glm_ctx, return_fitted = FALSE)$betas
  }

  # Initialize final values
  XtWXi_final <- initial_glm_ctx$proj$XtXinv
  dfres <- initial_glm_ctx$proj$dfres

  # Robust scale and per-observation scaled residual u (row median over
  # response columns of |r| / sigma). Takes |r| directly so the caller can
  # build it without intermediate copies; the division is by a per-column
  # vector expanded once (sweep() allocated three T x V blocks).
  compute_scale <- function(abs_resid, scope, sigma_fixed_value = NULL) {
    eps <- .Machine$double.eps

    if (!is.null(sigma_fixed_value)) {
      sigma_vec <- as.numeric(sigma_fixed_value)
      if (length(sigma_vec) == 1L) {
        sigma_vec <- rep(sigma_vec, ncol(abs_resid))
      } else if (length(sigma_vec) != ncol(abs_resid)) {
        stop("sigma_fixed must have length 1 or match number of response columns")
      }
      sigma_vec[!is.finite(sigma_vec) | sigma_vec <= eps] <- eps
    } else if (identical(scope, "voxel")) {
      sigma_vec <- 1.4826 * matrixStats::colMedians(abs_resid)
      sigma_vec[!is.finite(sigma_vec) | sigma_vec <= eps] <- eps
    } else {
      row_med <- matrixStats::rowMedians(abs_resid)
      sigma_scalar <- 1.4826 * stats::median(row_med)
      if (!is.finite(sigma_scalar) || sigma_scalar <= eps) {
        sigma_scalar <- eps
      }
      sigma_vec <- rep(sigma_scalar, ncol(abs_resid))
    }

    u <- if (length(unique(sigma_vec)) == 1L) {
      matrixStats::rowMedians(abs_resid / sigma_vec[1L])
    } else {
      matrixStats::rowMedians(abs_resid / rep(sigma_vec, each = nrow(abs_resid)))
    }
    sigma_scalar <- stats::median(sigma_vec)
    if (!is.finite(sigma_scalar) || sigma_scalar <= eps) {
      sigma_scalar <- eps
    }

    list(u = u, sigma_scalar = sigma_scalar, sigma_components = sigma_vec)
  }

  psi_weights <- function(u) {
    switch(psi_type,
           huber = pmin(1, k_huber / abs(u)),
           bisquare = ifelse(abs(u) <= c_tukey, (1 - (u / c_tukey)^2)^2, 0))
  }

  n_obs <- nrow(Y)
  active_prior <- if (is.null(obs_weights)) {
    rep(TRUE, n_obs)
  } else {
    as.numeric(obs_weights) > 0
  }

  # IRLS iterations. Everything returned below comes from the last weighted
  # solve: the weights it used, its betas, its (X'WX)^-1 and its residuals.
  #
  # The weights are one per observation (shared by all response columns), so
  # the WLS solve needs only the p x n pseudo-inverse of Xw = sqrt(W) X:
  # betas = Pinv(Xw) sqrt(W) Y = (Pinv(Xw) diag(sqrt(w))) Y. Scaling the
  # columns of that p x n matrix replaces the T x V weighted response
  # sqrt(W) Y and its QR projection.
  w <- rep(1, n_obs)
  scale_info <- NULL
  proj_w <- initial_glm_ctx$proj
  wls_solve <- function(w) {
    sqrtw <- sqrt(w)
    Xw <- X_orig_for_resid * sqrtw
    proj <- .fast_preproject(Xw)
    list(proj = proj, betas = .robust_wls_betas(proj, sqrtw, Y, colnames(Xw)))
  }
  if (!is.null(weights_fixed)) {
    # Weights estimated beforehand (e.g. once over all voxels): the fit is
    # the WLS solve the IRLS loop would end with for these weights.
    w <- weights_fixed
    sol <- wls_solve(w)
    proj_w <- sol$proj
    betas <- sol$betas
    XtWXi_final <- proj_w$XtXinv
  }
  for (it in seq_len(if (is.null(weights_fixed)) max_iter else 0L)) {
    betas_prev <- betas

    # |residuals| using the original (unweighted) X, in one expression so R
    # reuses the temporaries instead of allocating new T x V blocks.
    scale_info <- compute_scale(abs(Y - X_orig_for_resid %*% betas),
                                scale_scope, sigma_fixed)
    w <- psi_weights(scale_info$u)

    sol <- wls_solve(w)
    proj_w <- sol$proj
    betas <- sol$betas
    XtWXi_final <- proj_w$XtXinv

    # Convergence check on coefficient updates.
    if (it >= min_iter) {
      scale <- max(1, max(abs(betas_prev)))
      delta_beta <- max(abs(betas - betas_prev)) / scale
      if (is.finite(delta_beta) && delta_beta <= tol) {
        break
      }
    }
  }

  # Weighted RSS of the final solve, sum_i w_i r_i^2 per column (max_iter is
  # at least 1, so a weighted solve always took place).
  rss_w <- colSums(w * (Y - X_orig_for_resid %*% betas)^2)

  # Residual df of the WLS fit: a zero-weight row carries no information, as
  # in stats::lm.wfit().
  dfres <- sum(w > 0 & active_prior) - proj_w$rank
  sigma2 <- if (dfres > 0) rss_w / dfres else rep(NA_real_, ncol(Y))

  list(
    betas_robust = betas,
    XtWXi_final = XtWXi_final,
    sigma2 = sigma2,
    rss = rss_w,
    sigma_robust_scale_final = scale_info$sigma_scalar,
    sigma_robust_scale_components = scale_info$sigma_components,
    robust_weights_final = w,
    dfres = dfres,
    rank = proj_w$rank,
    scale_scope = scale_scope
  )
}

#' WLS coefficients for observation weights shared by all response columns
#'
#' With `Xw = sqrt(W) X` and `proj_w = .fast_preproject(Xw)`, returns
#' `Pinv(Xw) sqrt(W) Y` -- the coefficients `solve_glm_core()` gives for the
#' weighted response `sqrt(W) Y` -- without forming that n x V matrix: the
#' square-root weights scale the columns of the p x n pseudo-inverse instead.
#' Row names and rank attributes follow `solve_glm_core()`.
#' @keywords internal
#' @noRd
.robust_wls_betas <- function(proj_w, sqrtw, Y, varnames = NULL) {
  Pw <- .preproject_pinv(proj_w)
  betas <- (Pw * rep(sqrtw, each = nrow(Pw))) %*% Y
  attributes(betas) <- list(dim = dim(betas))
  aliased <- as.integer(proj_w$aliased %||% integer(0))
  if (length(aliased) > 0L) betas[aliased, ] <- 0
  if (!is.null(varnames) && nrow(betas) == length(varnames)) {
    rownames(betas) <- varnames
  }
  if (length(aliased) > 0L) {
    attr(betas, "rank_deficient") <- TRUE
    attr(betas, "aliased") <- aliased
    attr(betas, "estimable") <- as.integer(proj_w$estimable %||% integer(0))
  }
  betas
}
