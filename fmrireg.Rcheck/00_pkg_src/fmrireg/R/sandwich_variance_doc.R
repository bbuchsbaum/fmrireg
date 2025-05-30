#' Sandwich Variance Estimation in fmrireg
#'
#' @description
#' This documentation describes the sandwich variance estimation capabilities
#' in fmrireg, which provide robust standard errors for regression coefficients
#' when model assumptions are violated.
#'
#' @section Background:
#' The sandwich variance estimator (also known as the Huber-White estimator)
#' provides valid standard errors even when the residuals exhibit 
#' heteroscedasticity or other violations of the classical linear model
#' assumptions.
#'
#' @section Mathematical Details:
#' The sandwich estimator is computed as:
#' \deqn{V_{sandwich} = (X'X)^{-1} X' \Omega X (X'X)^{-1}}
#' 
#' where \eqn{\Omega} is a diagonal matrix with squared residuals on the diagonal.
#' For robust regression with weights \eqn{w_i}, the weighted version is:
#' \deqn{V_{sandwich} = (X'WX)^{-1} X'W \Omega WX (X'WX)^{-1}}
#'
#' @section Usage in fmrireg:
#' Sandwich variance estimation is automatically used when:
#' - Robust regression is enabled (using M-estimators)
#' - AR modeling is combined with robust regression
#' - Heteroscedasticity is suspected in the residuals
#'
#' @section Effective Degrees of Freedom:
#' When using robust regression and/or AR models, the effective degrees of
#' freedom are adjusted to account for:
#' - Downweighting of outliers in robust regression
#' - Loss of degrees of freedom due to AR parameter estimation
#' 
#' The adjustment formula is:
#' \deqn{df_{effective} = df_{base} \times \frac{\sum w_i}{n} \times \frac{n - p_{AR}}{n}}
#'
#' @examples
#' \dontrun{
#' # Fit model with robust regression
#' cfg <- fmri_lm_control(
#'   robust = list(
#'     type = "bisquare",
#'     c_tukey = 4.685
#'   )
#' )
#' 
#' fit <- fmri_lm(model, dataset, config = cfg)
#' 
#' # Standard errors in fit$betas will use sandwich variance
#' # P-values will use effective degrees of freedom
#' }
#'
#' @section Implementation Notes:
#' - Small sample corrections are applied (n/(n-p) factor)
#' - For multi-voxel data, computation is vectorized for efficiency
#' - Compatible with all contrast types (t, F, custom)
#'
#' @references
#' Huber, P. J. (1967). The behavior of maximum likelihood estimates under
#' nonstandard conditions. Proceedings of the Fifth Berkeley Symposium on
#' Mathematical Statistics and Probability.
#'
#' White, H. (1980). A heteroskedasticity-consistent covariance matrix
#' estimator and a direct test for heteroskedasticity. Econometrica, 48(4),
#' 817-838.
#'
#' @name sandwich_variance
#' @keywords internal
NULL