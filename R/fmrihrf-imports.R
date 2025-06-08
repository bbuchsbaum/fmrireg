# Auto-generated file for fmrihrf imports

#' @import fmrihrf
#' @importFrom fmrihrf HRF as_hrf gen_hrf bind_basis
#' @importFrom fmrihrf hrf_spmg1 hrf_gamma hrf_gaussian
#' @importFrom fmrihrf hrf_bspline hrf_mexhat hrf_sine hrf_inv_logit hrf_half_cosine
#' @importFrom fmrihrf lag_hrf block_hrf normalise_hrf
#' @importFrom fmrihrf regressor sampling_frame
#' @importFrom fmrihrf evaluate nbasis penalty_matrix reconstruction_matrix
#' @importFrom fmrihrf blocklens blockids samples global_onsets
#' @importFrom fmrihrf HRF_SPMG1 HRF_SPMG2 HRF_SPMG3 HRF_GAMMA HRF_GAUSSIAN HRF_BSPLINE HRF_FIR
NULL

# Re-export commonly used functions for backward compatibility
# Note: We don't re-export generics that conflict with fmrireg's own generics
#' @export
fmrihrf::HRF

#' @export
fmrihrf::as_hrf

#' @export
fmrihrf::gen_hrf

#' @export
fmrihrf::sampling_frame

#' @export
fmrihrf::hrf_spmg1

#' @export
fmrihrf::regressor

#' @export
fmrihrf::blocklens

#' @export
fmrihrf::blockids

#' @export
fmrihrf::samples

#' @export
fmrihrf::global_onsets

#' @export
fmrihrf::evaluate

#' @export
fmrihrf::nbasis

