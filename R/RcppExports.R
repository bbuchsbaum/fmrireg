# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

instantaneous_correlation_rcpp <- function(x, y, eta, offset = 0L) {
    .Call('_fmrireg_instantaneous_correlation_rcpp', PACKAGE = 'fmrireg', x, y, eta, offset)
}

compute_residuals_cpp <- function(X_base_fixed, data_matrix, dmat_ran) {
    .Call('_fmrireg_compute_residuals_cpp', PACKAGE = 'fmrireg', X_base_fixed, data_matrix, dmat_ran)
}

lss_compute_cpp <- function(Q_dmat_ran, residual_data) {
    .Call('_fmrireg_lss_compute_cpp', PACKAGE = 'fmrireg', Q_dmat_ran, residual_data)
}

mixed_solve_internal <- function(y_in, Z_in = NULL, K_in = NULL, X_in = NULL, method = "REML", bounds = as.numeric( c(1e-9, 1e9)), SE = FALSE, return_Hinv = FALSE) {
    .Call('_fmrireg_mixed_solve_internal', PACKAGE = 'fmrireg', y_in, Z_in, K_in, X_in, method, bounds, SE, return_Hinv)
}

neural_input_rcpp <- function(x, from, to, resolution) {
    .Call('_fmrireg_neural_input_rcpp', PACKAGE = 'fmrireg', x, from, to, resolution)
}

evaluate_regressor_convolution <- function(grid, onsets, durations, amplitudes, hrf_values, hrf_span, start, end, precision) {
    .Call('_fmrireg_evaluate_regressor_convolution', PACKAGE = 'fmrireg', grid, onsets, durations, amplitudes, hrf_values, hrf_span, start, end, precision)
}

