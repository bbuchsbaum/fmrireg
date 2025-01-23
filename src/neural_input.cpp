#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List neural_input_rcpp(List x, double from, double to, double resolution) {
  int n = (to - from) / resolution;
  NumericVector time(n);
  NumericVector out(n);
  NumericVector ons = x["onsets"];
  NumericVector dur = x["duration"];
  NumericVector amp = x["amplitude"];
  
  for (int i = 0; i < ons.length(); i++) {
    double on = ons[i];
    double d = dur[i];
    int startbin = (int) ((on - from) / resolution) + 1;
    if (d > 0) {
      int endbin = (int) ((on - from) / resolution + d / resolution) + 1;
      for (int j = startbin; j <= endbin; j++) {
        out[j-1] += amp[i];
      }
    } else {
      out[startbin-1] += amp[i];
    }
  }
  
  for (int i = 0; i < n; i++) {
    time[i] = from + (i + 0.5) * resolution;
  }
  
  List result;
  result["time"] = time;
  result["neural_input"] = out;
  return result;
}


using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix evaluate_regressor_convolution(NumericVector grid,
                                             NumericVector onsets,
                                             NumericVector durations,
                                             NumericVector amplitudes,
                                             NumericMatrix hrf_values,
                                             double hrf_span,
                                             double start,
                                             double end,
                                             double precision) {
  int ngrid = grid.size();
  int nonsets = onsets.size();
  int n_hrf = hrf_values.nrow();
  int nfine = int((end - start) / precision) + 1;
  NumericVector finegrid(nfine);
  for (int i = 0; i < nfine; i++) {
    finegrid[i] = start + i * precision;
  }
  arma::vec neural_input(nfine, arma::fill::zeros);
  // Generate neural input signal
  for (int i = 0; i < nonsets; i++) {
    double onset = onsets[i];
    double duration = durations[i];
    double amplitude = amplitudes[i];
    int start_idx = int((onset - start) / precision);
    int end_idx = int((onset + duration - start) / precision);
    if (start_idx < 0) start_idx = 0;
    if (end_idx >= nfine) end_idx = nfine - 1;
    for (int j = start_idx; j <= end_idx; j++) {
      neural_input[j] += amplitude;
    }
  }
  // Convolve neural input with HRF for each basis
  int nbasis = hrf_values.ncol();
  arma::mat conv_result(nfine, nbasis);
  for (int b = 0; b < nbasis; b++) {
    arma::vec hrf_b = hrf_values(_, b);
    arma::vec conv_b = arma::conv(neural_input, hrf_b);
    // Trim the convolution result to match the fine grid size
    conv_result.col(b) = conv_b.subvec(0, nfine - 1);
  }
  // Interpolate conv_result to grid
  NumericMatrix outmat(ngrid, nbasis);
  for (int b = 0; b < nbasis; b++) {
    arma::vec conv_b = conv_result.col(b);
    for (int i = 0; i < ngrid; i++) {
      double t = grid[i];
      if (t <= finegrid[0]) {
        outmat(i, b) = conv_b[0];
      } else if (t >= finegrid[nfine - 1]) {
        outmat(i, b) = conv_b[nfine - 1];
      } else {
        // Linear interpolation
        int idx = int((t - start) / precision);
        if (idx >= nfine - 1) idx = nfine - 2;
        double t1 = finegrid[idx];
        double t2 = finegrid[idx + 1];
        double y1 = conv_b[idx];
        double y2 = conv_b[idx + 1];
        double alpha = (t - t1) / (t2 - t1);
        outmat(i, b) = y1 + alpha * (y2 - y1);
      }
    }
  }
  return outmat;
}