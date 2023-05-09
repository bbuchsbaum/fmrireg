#include <Rcpp.h>
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
