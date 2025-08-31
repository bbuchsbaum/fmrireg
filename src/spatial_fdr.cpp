// Spatial FDR: Structure-adaptive weighted BH for multiple comparisons
// SABHA-style implementation for fmrireg
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <unordered_map>

using namespace Rcpp;
using namespace arma;

// Ensure groups are 1..G (compress factor/integers) ---------------------------
static inline IntegerVector compress_groups(const IntegerVector& g, int& G) {
  int n = g.size();
  std::vector<int> uniq;
  uniq.reserve(n);
  std::unordered_map<int,int> map;
  IntegerVector out(n);
  int next = 1;
  
  for (int i = 0; i < n; ++i) {
    // Handle NA values
    if (IntegerVector::is_na(g[i])) {
      out[i] = NA_INTEGER;
      continue;
    }
    
    int gi = g[i];
    auto it = map.find(gi);
    if (it == map.end()) {
      map.emplace(gi, next);
      out[i] = next;
      ++next;
    } else {
      out[i] = it->second;
    }
  }
  G = next - 1;
  return out;
}

// Groupwise pi0 counts (#p>tau) and sizes ------------------------------------
// [[Rcpp::export]]
List group_pi0_counts_cpp(NumericVector p, IntegerVector group, double tau) {
  int n = p.size();
  if (group.size() != n) {
    stop("group length must match p length.");
  }
  
  // Validate tau
  if (tau <= 0.0 || tau >= 1.0) {
    stop("tau must be in (0, 1)");
  }
  
  int G;
  IntegerVector g = compress_groups(group, G);
  
  if (G == 0) {
    stop("No valid groups found");
  }
  
  IntegerVector m_g(G);       // sizes
  IntegerVector tail_g(G);    // # p > tau
  
  for (int i = 0; i < n; ++i) {
    if (IntegerVector::is_na(g[i])) continue;
    
    int idx = g[i] - 1;
    if (idx < 0 || idx >= G) continue;
    
    ++m_g[idx];
    if (R_finite(p[i]) && p[i] > tau) {
      ++tail_g[idx];
    }
  }
  
  NumericVector pi0_raw(G);
  for (int j = 0; j < G; ++j) {
    if (m_g[j] == 0) { 
      pi0_raw[j] = NA_REAL; 
      continue; 
    }
    
    double est = (double)tail_g[j] / ((1.0 - tau) * (double)m_g[j]);
    
    // Bound estimate to [0, 1]
    if (!R_finite(est) || est < 0.0) est = 0.0;
    if (est > 1.0) est = 1.0;
    
    pi0_raw[j] = est;
  }
  
  return List::create(
    _["G"] = G,
    _["m_g"] = m_g,
    _["tail_g"] = tail_g,
    _["pi0_raw"] = pi0_raw,
    _["groups"] = g
  );
}

// Simple neighbor-mean smoothing of pi0 over adjacency -----------------------
// neighbors: list of integer vectors (1-based neighbor ids) length G
// [[Rcpp::export]]
NumericVector pi0_smooth_cpp(NumericVector pi0_raw, List neighbors,
                             double lambda, int iters) {
  int G = pi0_raw.size();
  
  if (neighbors.size() != G) {
    stop("neighbors list must have length G");
  }
  
  if (lambda < 0.0) {
    warning("lambda < 0 is not recommended; setting to 0");
    lambda = 0.0;
  }
  
  if (iters < 1) {
    warning("iters must be >= 1; setting to 1");
    iters = 1;
  }
  
  NumericVector cur = clone(pi0_raw);
  NumericVector nxt(G);
  
  for (int t = 0; t < iters; ++t) {
    for (int g = 0; g < G; ++g) {
      // Handle NA values
      if (!R_finite(cur[g])) {
        nxt[g] = cur[g];
        continue;
      }
      
      // Get neighbors for group g
      SEXP nb_sexp = neighbors[g];
      if (Rf_isNull(nb_sexp)) {
        nxt[g] = cur[g];
        continue;
      }
      
      NumericVector nb = as<NumericVector>(nb_sexp);
      double mean_nb = 0.0;
      int cnt = 0;
      
      for (int k = 0; k < nb.size(); ++k) {
        int j = nb[k] - 1;  // Convert to 0-based
        if (j >= 0 && j < G && j != g && R_finite(cur[j])) {
          mean_nb += cur[j];
          ++cnt;
        }
      }
      
      if (cnt > 0) {
        mean_nb /= (double)cnt;
        nxt[g] = (cur[g] + lambda * mean_nb) / (1.0 + lambda);
      } else {
        nxt[g] = cur[g];  // No valid neighbors
      }
    }
    cur = clone(nxt);
  }
  
  return cur;
}

// Weighted BH: returns reject mask + threshold + k ---------------------------
// We normalize weights inside so sum(w)=m, as required for weighted BH.
// [[Rcpp::export]]
List weighted_bh_cpp(NumericVector p, NumericVector w, double alpha) {
  int n = p.size();
  
  if (w.size() != n) {
    stop("w length must match p length.");
  }
  
  if (alpha <= 0.0 || alpha >= 1.0) {
    stop("alpha must be in (0, 1)");
  }
  
  // Clone weights to avoid modifying input
  NumericVector weights = clone(w);
  
  // Normalize weights to sum = n
  double sw = 0.0;
  int n_valid = 0;
  
  for (int i = 0; i < n; ++i) {
    if (!R_finite(p[i])) {
      weights[i] = 0.0;  // Zero weight for NA p-values
      continue;
    }
    
    double wi = weights[i];
    if (!R_finite(wi) || wi <= 0) {
      wi = 1e-8;  // Small positive weight
    }
    weights[i] = wi;
    sw += wi;
    n_valid++;
  }
  
  if (n_valid == 0) {
    // No valid p-values
    return List::create(
      _["reject"] = LogicalVector(n, false),
      _["threshold"] = 0.0,
      _["k"] = 0,
      _["w_norm"] = weights
    );
  }
  
  // Scale weights so sum = n_valid
  double scale = (sw > 0.0) ? ((double)n_valid / sw) : 1.0;
  for (int i = 0; i < n; ++i) {
    weights[i] *= scale;
  }
  
  // Build q = p / w and sort
  std::vector<std::pair<double, int>> arr;
  arr.reserve(n_valid);
  
  for (int i = 0; i < n; ++i) {
    if (!R_finite(p[i])) continue;
    
    double qi = (weights[i] > 0) ? (p[i] / weights[i]) : R_PosInf;
    if (!R_finite(qi)) qi = R_PosInf;
    
    arr.emplace_back(qi, i);
  }
  
  std::sort(arr.begin(), arr.end(),
            [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
              return a.first < b.first;
            });
  
  // Find largest k with q_(k) <= alpha * k / n_valid
  int kmax = 0;
  double tstar = 0.0;
  
  for (int k = 1; k <= n_valid; ++k) {
    double thresh = alpha * ((double)k / (double)n_valid);
    if (arr[k-1].first <= thresh) {
      kmax = k;
      tstar = thresh;
    }
  }
  
  // Build rejection vector
  LogicalVector reject(n, false);
  
  if (kmax > 0) {
    // Reject all i with p[i] <= tstar * w[i]
    for (int i = 0; i < n; ++i) {
      if (R_finite(p[i]) && weights[i] > 0 && p[i] <= tstar * weights[i]) {
        reject[i] = true;
      }
    }
  }
  
  return List::create(
    _["reject"] = reject,
    _["threshold"] = tstar,
    _["k"] = kmax,
    _["w_norm"] = weights
  );
}

// BH q-values on scaled p': q = p / w (step-up estimator) --------------------
// [[Rcpp::export]]
NumericVector bh_qvalues_scaled_cpp(NumericVector q) {
  int n = q.size();
  
  // Count valid values
  int n_valid = 0;
  for (int i = 0; i < n; ++i) {
    if (R_finite(q[i])) n_valid++;
  }
  
  if (n_valid == 0) {
    return NumericVector(n, NA_REAL);
  }
  
  // Build array of valid values with indices
  std::vector<std::pair<double, int>> arr;
  arr.reserve(n_valid);
  
  for (int i = 0; i < n; ++i) {
    if (!R_finite(q[i])) continue;
    
    double qi = q[i];
    if (qi < 0.0) qi = 0.0;
    if (qi > 1.0) qi = 1.0;
    
    arr.emplace_back(qi, i);
  }
  
  // Sort by q-value
  std::sort(arr.begin(), arr.end(),
            [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
              return a.first < b.first;
            });
  
  // Compute q-values
  NumericVector qv(n, NA_REAL);
  double minv = 1.0;
  
  for (int k = n_valid; k >= 1; --k) {
    double val = (double)n_valid * arr[k-1].first / (double)k;
    
    if (!R_finite(val)) val = 1.0;
    if (val < minv) minv = val;
    
    qv[arr[k-1].second] = std::min(1.0, std::max(0.0, minv));
  }
  
  return qv;
}