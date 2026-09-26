// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#ifdef _OPENMP
  #include <omp.h>
#endif
using namespace arma;

// In-place unnormalised Walsh-Hadamard transform of one contiguous vector of
// length n (a power of two). H has entries +-1, H H' = n I and H = H'.
// The butterflies run in the classic stage order, so every output is formed
// by the same sequence of additions as a row-wise transform of a matrix.
static inline void fwht_inplace(double* x, const int n) {
  for (int len = 1; len < n; len <<= 1) {
    const int step = len << 1;
    for (int i = 0; i < n; i += step) {
      double* a = x + i;
      double* b = a + len;
      for (int j = 0; j < len; ++j) {
        const double u = a[j];
        const double v = b[j];
        a[j] = u + v;
        b[j] = u - v;
      }
    }
  }
}

static int next_pow2(int T) {
  int T2 = 1;
  while (T2 < T) T2 <<= 1;
  return T2;
}

// Threads for a column-parallel loop: n_threads <= 0 means the OpenMP
// default (itself capped by OMP_THREAD_LIMIT / OMP_NUM_THREADS); never more
// threads than columns.
static int sketch_threads(const int n_threads, const int ncol) {
#ifdef _OPENMP
  int nth = n_threads > 0 ? n_threads : omp_get_max_threads();
  if (nth > ncol) nth = ncol;
  return nth < 1 ? 1 : nth;
#else
  (void)n_threads;
  (void)ncol;
  return 1;
#endif
}

static void srht_check_plan(int T, const uvec& rows, const vec& signs,
                            const uvec& perm, const char* who) {
  if ((int)signs.n_elem != T) {
    Rcpp::stop("%s: length(signs) must equal the series length.", who);
  }
  if ((int)perm.n_elem != T) {
    Rcpp::stop("%s: length(perm) must equal the series length.", who);
  }
  if (perm.n_elem > 0 && perm.max() >= (uword)T) {
    Rcpp::stop("%s: perm contains out-of-bounds indices.", who);
  }
  if (rows.n_elem > 0 && rows.max() >= (uword)T) {
    Rcpp::stop("%s: rows contains out-of-bounds indices.", who);
  }
}

// SRHT apply: S M for M (T x k), with S = scale * R P H D (m x T).
//   D: random signs (T), H: unnormalised Walsh-Hadamard transform on the
//   zero-padded power-of-two length, P: permutation of the first T outputs,
//   R: selection of m of them.
// Every row of H D has squared norm T over the unpadded coordinates and
// E_D |(H D r)_i|^2 = ||r||^2, so scale = 1 / sqrt(m) gives the normalised
// sketch E ||S r||^2 = ||r||^2 (make_srht_plan() uses that scale). When T is
// a power of two, S S' = (T / m) I_m exactly.
//
// Each column is transformed independently in a thread-private contiguous
// buffer of the padded length (sign flip on load, row selection on store),
// in parallel over columns.
// [[Rcpp::export]]
arma::mat cpp_srht_apply(const arma::mat& M,
                         const arma::uvec& rows,
                         const arma::vec& signs,
                         const arma::uvec& perm,
                         const double scale,
                         const int n_threads = 0) {
  const int T = M.n_rows, K = M.n_cols;
  srht_check_plan(T, rows, signs, perm, "cpp_srht_apply");
  const int T2 = next_pow2(T);
  const int m = rows.n_elem;
  // Output row i reads transformed coordinate perm(rows(i)).
  std::vector<uword> src(m);
  for (int i = 0; i < m; ++i) src[i] = perm(rows(i));
  mat out(m, K);
  const double* sg = signs.memptr();
  const int nth = sketch_threads(n_threads, K);
#ifdef _OPENMP
#pragma omp parallel num_threads(nth)
#endif
  {
    std::vector<double> buf(T2);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int k = 0; k < K; ++k) {
      const double* x = M.colptr(k);
      for (int t = 0; t < T; ++t) buf[t] = x[t] * sg[t];
      for (int t = T; t < T2; ++t) buf[t] = 0.0;
      fwht_inplace(buf.data(), T2);
      double* o = out.colptr(k);
      for (int i = 0; i < m; ++i) o[i] = buf[src[i]] * scale;
    }
  }
  (void)nth;
  return out;
}

// Adjoint of cpp_srht_apply: S' B for B (m x k), returning T x k. With
// cpp_srht_apply it gives S S' = S (S' I_m), the Gram matrix of the sketch
// rows needed for the conditional variance of sketch-and-solve estimates.
// perm and rows are injective, so each sketch row scatters to a distinct
// coordinate of the padded buffer; H is symmetric.
// [[Rcpp::export]]
arma::mat cpp_srht_adjoint(const arma::mat& B,
                           const arma::uvec& rows,
                           const arma::vec& signs,
                           const arma::uvec& perm,
                           const double scale,
                           const int n_threads = 0) {
  const int T = signs.n_elem, K = B.n_cols;
  srht_check_plan(T, rows, signs, perm, "cpp_srht_adjoint");
  if (B.n_rows != rows.n_elem) {
    Rcpp::stop("cpp_srht_adjoint: nrow(B) must equal length(rows).");
  }
  const int T2 = next_pow2(T);
  const int m = rows.n_elem;
  std::vector<uword> dst(m);
  for (int i = 0; i < m; ++i) dst[i] = perm(rows(i));
  mat out(T, K);
  const double* sg = signs.memptr();
  const int nth = sketch_threads(n_threads, K);
#ifdef _OPENMP
#pragma omp parallel num_threads(nth)
#endif
  {
    std::vector<double> buf(T2);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int k = 0; k < K; ++k) {
      std::fill(buf.begin(), buf.end(), 0.0);
      const double* b = B.colptr(k);
      for (int i = 0; i < m; ++i) buf[dst[i]] += b[i] * scale;
      fwht_inplace(buf.data(), T2);
      double* o = out.colptr(k);
      for (int t = 0; t < T; ++t) o[t] = buf[t] * sg[t];
    }
  }
  (void)nth;
  return out;
}
