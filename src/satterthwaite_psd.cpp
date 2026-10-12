// src/satterthwaite_psd.cpp
//
// Exact PSD-projected traces for the model-based Satterthwaite df
// (R/fmri_lm_variance.R, .fmri_lm_satterthwaite_banded()).
//
// For one voxel the residual covariance V is block diagonal over segments,
// each block V_s the symmetric banded Toeplitz matrix of tapered lag
// autocovariances (bandwidth L_s). The reference projects V onto the PSD
// cone: eigenvalues below tol = max(max|lambda(V)|, 1) * sqrt(eps) are set to
// 0, i.e. V+ = V - U D U' where (D, U) are the clipped eigenpairs. The
// Satterthwaite traces then follow from the unprojected traces with a
// low-rank correction (R = I - X G X', G = (X'X)^-1, Yz = U'X):
//
//   tr(R V+)      = tr(RV) - tr(D U'RU)
//   tr(R V+ R V+) = tr(RVRV) - 2 tr(D U'RVRU) + tr(D U'RU D U'RU)
//   U'RU   = I - F,                  F = Yz G Yz'
//   U'RVRU = D - F D - D F + Yz G (X'VX) G Yz'      (V U = U D)
//
// or, when more than half of the spectrum is clipped, directly from the kept
// eigenpairs: tr(R V+) = tr(D U'RU), tr(R V+ R V+) = tr((D U'RU)^2).
//
// Per voxel, without forming an m x m matrix:
//   0. Certificate: a banded Cholesky of V_s - c I for every segment, with c
//      above an upper bound of tol plus the Cholesky backward-error bound,
//      proves that nothing is clipped (O(m L^2)); the unprojected traces are
//      then exact. Otherwise:
//   1. V_s (symmetric Toeplitz, so centrosymmetric) is split exactly into its
//      symmetric and skew halves, each banded (see Block), and each half is
//      reduced to tridiagonal form T = Q' A Q by Givens rotations with bulge
//      chasing (Schwarz's band reduction, as in LAPACK dsbtrd), O(m^2 L)
//      flops on (L + 2) x m band storage. The rotations are recorded so that
//      Q' X_s is formed only if the voxel has eigenvalues to clip.
//   2. LAPACK dstebz (bisection) gives the extreme eigenvalues of every T,
//      hence tol over the whole of V as the reference defines it, and then
//      the eigenvalues on the needed side of tol (a Sturm count at tol and
//      bisection on (-inf, largest double below tol], i.e. lambda < tol).
//   3. LAPACK dstein (inverse iteration, reorthogonalised within clusters)
//      gives their eigenvectors Z of T, so U'X_s = Z' (Q' X_s).
// The band reduction, bisection and inverse iteration are backward stable, so
// the eigenvalues and the clipped invariant subspace are as accurate as the
// dense eigen() of the reference; clip decisions can differ only for an
// eigenvalue within rounding (~ eps ||V||) of tol, where the reference's own
// decision is equally arbitrary.
//
// status != 0 marks a voxel to recompute with the dense per-block helper in
// R: 1 = a LAPACK call reported info != 0, 2 = non-finite autocovariances,
// 3 = the Sturm count and the bisection disagree on the number clipped.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#ifndef FCONE
#define FCONE
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Rot {
  int p;       // rotation acts on rows/columns p and p + 1
  double c, s;
};

// Lower band storage with room for one bulge: A(i, k) (i >= k, i - k <= L + 1)
// at a[(i - k) + k * ld], ld = L + 2.
struct BandReducer {
  int m, L, ld;
  double* a;
  double* Y = nullptr;  // optional p x m matrix whose columns are rotated too
  int p = 0;

  BandReducer(int m_, int L_, double* a_) : m(m_), L(L_), ld(L_ + 2), a(a_) {}

  double& at(int i, int k) { return a[(i - k) + k * ld]; }

  // A <- P A P' with P = [[c, s], [-s, c]] on rows/columns (p, p + 1).
  void rotate(int p, double c, double s) {
    const int q = p + 1;
    const int W = L + 1;
    // Row parts left of the diagonal block: (A(p, k), A(q, k)), k < p.
    for (int k = std::max(0, q - W); k < p; ++k) {
      double* col = a + k * ld;
      double x = col[p - k], y = col[q - k];
      col[p - k] = c * x + s * y;
      col[q - k] = -s * x + c * y;
    }
    // Column parts below: (A(k, p), A(k, q)), k > q.
    const int kend = std::min(m - 1, p + W);
    double* cp = a + p * ld;
    double* cq = a + q * ld;
    for (int k = q + 1; k <= kend; ++k) {
      double x = cp[k - p], y = cq[k - q];
      cp[k - p] = c * x + s * y;
      cq[k - q] = -s * x + c * y;
    }
    const double app = cp[0], aqq = cq[0], apq = cp[1];
    const double cc = c * c, ss = s * s, cs = c * s;
    cp[0] = cc * app + 2.0 * cs * apq + ss * aqq;
    cq[0] = ss * app - 2.0 * cs * apq + cc * aqq;
    cp[1] = cs * (aqq - app) + (cc - ss) * apq;
  }

  // Zero A(q, k) (q = p + 1) against A(p, k). Returns false if already zero.
  bool annihilate(int p, int k, std::vector<Rot>* rots) {
    const int q = p + 1;
    double y = at(q, k);
    if (y == 0.0) return false;
    double x = at(p, k);
    double r = std::hypot(x, y);
    double c = x / r, s = y / r;
    rotate(p, c, s);
    at(p, k) = r;
    at(q, k) = 0.0;
    if (rots) rots->push_back(Rot{p, c, s});
    if (Y) {
      double* yp = Y + static_cast<size_t>(p) * this->p;
      double* yq = yp + this->p;
      for (int j = 0; j < this->p; ++j) {
        double u = yp[j], v = yq[j];
        yp[j] = c * u + s * v;
        yq[j] = -s * u + c * v;
      }
    }
    return true;
  }

  void reduce(std::vector<Rot>* rots) {
    if (L < 2) return;
    for (int j = 0; j + 2 < m; ++j) {
      for (int i = std::min(j + L, m - 1); i >= j + 2; --i) {
        if (!annihilate(i - 1, j, rots)) continue;
        // The rotation filled A(i + L, i - 1); chase it off the band.
        int br = i + L, bc = i - 1;
        while (br <= m - 1) {
          if (!annihilate(br - 1, bc, rots)) break;
          bc = br - 1;
          br += L;
        }
      }
    }
  }
};

struct Design {
  int n, p, nseg;
  const double* X;     // n x p, column major
  const double* G;     // p x p
  const int* rows;     // concatenated 0-based rows of the segments
  const int* seg_off;  // start of segment s in rows
  const int* seg_len;
  const int* nlag;     // L_s + 1
  const int* goff;     // start of segment s in gamma
};

// A segment's block V_s is symmetric Toeplitz, hence centrosymmetric
// (J V_s J = V_s, J the exchange matrix), and the orthogonal
// Q = [[I, I], [J, -J]] / sqrt(2) (with the middle row kept, scaled by
// sqrt(2) in the symmetric half, when m is odd) splits it exactly into a
// "symmetric" block of size ceil(m/2) and a "skew" block of size floor(m/2):
//   A_sym(i, j)  = g_|i-j| + g_(m-1-i-j),   A_skew(i, j) = g_|i-j| - g_(m-1-i-j)
//   (odd m) A_sym(q, j) = sqrt(2) g_(q-j),  A_sym(q, q) = g_0,   q = floor(m/2)
// with g_h = 0 for h > L. Both are banded with bandwidth L (the Hankel term
// only touches the trailing corner), so reducing the two halves costs half of
// reducing V_s, and the eigenpairs of V_s are those of the halves mapped back
// by Q. kind: 0 = whole segment (m = 1), 1 = symmetric half, 2 = skew half.
struct Block {
  int seg = 0, kind = 0, m = 0, L = 0;
  std::vector<double> d, e;
  std::vector<Rot> rots;
  bool recorded = false;
  double lo = 0.0, hi = 0.0;
};

const int kRecordLimit = 1024;  // record rotations for blocks up to this size

void fill_band(const Design& D, const Block& b, const double* gamma,
               std::vector<double>& band) {
  const int ms = D.seg_len[b.seg];
  const int Ls = D.nlag[b.seg] - 1;
  const int L = b.L;
  const int ld = L + 2;
  const int q = ms / 2;
  const bool odd = (ms % 2) == 1;
  const double* g = gamma + D.goff[b.seg];
  auto gf = [&](int h) { return h <= Ls ? g[h] : 0.0; };
  band.assign(static_cast<size_t>(ld) * b.m, 0.0);
  for (int j = 0; j < b.m; ++j) {
    double* col = &band[static_cast<size_t>(j) * ld];
    for (int d = 0; d <= L && j + d < b.m; ++d) {
      const int i = j + d;
      double v;
      if (b.kind == 0) {
        v = gf(d);
      } else if (b.kind == 1 && odd && i == q) {
        v = (j == q) ? g[0] : M_SQRT2 * gf(q - j);
      } else {
        const double hank = gf(ms - 1 - i - j);
        v = b.kind == 1 ? gf(d) + hank : gf(d) - hank;
      }
      col[d] = v;
    }
  }
}

// Q_b' X_s for the block (Q_b: the half of Q, then the band reduction),
// stored p x m (column i holds row i).
void rotated_design(const Design& D, const double* gamma, Block& b,
                    std::vector<double>& band, std::vector<double>& Y) {
  const int m = b.m, p = D.p;
  const int ms = D.seg_len[b.seg];
  const int q = ms / 2;
  const bool odd = (ms % 2) == 1;
  Y.assign(static_cast<size_t>(p) * m, 0.0);
  const int* rows = D.rows + D.seg_off[b.seg];
  const double r2 = M_SQRT1_2;
  for (int i = 0; i < m; ++i) {
    double* yi = &Y[static_cast<size_t>(i) * p];
    for (int j = 0; j < p; ++j) {
      const double* xj = D.X + static_cast<size_t>(j) * D.n;
      if (b.kind == 0) {
        yi[j] = xj[rows[i]];
      } else if (b.kind == 1 && odd && i == q) {
        yi[j] = xj[rows[q]];
      } else {
        const double a = xj[rows[i]], c = xj[rows[ms - 1 - i]];
        yi[j] = b.kind == 1 ? r2 * (a + c) : r2 * (a - c);
      }
    }
  }
  if (!b.recorded) {
    // Repeat the (deterministic) reduction, rotating Y as it goes.
    fill_band(D, b, gamma, band);
    BandReducer br(m, b.L, band.data());
    br.Y = Y.data();
    br.p = p;
    br.reduce(nullptr);
    return;
  }
  for (const Rot& r : b.rots) {
    double* yp = &Y[static_cast<size_t>(r.p) * p];
    double* yq = yp + p;
    for (int j = 0; j < p; ++j) {
      double x = yp[j], y = yq[j];
      yp[j] = r.c * x + r.s * y;
      yq[j] = -r.s * x + r.c * y;
    }
  }
}

int stebz(const char* range, int m, double vl, double vu, int il, int iu,
          const std::vector<double>& d, const std::vector<double>& e,
          std::vector<double>& w, std::vector<int>& iblock,
          std::vector<int>& isplit, int& found) {
  const double abstol = 2.0 * DBL_MIN;
  int nsplit = 0, info = 0;
  w.assign(m, 0.0);
  iblock.assign(m, 0);
  isplit.assign(m, 0);
  std::vector<double> work(4 * static_cast<size_t>(m));
  std::vector<int> iwork(3 * static_cast<size_t>(m));
  double edummy = 0.0;
  const double* ep = m > 1 ? e.data() : &edummy;
  F77_CALL(dstebz)(range, "B", &m, &vl, &vu, &il, &iu, &abstol, d.data(), ep,
                   &found, &nsplit, w.data(), iblock.data(), isplit.data(),
                   work.data(), iwork.data(), &info FCONE FCONE);
  return info;
}

// Sturm count: number of eigenvalues of the tridiagonal (d, e) below x.
int sturm_below(const std::vector<double>& d, const std::vector<double>& e,
                double x) {
  const int m = static_cast<int>(d.size());
  int count = 0;
  double q = d[0] - x;
  const double pivmin = DBL_MIN;
  if (std::fabs(q) < pivmin) q = -pivmin;
  if (q < 0) ++count;
  for (int i = 1; i < m; ++i) {
    q = d[i] - x - e[i - 1] * e[i - 1] / q;
    if (std::fabs(q) < pivmin) q = -pivmin;
    if (q < 0) ++count;
  }
  return count;
}

// Banded Cholesky of the segment's Toeplitz block minus shift * I; true if it
// completes with positive pivots. Lower band storage, ld = L + 1.
bool band_cholesky_ok(int m, int L, const double* g, double shift,
                      std::vector<double>& c) {
  const int ld = L + 1;
  c.assign(static_cast<size_t>(ld) * m, 0.0);
  for (int j = 0; j < m; ++j) {
    double* cj = &c[static_cast<size_t>(j) * ld];
    for (int d = 0; d <= L && j + d < m; ++d) cj[d] = g[d];
    cj[0] -= shift;
  }
  for (int j = 0; j < m; ++j) {
    double* cj = &c[static_cast<size_t>(j) * ld];
    const int k0 = std::max(0, j - L);
    // Pivot: A(j, j) - sum_k R(j, k)^2, R(j, k) at c[(j - k) + k ld].
    double s = cj[0];
    for (int k = k0; k < j; ++k) {
      const double r = c[(j - k) + static_cast<size_t>(k) * ld];
      s -= r * r;
    }
    if (!(s > 0.0) || !std::isfinite(s)) return false;
    const double r = std::sqrt(s);
    cj[0] = r;
    const int iend = std::min(m - 1, j + L);
    for (int i = j + 1; i <= iend; ++i) {
      double t = cj[i - j];
      for (int k = std::max(k0, i - L); k < j; ++k)
        t -= c[(i - k) + static_cast<size_t>(k) * ld] *
             c[(j - k) + static_cast<size_t>(k) * ld];
      cj[i - j] = t / r;
    }
  }
  return true;
}

// Projected traces for one voxel. Returns 0 on success.
// path: 0 = Cholesky certificate (nothing clipped), 1 = eigenvalues computed,
// nothing clipped, 2 = clipped-part correction, 3 = rebuilt from kept part.
int voxel_traces(const Design& D, const double* gamma, double& tr_RV,
                 double& tr_RVRV, int& nclip, int& path) {
  path = -1;
  const int p = D.p;
  // Certificate that nothing is clipped. ub >= ||V_s||_2 >= max|lambda| for
  // every block, so tol <= tol_hi. If V_s - c I has a Cholesky factor in
  // floating point, then lambda_min(V_s) > c - ||dA|| with the backward error
  // ||dA||_2 <= gamma_{L+2} (2L + 1) max_i a_ii (each entry of |R'||R| is at
  // most sqrt(a_ii a_jj), at most 2L + 1 per row); c exceeds tol_hi by that
  // bound plus 8 n eps ub, which also covers the eigenvalue error of the
  // dense reference, so neither the reference nor the exact path clips.
  double ub = 0.0;
  int n_all = 0, Lmax = 0;
  for (int s = 0; s < D.nseg; ++s) {
    const int L = D.nlag[s] - 1;
    const double* g = gamma + D.goff[s];
    double b = std::fabs(g[0]);
    if (!std::isfinite(g[0])) return 2;
    for (int h = 1; h <= L; ++h) {
      if (!std::isfinite(g[h])) return 2;
      b += 2.0 * std::fabs(g[h]);
    }
    ub = std::max(ub, b);
    n_all += D.seg_len[s];
    Lmax = std::max(Lmax, L);
  }
  {
    const double eps = DBL_EPSILON;
    const double tol_hi = std::max(ub, 1.0) * std::sqrt(eps) * (1.0 + 1e-8);
    const double margin =
        (4.0 * (Lmax + 2.0) * (2.0 * Lmax + 1.0) + 8.0 * n_all) * eps * ub;
    std::vector<double> chol;
    bool certified = true;
    for (int s = 0; s < D.nseg && certified; ++s) {
      certified = band_cholesky_ok(D.seg_len[s], D.nlag[s] - 1,
                                   gamma + D.goff[s], tol_hi + margin, chol);
    }
    if (certified) {
      nclip = 0;
      path = 0;
      return 0;
    }
  }

  std::vector<Block> blocks;
  for (int s = 0; s < D.nseg; ++s) {
    const int ms = D.seg_len[s];
    const int Ls = D.nlag[s] - 1;
    if (ms == 1) {
      Block b;
      b.seg = s; b.kind = 0; b.m = 1; b.L = 0;
      blocks.push_back(b);
      continue;
    }
    for (int kind = 1; kind <= 2; ++kind) {
      Block b;
      b.seg = s;
      b.kind = kind;
      b.m = kind == 1 ? ms - ms / 2 : ms / 2;
      b.L = std::min(Ls, b.m - 1);
      blocks.push_back(b);
    }
  }
  std::vector<double> band;
  double maxabs = 0.0;
  for (Block& sg : blocks) {
    fill_band(D, sg, gamma, band);
    BandReducer br(sg.m, sg.L, band.data());
    sg.recorded = sg.m <= kRecordLimit;
    br.reduce(sg.recorded ? &sg.rots : nullptr);
    sg.d.resize(sg.m);
    sg.e.resize(std::max(sg.m - 1, 0));
    for (int k = 0; k < sg.m; ++k) sg.d[k] = br.at(k, k);
    for (int k = 0; k + 1 < sg.m; ++k) sg.e[k] = br.at(k + 1, k);
    std::vector<double> w;
    std::vector<int> ib, is;
    int found = 0;
    if (stebz("I", sg.m, 0.0, 0.0, 1, 1, sg.d, sg.e, w, ib, is, found) != 0 ||
        found != 1) return 1;
    sg.lo = w[0];
    if (stebz("I", sg.m, 0.0, 0.0, sg.m, sg.m, sg.d, sg.e, w, ib, is, found) != 0 ||
        found != 1) return 1;
    sg.hi = w[0];
    maxabs = std::max(maxabs, std::max(std::fabs(sg.lo), std::fabs(sg.hi)));
  }
  const double tol = std::max(maxabs, 1.0) * std::sqrt(DBL_EPSILON);
  // lambda < tol  <=>  lambda <= below_tol.
  const double below_tol = std::nextafter(tol, -HUGE_VAL);

  int n_clip = 0;
  for (const Block& sg : blocks) {
    if (sg.lo < tol) n_clip += sturm_below(sg.d, sg.e, tol);
  }
  nclip = n_clip;
  path = 1;
  if (n_clip == 0) return 0;
  const bool correction = 2 * n_clip <= n_all;

  // Eigenpairs on the needed side of tol; Yz = U'X (k x p, row major).
  std::vector<double> lam, Yz;
  std::vector<double> w, Z, Y, work;
  std::vector<int> ib, is, iwork, ifail;
  int k_clip_found = 0;
  for (Block& sg : blocks) {
    const int m = sg.m;
    double vl, vu;
    if (correction) {
      if (!(sg.lo <= below_tol)) continue;
      vl = sg.lo - (std::fabs(sg.lo) + 1.0);
      vu = below_tol;
    } else {
      if (!(sg.hi >= tol)) {
        k_clip_found += m;
        continue;
      }
      vl = below_tol;
      vu = sg.hi + (std::fabs(sg.hi) + 1.0);
    }
    int found = 0;
    if (stebz("V", m, vl, vu, 0, 0, sg.d, sg.e, w, ib, is, found) != 0) return 1;
    if (correction) k_clip_found += found; else k_clip_found += m - found;
    if (found == 0) continue;
    Z.assign(static_cast<size_t>(m) * found, 0.0);
    work.assign(5 * static_cast<size_t>(m), 0.0);
    iwork.assign(m, 0);
    ifail.assign(found, 0);
    int info = 0;
    double edummy = 0.0;
    F77_CALL(dstein)(&m, sg.d.data(), m > 1 ? sg.e.data() : &edummy, &found,
                     w.data(), ib.data(), is.data(), Z.data(), &m, work.data(),
                     iwork.data(), ifail.data(), &info);
    if (info != 0) return 1;
    rotated_design(D, gamma, sg, band, Y);
    const size_t base = lam.size();
    lam.insert(lam.end(), w.begin(), w.begin() + found);
    Yz.resize((base + found) * p, 0.0);
    for (int i = 0; i < found; ++i) {
      double* out = &Yz[(base + i) * p];
      const double* z = &Z[static_cast<size_t>(i) * m];
      for (int r = 0; r < m; ++r) {
        const double zr = z[r];
        if (zr == 0.0) continue;
        const double* y = &Y[static_cast<size_t>(r) * p];
        for (int j = 0; j < p; ++j) out[j] += zr * y[j];
      }
    }
  }
  // Every count comes from a Sturm count at the same point; disagreement
  // means the bisection and the count saw different matrices.
  if (k_clip_found != n_clip) return 3;

  const int k = static_cast<int>(lam.size());
  // W1 = Yz G (k x p), F = W1 Yz'.
  std::vector<double> W1(static_cast<size_t>(k) * p, 0.0);
  for (int i = 0; i < k; ++i)
    for (int b = 0; b < p; ++b) {
      double acc = 0.0;
      for (int a = 0; a < p; ++a) acc += Yz[i * p + a] * D.G[a + b * p];
      W1[i * p + b] = acc;
    }
  std::vector<double> M1(static_cast<size_t>(k) * k);
  for (int i = 0; i < k; ++i)
    for (int j = 0; j < k; ++j) {
      double acc = 0.0;
      for (int b = 0; b < p; ++b) acc += W1[i * p + b] * Yz[j * p + b];
      M1[i + static_cast<size_t>(j) * k] = (i == j ? 1.0 : 0.0) - acc;
    }
  double t1 = 0.0, t3 = 0.0;
  for (int i = 0; i < k; ++i) {
    t1 += lam[i] * M1[i + static_cast<size_t>(i) * k];
    for (int j = 0; j < k; ++j)
      t3 += lam[i] * lam[j] * M1[i + static_cast<size_t>(j) * k] *
            M1[j + static_cast<size_t>(i) * k];
  }
  if (!correction) {
    tr_RV = t1;
    tr_RVRV = t3;
    path = 3;
    return 0;
  }

  // S = X'VX (p x p) from banded products.
  std::vector<double> S(static_cast<size_t>(p) * p, 0.0), VX;
  for (int s = 0; s < D.nseg; ++s) {
    const int m = D.seg_len[s], L = D.nlag[s] - 1;
    const int* rows = D.rows + D.seg_off[s];
    const double* g = gamma + D.goff[s];
    VX.assign(static_cast<size_t>(m) * p, 0.0);
    for (int j = 0; j < p; ++j) {
      const double* xj = D.X + static_cast<size_t>(j) * D.n;
      double* out = &VX[static_cast<size_t>(j) * m];
      for (int i = 0; i < m; ++i) {
        double acc = g[0] * xj[rows[i]];
        for (int h = 1; h <= L; ++h) {
          if (i + h < m) acc += g[h] * xj[rows[i + h]];
          if (i - h >= 0) acc += g[h] * xj[rows[i - h]];
        }
        out[i] = acc;
      }
    }
    for (int a = 0; a < p; ++a) {
      const double* xa = D.X + static_cast<size_t>(a) * D.n;
      for (int b = 0; b < p; ++b) {
        const double* vb = &VX[static_cast<size_t>(b) * m];
        double acc = 0.0;
        for (int i = 0; i < m; ++i) acc += xa[rows[i]] * vb[i];
        S[a + static_cast<size_t>(b) * p] += acc;
      }
    }
  }
  // SG = S G; diag_i of Yz G S G Yz' = sum_b (W1 SG)[i, b] Yz[i, b].
  std::vector<double> SG(static_cast<size_t>(p) * p, 0.0);
  for (int a = 0; a < p; ++a)
    for (int b = 0; b < p; ++b) {
      double acc = 0.0;
      for (int c = 0; c < p; ++c) acc += S[a + c * p] * D.G[c + b * p];
      SG[a + b * p] = acc;
    }
  double t2 = 0.0;
  for (int i = 0; i < k; ++i) {
    double hvh = 0.0;
    for (int b = 0; b < p; ++b) {
      double acc = 0.0;
      for (int a = 0; a < p; ++a) acc += W1[i * p + a] * SG[a + b * p];
      hvh += acc * Yz[i * p + b];
    }
    const double Fii = 1.0 - M1[i + static_cast<size_t>(i) * k];
    t2 += lam[i] * (lam[i] - 2.0 * lam[i] * Fii + hvh);
  }
  tr_RV = tr_RV - t1;
  tr_RVRV = tr_RVRV - 2.0 * t2 + t3;
  path = 2;
  return 0;
}

}  // namespace

//' PSD-projected Satterthwaite traces for voxels whose projection is active
//'
//' @param Gamma K x V tapered lag autocovariances (rows: segment s, lag h at
//'   goff[s] + h).
//' @param rows Concatenated 0-based rows of all segments; seg_off/seg_len
//'   locate segment s within them.
//' @param nlag L_s + 1 per segment; goff its first row in Gamma.
//' @param tr_RV,tr_RVRV Unprojected traces per voxel.
//' @return list(tr_RV, tr_RVRV, status, n_clipped, path); status 0 = exact
//'   result; path as documented at voxel_traces().
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List satterthwaite_psd_traces_cpp(Rcpp::NumericMatrix Gamma,
                                        Rcpp::NumericMatrix X,
                                        Rcpp::NumericMatrix G,
                                        Rcpp::IntegerVector rows,
                                        Rcpp::IntegerVector seg_off,
                                        Rcpp::IntegerVector seg_len,
                                        Rcpp::IntegerVector nlag,
                                        Rcpp::IntegerVector goff,
                                        Rcpp::NumericVector tr_RV,
                                        Rcpp::NumericVector tr_RVRV,
                                        int n_threads = 1) {
  const int nv = Gamma.ncol();
  const int K = Gamma.nrow();
  if (tr_RV.size() != nv || tr_RVRV.size() != nv)
    Rcpp::stop("tr_RV / tr_RVRV must have one value per voxel.");
  if (G.nrow() != X.ncol() || G.ncol() != X.ncol())
    Rcpp::stop("G must be p x p.");
  Design D{X.nrow(), X.ncol(), static_cast<int>(seg_len.size()), X.begin(),
           G.begin(), rows.begin(), seg_off.begin(), seg_len.begin(),
           nlag.begin(), goff.begin()};
  for (int s = 0; s < D.nseg; ++s) {
    if (seg_len[s] < 1 || nlag[s] < 1 || nlag[s] > seg_len[s] ||
        goff[s] + nlag[s] > K || seg_off[s] + seg_len[s] > rows.size())
      Rcpp::stop("Inconsistent segment description.");
  }
  for (int i = 0; i < rows.size(); ++i)
    if (rows[i] < 0 || rows[i] >= D.n) Rcpp::stop("Segment row out of range.");

  Rcpp::NumericVector out_rv = Rcpp::clone(tr_RV);
  Rcpp::NumericVector out_rvrv = Rcpp::clone(tr_RVRV);
  Rcpp::IntegerVector status(nv), nclip(nv), path(nv);
  double* prv = out_rv.begin();
  double* prvrv = out_rvrv.begin();
  int* pst = status.begin();
  int* pnc = nclip.begin();
  int* ppa = path.begin();
  const double* pg = Gamma.begin();

  int nth = n_threads;
#ifdef _OPENMP
  if (nth <= 0) nth = omp_get_max_threads();
#else
  nth = 1;
#endif
#pragma omp parallel for num_threads(nth) schedule(dynamic, 1)
  for (int v = 0; v < nv; ++v) {
    double a = prv[v], b = prvrv[v];
    int nc = 0, pa = -1;
    int st = voxel_traces(D, pg + static_cast<size_t>(v) * K, a, b, nc, pa);
    if (st == 0) {
      prv[v] = a;
      prvrv[v] = b;
    }
    pst[v] = st;
    pnc[v] = nc;
    ppa[v] = pa;
  }
  return Rcpp::List::create(Rcpp::Named("tr_RV") = out_rv,
                            Rcpp::Named("tr_RVRV") = out_rvrv,
                            Rcpp::Named("status") = status,
                            Rcpp::Named("n_clipped") = nclip,
                            Rcpp::Named("path") = path);
}
