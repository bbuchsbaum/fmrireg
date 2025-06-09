#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// internal helper operating on Armadillo views of R matrices
static void ar_whiten_matrix(arma::mat& M, const arma::vec& phi, bool exact_first_ar1) {
    const arma::uword n_time = M.n_rows;
    const arma::uword n_cols = M.n_cols;
    const arma::uword p = phi.n_elem;

    #pragma omp parallel for
    for (arma::uword col = 0; col < n_cols; ++col) {
        std::vector<double> prev(p, 0.0);
        double scale = 1.0;
        if (exact_first_ar1 && p == 1) {
            scale = std::sqrt(1.0 - phi[0] * phi[0]);
        }
        for (arma::uword t = 0; t < n_time; ++t) {
            double orig = M(t, col);
            double val = orig;
            for (arma::uword k = 0; k < p; ++k) {
                val -= phi[k] * prev[k];
            }
            if (t == 0 && exact_first_ar1 && p == 1) {
                val *= scale;
            }
            if (p > 0) {
                for (int k = p - 1; k > 0; --k) {
                    prev[k] = prev[k - 1];
                }
                prev[0] = orig;
            }
            M(t, col) = val;
        }
    }
}

//' AR(p) whitening of data and design matrices
//'
//' Applies a causal AR filter defined by `phi_coeffs` to both `Y` and `X`
//' matrices in place.
//'
//' @param Y Numeric matrix of data (time x voxels)
//' @param X Numeric matrix of design (time x predictors)
//' @param phi_coeffs Numeric vector of AR coefficients (length p)
//' @param exact_first_ar1 Logical, apply exact scaling of first sample for AR(1)
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List ar_whiten_inplace(Rcpp::NumericMatrix Y, Rcpp::NumericMatrix X,
                       const arma::vec& phi_coeffs, bool exact_first_ar1 = false) {

    // Create Armadillo views on the underlying R matrices without copying
    arma::mat Ymat(Y.begin(), Y.nrow(), Y.ncol(), false);
    arma::mat Xmat(X.begin(), X.nrow(), X.ncol(), false);

    // Use the OpenMP-enabled helper to whiten both matrices column-wise
    ar_whiten_matrix(Ymat, phi_coeffs, exact_first_ar1);
    ar_whiten_matrix(Xmat, phi_coeffs, exact_first_ar1);

    // Return both matrices as a list
    return Rcpp::List::create(
        Rcpp::Named("Y") = Y,
        Rcpp::Named("X") = X
    );
}
