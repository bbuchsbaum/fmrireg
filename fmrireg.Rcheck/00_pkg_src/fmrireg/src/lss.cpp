#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List compute_residuals_cpp(const arma::mat& X_base_fixed, const arma::mat& data_matrix, const arma::mat& dmat_ran) {
    // Compute pseudoinverse using Armadillo's pinv (more efficient than MASS::ginv)
    arma::mat P_base_fixed = pinv(X_base_fixed);
    
    // Compute projection matrix
    arma::mat Q_base_fixed = arma::eye(X_base_fixed.n_rows, X_base_fixed.n_rows) - X_base_fixed * P_base_fixed;
    
    // Compute residuals and projected random effects design
    arma::mat residual_data = Q_base_fixed * data_matrix;
    arma::mat Q_dmat_ran = Q_base_fixed * dmat_ran;
    
    return List::create(
        Named("Q_dmat_ran") = Q_dmat_ran,
        Named("residual_data") = residual_data
    );
}

// [[Rcpp::export]]
arma::mat lss_compute_cpp(const arma::mat& Q_dmat_ran, const arma::mat& residual_data) {
    const uword n_timepoints = Q_dmat_ran.n_rows;
    const uword n_events = Q_dmat_ran.n_cols;
    const uword n_voxels = residual_data.n_cols;
    
    // Precompute sum of all regressors
    arma::vec total_sum = sum(Q_dmat_ran, 1);
    arma::mat beta_matrix(n_events, n_voxels);
    
    #pragma omp parallel for
    for (uword i = 0; i < n_events; ++i) {
        const arma::vec& c = Q_dmat_ran.col(i);
        // Subtract current column from total sum to get b
        arma::vec b = total_sum - c;
        
        // Compute v = c - b(b'b)^(-1)b'c
        const double b_norm = dot(b, b);
        arma::vec v;
        
        if (b_norm > 1e-10) {
            const double bc = dot(b, c);
            v = c - (bc/b_norm) * b;
        } else {
            v = c;
        }
        
        const double cvdot = dot(c, v);
        const double c_norm = dot(c, c);
        
        // Compute beta for current trial
        if (std::abs(cvdot) < 1e-5 * c_norm) {
            beta_matrix.row(i).zeros();
        } else {
            beta_matrix.row(i) = (v.t() * residual_data) / cvdot;
        }
    }
    
    return beta_matrix;
}
