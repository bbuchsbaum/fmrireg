#define ARMA_DONT_PRINT_FAST_MATH_WARNING
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Templated helper to eliminate code duplication between Y and X processing
template<typename Matrix>
inline void whiten_matrix_impl(Matrix& M, const arma::vec& phi, 
                               bool exact_first_ar1, bool parallel) {
    const int n_time = M.nrow();
    const int n_cols = M.ncol();
    const int p = phi.n_elem;
    
    if (p == 0) return;  // No AR terms, nothing to do
    
    const double scale = (exact_first_ar1 && p == 1) ? 
                         std::sqrt(1.0 - phi[0] * phi[0]) : 1.0;
    
    #pragma omp parallel for if(parallel) schedule(static)
    for (int col = 0; col < n_cols; ++col) {
        // Use pointer arithmetic for better performance
        double* colPtr = &M(0, col);
        std::vector<double> prev(p, 0.0);
        
        for (int t = 0; t < n_time; ++t) {
            double orig = colPtr[t];
            double val = orig;
            
            // Apply AR filter: val = orig - sum(phi[k] * prev[k])
            for (int k = 0; k < p; ++k) {
                val -= phi[k] * prev[k];
            }
            
            // Apply exact first-sample scaling for AR(1) if requested
            if (t == 0 && exact_first_ar1 && p == 1) {
                val *= scale;
            }
            
            // Update the lag buffer: shift right and insert current original value
            for (int k = p - 1; k > 0; --k) {
                prev[k] = prev[k - 1];
            }
            prev[0] = orig;
            
            // Store the whitened value
            colPtr[t] = val;
        }
    }
}

//' AR(p) whitening of data and design matrices
//'
//' Applies a causal AR filter defined by `phi_coeffs` to both `Y` and `X`
//' matrices in place. The filter equation is:
//' 
//' v_t = y_t - sum(phi_k * y_{t-k}, k=1 to p)
//'
//' @param Y Numeric matrix of data (time x voxels)
//' @param X Numeric matrix of design (time x predictors)  
//' @param phi_coeffs Numeric vector of AR coefficients (length p)
//' @param exact_first_ar1 Logical, apply exact variance-normalizing scaling 
//'   of first sample for AR(1). For p > 1, no scaling is applied.
//' @param parallel Logical, enable OpenMP parallelization across columns
//'
//' @details 
//' The function assumes valid (stationary) AR coefficients are provided.
//' No checks for stationarity are performed. 
//' 
//' For exact_first_ar1 = TRUE and p = 1, the first residual is multiplied 
//' by sqrt(1 - phi^2) for proper variance normalization. This scaling is 
//' only applied for AR(1) models.
//'
//' @return List with components 'Y' and 'X' containing the whitened matrices
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List ar_whiten_inplace(Rcpp::NumericMatrix Y, Rcpp::NumericMatrix X,
                             const arma::vec& phi_coeffs, 
                             bool exact_first_ar1 = false,
                             bool parallel = true) {
    
    // Apply whitening to both matrices using the templated helper
    whiten_matrix_impl(Y, phi_coeffs, exact_first_ar1, parallel);
    whiten_matrix_impl(X, phi_coeffs, exact_first_ar1, parallel);
    
    // Return both matrices as a list for API compatibility
    return Rcpp::List::create(
        Rcpp::Named("Y") = Y,
        Rcpp::Named("X") = X
    );
}

//' AR(p) whitening with void return (no-copy version)
//'
//' More efficient version that modifies matrices in place without returning
//' copies. Use when you don't need the return values.
//'
//' @param Y Numeric matrix of data (time x voxels) - modified in place
//' @param X Numeric matrix of design (time x predictors) - modified in place
//' @param phi_coeffs Numeric vector of AR coefficients (length p)
//' @param exact_first_ar1 Logical, apply exact scaling of first sample for AR(1)
//' @param parallel Logical, enable OpenMP parallelization across columns
//'
//' @keywords internal
// [[Rcpp::export]]
void ar_whiten_void(Rcpp::NumericMatrix Y, Rcpp::NumericMatrix X,
                    const arma::vec& phi_coeffs, 
                    bool exact_first_ar1 = false,
                    bool parallel = true) {
    
    whiten_matrix_impl(Y, phi_coeffs, exact_first_ar1, parallel);
    whiten_matrix_impl(X, phi_coeffs, exact_first_ar1, parallel);
}
