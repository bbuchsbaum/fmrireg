#include <RcppArmadillo.h>
#include <roptim.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(roptim)]]

using namespace Rcpp;

// [[Rcpp::export]]
List mixed_solve_internal(NumericVector y_in,
                     Nullable<NumericMatrix> Z_in = R_NilValue,
                     Nullable<NumericMatrix> K_in = R_NilValue,
                     Nullable<NumericMatrix> X_in = R_NilValue,
                     std::string method = "REML",
                     NumericVector bounds = NumericVector::create(1e-9, 1e9),
                     bool SE = false,
                     bool return_Hinv = false) {
  
  // Convert y to arma::vec
  arma::vec y_full(y_in.begin(), y_in.size(), false);
  
  // Identify non-NA indices
  std::vector<int> not_NA_indices;
  for (size_t i = 0; i < y_full.n_elem; ++i) {
    if (!arma::is_finite(y_full[i])) continue;
    not_NA_indices.push_back(i);
  }
  
  int n_filtered = not_NA_indices.size();
  if (n_filtered == 0) stop("All y values are NA.");
  
  // Convert std::vector<int> to arma::uvec
  arma::uvec arma_not_NA_indices(not_NA_indices.size());
  for (size_t i = 0; i < not_NA_indices.size(); ++i) {
    arma_not_NA_indices[i] = static_cast<arma::uword>(not_NA_indices[i]);
  }
  
  // Set up X based on non-NA indices
  arma::mat X;
  if (X_in.isNull()) {
    X = arma::ones<arma::mat>(n_filtered, 1);
  } else {
    arma::mat X_temp = as<arma::mat>(X_in);
    X = X_temp.rows(arma_not_NA_indices);
  }
  
  int p = X.n_cols;
  if (p == 0) {
    p = 1;
    X = arma::ones<arma::mat>(n_filtered, 1);
  }
  
  // Set up Z based on non-NA indices
  arma::mat Z;
  if (Z_in.isNull()) {
    Z = arma::eye<arma::mat>(n_filtered, n_filtered);
  } else {
    arma::mat Z_temp = as<arma::mat>(Z_in);
    Z = Z_temp.rows(arma_not_NA_indices);
  }
  
  int m = Z.n_cols;
  if (m == 0) {
    m = 1;
    Z = arma::ones<arma::mat>(n_filtered, 1);
  }
  
  // Check dimensions
  if (Z.n_rows != n_filtered) stop("nrow(Z) != length of non-NA y");
  if (X.n_rows != n_filtered) stop("nrow(X) != length of non-NA y");
  
  // Handle K
  arma::mat K;
  bool K_null = false;
  if (K_in.isNull()) {
    K_null = true;
  } else {
    K = as<arma::mat>(K_in);
    if (K.n_rows != m || K.n_cols != m) stop("K dimensions mismatch with Z");
  }
  
  // Compute XtX
  arma::mat XtX = X.t() * X;
  int rank_X = arma::rank(XtX);
  if (rank_X < p) stop("X not full rank");
  
  arma::mat XtXinv = arma::inv_sympd(XtX);
  arma::mat S = arma::eye<arma::mat>(n_filtered, n_filtered) - X * XtXinv * X.t();
  
  // Determine spectral method
  std::string spectral_method;
  if (n_filtered <= m + p) {
    spectral_method = "eigen";
  } else {
    spectral_method = "cholesky";
    if (!K_null) {
      K.diag() += 1e-6;
      arma::mat cholK;
      bool chol_success = arma::chol(cholK, K);
      if (!chol_success) stop("K not positive semi-definite.");
    }
  }
  
  // Initialize variables for spectral decomposition
  arma::vec phi;
  arma::mat U;
  arma::vec theta;
  arma::mat Q;
  
  if (spectral_method == "cholesky") {
    arma::mat ZBt;
    if (K_null) {
      ZBt = Z;
    } else {
      arma::mat B;
      arma::chol(B, K);
      ZBt = Z * B.t();
    }
    
    arma::mat U_svd;
    arma::vec s_svd;
    arma::mat V_svd;
    arma::svd(U_svd, s_svd, V_svd, ZBt);
    
    U = U_svd;
    phi = arma::vec(n_filtered);
    phi.subvec(0, s_svd.n_elem - 1) = arma::square(s_svd);
    phi.subvec(s_svd.n_elem, n_filtered - 1).zeros();
    
    arma::mat SZBt = S * ZBt;
    arma::mat U_svd_SZBt;
    arma::vec s_svd_SZBt;
    arma::mat V_svd_SZBt;
    bool svd_success = arma::svd(U_svd_SZBt, s_svd_SZBt, V_svd_SZBt, SZBt);
    
    if (!svd_success) {
      SZBt += arma::randn<arma::mat>(SZBt.n_rows, SZBt.n_cols) * 1e-10;
      arma::svd(U_svd_SZBt, s_svd_SZBt, V_svd_SZBt, SZBt);
    }
    
    arma::mat QR_mat = arma::join_rows(X, U_svd_SZBt);
    arma::mat Q_full;
    arma::mat R_full;
    arma::qr(Q_full, R_full, QR_mat);
    
    Q = Q_full.cols(p, Q_full.n_cols - 1);
    
    int R_start = p;
    int R_end = std::min(p + m - 1, static_cast<int>(R_full.n_cols) - 1);
    arma::mat R = R_full.submat(R_start, R_start, R_end, R_end);
    
    arma::mat R_squared = arma::square(R);
    arma::mat R_squared_t = R_squared.t();
    
    arma::vec ans;
    bool solve_success = arma::solve(ans, R_squared_t, arma::square(s_svd_SZBt));
    
    if (!solve_success) {
      spectral_method = "eigen";
    } else {
      theta = arma::vec(n_filtered - p);
      // Replace those lines with:
      if (ans.n_elem > 0) {
        arma::uword end_idx = std::min<arma::uword>(ans.n_elem - 1, n_filtered - p - 1);
        theta.subvec(0, end_idx) = ans;
      }
      
      arma::uword start_idx = std::min<arma::uword>(ans.n_elem, n_filtered - p);
      theta.subvec(start_idx, n_filtered - p - 1).zeros();
    }
  }
  
  if (spectral_method == "eigen") {
    double offset = std::sqrt(static_cast<double>(n_filtered));
    arma::mat Hb;
    if (K_null) {
      Hb = Z * Z.t() + offset * arma::eye<arma::mat>(n_filtered, n_filtered);
    } else {
      Hb = Z * K * Z.t() + offset * arma::eye<arma::mat>(n_filtered, n_filtered);
    }
    
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, Hb);
    
    eigval = arma::reverse(eigval);
    eigvec = arma::fliplr(eigvec);
    
    phi = eigval - offset;
    if (phi.min() < -1e-6) stop("K not positive semi-definite.");
    U = eigvec;
    
    arma::mat SHbS = S * Hb * S;
    arma::vec eigval_SHbS;
    arma::mat eigvec_SHbS;
    arma::eig_sym(eigval_SHbS, eigvec_SHbS, SHbS);
    
    eigval_SHbS = arma::reverse(eigval_SHbS);
    eigvec_SHbS = arma::fliplr(eigvec_SHbS);
    
    theta = eigval_SHbS.subvec(0, n_filtered - p - 1) - offset;
    Q = eigvec_SHbS.cols(0, n_filtered - p - 1);
  }
  
  arma::vec y = y_full.elem(arma_not_NA_indices);
  arma::vec omega = Q.t() * y;
  arma::vec omega_sq = arma::square(omega);
  
  // Optimization
  double lambda_opt;
  double LL;
  int df;
  
  if (method == "ML") {
    struct ML_obj_fun : public roptim::Functor {
      int n;
      arma::vec theta;
      arma::vec omega_sq;
      arma::vec phi;
      
      ML_obj_fun(int n_, arma::vec theta_, arma::vec omega_sq_, arma::vec phi_) 
        : n(n_), theta(theta_), omega_sq(omega_sq_), phi(phi_) {}
      
      double operator()(const arma::vec &lambda_vec) override {
        double lambda = lambda_vec[0];
        arma::vec theta_lambda = theta + lambda;
        arma::vec phi_lambda = phi + lambda;
        
        if (arma::any(theta_lambda <= 0) || arma::any(phi_lambda <= 0)) {
          return std::numeric_limits<double>::infinity();
        }
        
        double sum_omega_sq_div_theta_lambda = arma::sum(omega_sq / theta_lambda);
        double obj = n * std::log(sum_omega_sq_div_theta_lambda) + arma::sum(arma::log(phi_lambda));
        return obj;
      }
    };
    
    ML_obj_fun obj_fun(n_filtered, theta, omega_sq, phi);
    roptim::Roptim<ML_obj_fun> opt("L-BFGS-B");
    
    opt.control.trace = 1;
    opt.control.maxit = 1000;
    opt.control.abstol = 1e-8;
    opt.control.reltol = 1e-8;
    
    arma::vec lower(1);
    lower[0] = bounds[0];
    opt.set_lower(lower);
    
    arma::vec upper(1);
    upper[0] = bounds[1];
    opt.set_upper(upper);
    
    arma::vec lambda_init = arma::vec(1);
    lambda_init[0] = std::max(bounds[0], std::min(arma::mean(phi), bounds[1]));
    
    opt.minimize(obj_fun, lambda_init);
    lambda_opt = lambda_init[0];
    df = n_filtered;
    LL = -0.5 * (opt.value() + df + df * std::log(2 * M_PI / df));
    
  } else {
    struct REML_obj_fun : public roptim::Functor {
      int n_p;
      arma::vec theta;
      arma::vec omega_sq;
      
      REML_obj_fun(int n_p_, arma::vec theta_, arma::vec omega_sq_) 
        : n_p(n_p_), theta(theta_), omega_sq(omega_sq_) {}
      
      double operator()(const arma::vec &lambda_vec) override {
        double lambda = lambda_vec[0];
        arma::vec theta_lambda = theta + lambda;
        
        if (arma::any(theta_lambda <= 0)) {
          return std::numeric_limits<double>::infinity();
        }
        
        double sum_omega_sq_div_theta_lambda = arma::sum(omega_sq / theta_lambda);
        double obj = n_p * std::log(sum_omega_sq_div_theta_lambda) + arma::sum(arma::log(theta_lambda));
        return obj;
      }
    };
    
    REML_obj_fun obj_fun(n_filtered - p, theta, omega_sq);
    roptim::Roptim<REML_obj_fun> opt("L-BFGS-B");
    
    opt.control.trace = 1;
    opt.control.maxit = 1000;
    opt.control.abstol = 1e-8;
    opt.control.reltol = 1e-8;
    
    arma::vec lower(1);
    lower[0] = bounds[0];
    opt.set_lower(lower);
    
    arma::vec upper(1);
    upper[0] = bounds[1];
    opt.set_upper(upper);
    
    arma::vec lambda_init = arma::vec(1);
    lambda_init[0] = std::max(bounds[0], std::min(arma::mean(theta), bounds[1]));
    
    opt.minimize(obj_fun, lambda_init);
    lambda_opt = lambda_init[0];
    df = n_filtered - p;
    LL = -0.5 * (opt.value() + df + df * std::log(2 * M_PI / df));
  }
  
  double Vu_opt = arma::sum(omega_sq / (theta + lambda_opt)) / df;
  double Ve_opt = lambda_opt * Vu_opt;
  
  arma::vec phi_lambda = phi + lambda_opt;
  arma::mat Ut_scaled = U.t();
  Ut_scaled.each_col() /= phi_lambda;
  arma::mat Hinv = U * Ut_scaled;
  
  arma::mat W = X.t() * Hinv * X;
  arma::vec beta = arma::solve(W, X.t() * Hinv * y);
  
  arma::mat KZt;
  if (K_null) {
    KZt = Z.t();
  } else {
    K.diag() += 1e-6;
    KZt = K * Z.t();
  }
  
  arma::mat KZt_Hinv = KZt * Hinv;
  arma::vec residual = y - X * beta;
  arma::vec u = KZt_Hinv * residual;
  
  List result;
  result["Vu"] = Vu_opt;
  result["Ve"] = Ve_opt;
  result["beta"] = beta;
  result["u"] = u;
  result["LL"] = LL;
  
  if (return_Hinv) {
    result["Hinv"] = Hinv;
  }
  
  if (SE) {
    arma::mat Winv = arma::inv_sympd(W);
    arma::vec beta_SE = arma::sqrt(Vu_opt * Winv.diag());
    arma::mat WW = KZt_Hinv * KZt.t();
    arma::mat WWW = KZt_Hinv * X;
    
    arma::vec u_SE;
    if (K_null) {
      arma::mat temp = WWW * Winv * WWW.t();
      u_SE = arma::sqrt(Vu_opt * (arma::ones<arma::vec>(m) - WW.diag() + temp.diag()));
    } else {
      arma::mat temp = WWW * Winv * WWW.t();
      u_SE = arma::sqrt(Vu_opt * (K.diag() - WW.diag() + arma::diagvec(temp)));
    }
    
    result["beta.SE"] = beta_SE;
    result["u.SE"] = u_SE;
  }
  
  return result;
}