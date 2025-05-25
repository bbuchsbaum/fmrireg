This is an excellent, detailed, and actionable proposal for a high-performance robust GLM. The C++/Armadillo kernel sketch is particularly valuable and directly addresses the performance bottlenecks of the current `lmrob`-per-voxel approach. The focus on row-wise weighting is spot-on for fMRI data.

**Review of "Row-wise Huber IRLS" Proposal:**

**1. Statistical Correctness & Tweaks:**
*   **OLS Warm-start:** Correct. Adding a tiny ridge `+ 1e-8 * I` to `XᵀX` if `rcond(XᵀX)` is low is a good numerical safety net, though often not strictly needed if `X` is well-conditioned (e.g., after confound projection and if no extreme collinearity in task regressors).
*   **Scale `σ = MAD(vec(r))/0.6745`:**
    *   "One global `σ` is fine": Yes.
    *   "Use fastMAD (collapse) or `matrixStats::rowMedians` to avoid full vectorisation": `arma::median(arma::abs(arma::vectorise(residuals)))` in C++ is efficient as Armadillo handles vectorization well. `fastMAD` or similar might offer marginal gains in R but Armadillo's direct approach is good.
*   **`ρᵢ = rowMeans(r²)/σ²` (Squared Frame Energy):**
    *   "Robust against single-voxel noise": Yes.
    *   "If you fear localised artefacts use `rowMedians(abs(r)) / σ` instead": This is a good alternative for `u` directly (see below). The current `arma::mean(arma::square(R),1)` for `rtmp` and then `u = arma::sqrt(rtmp) / sigma` is fine and standard.
*   **Huber Weight `wᵢ = min(1, k / √ρᵢ)` (or `w = clamp(huber_k / u, 0.0, 1.0)`):**
    *   Correct. `k_huber = 1.345` for 95% efficiency is standard. Offering `k_huber = 2.0` as an option for less aggressive down-weighting is a good idea for users.
    *   `w.replace(arma::datum::nan, 1.0); // u=0 → w=1`: This is crucial. If a row's (standardized) residual energy `u` is zero, it's a perfect fit, so it should get full weight.
*   **Weighted Solve `(XᵀWX)β = XᵀWY`:** Correct.
*   **Update `σ` each iteration:** Standard. "Stop updating `σ` once it changes < 1 %" is a good heuristic for potentially faster convergence if `σ` stabilizes before `β`.
*   **Convergence `max|βᵗ - βᵗ⁻¹| < tol * (1 + max|βᵗ⁻¹|)`:** Relative tolerance is good practice.

**2. Why Existing R Code is Slow:**
*   This analysis is entirely correct. The overhead of `v` separate calls to a complex robust estimator like `lmrob()` (which does its own S-estimation, scale tuning, multiple IRLS iterations, etc.) is the killer.

**3. Drop-in C++/Armadillo Kernel (`cpp_fast_rlm`):**
*   **Signature:** `Rcpp::List cpp_fast_rlm(const mat& X, const mat& Y, double huber_k, int maxit, double tol)` - Perfect.
*   **OLS Warm-start:** `beta = arma::solve( X.t()*X, X.t()*Y )`. Correct. (Again, tiny ridge to `X.t()*X` if `X` could be rank-deficient).
*   **Initial `sigma`:** `arma::median( arma::abs(arma::vectorise(Y - X*beta)) ) / .6745;`. Correct and efficient.
*   **IRLS Loop:**
    *   `R = Y - X*beta;`: Efficient GEMM.
    *   `rtmp = arma::mean(arma::square(R),1);`: Row energies. `rtmp` is `n x 1`.
    *   `u = arma::sqrt(rtmp) / sigma;`: Standardized row "residual magnitude."
    *   `w = arma::clamp(huber_k / u, 0.0, 1.0);`: Huber weights.
    *   `w.replace(arma::datum::nan, 1.0);`: Correct handling for `u=0`.
    *   `sqw = arma::sqrt(w);`.
    *   `WX = X.each_col() % sqw;`: Efficient element-wise scaling of columns of `X` by `sqw`. *(Note: Armadillo's syntax for broadcasting a vector to scale rows would be `X.each_row() % sqw.t()` if `sqw` is a column vector, or more directly `mat WX = X; for(uword r=0; r<n; ++r) WX.row(r) *= sqw(r);` as sketched in my previous review. The `.each_col() % vec` syntax usually means each column of `X` is element-wise multiplied by `vec` if `vec` has `n` rows. This is correct if `sqw` is treated as a column vector and broadcast.)*
    *   `WY = Y.each_col() % sqw;`: Same logic for `Y`.
    *   `XtWX = WX.t()*WX;`.
    *   `XtWY = WX.t()*WY;`.
    *   `beta_new = arma::solve( XtWX, XtWY );`. (Consider adding `arma::solve_opts::fast` or `arma::solve_opts::likely_sympd` if `XtWX` is guaranteed positive definite).
    *   **Convergence Check:** `if( arma::abs(beta_new - beta).max() < tol * (1+arma::abs(beta).max()) )`. Excellent relative tolerance check.
    *   `beta = std::move(beta_new);`. Good C++ practice for efficiency.
    *   Update `sigma`.
*   **Return List:** `beta`, `weights`. Good. `sigma` should also be returned.

**4. R Wrapper (`fast_rlm`) and Integration:**
*   Simple R wrapper `fast_rlm(...) { cpp_fast_rlm(...) }` is perfect.
*   Integration into `runwise_rlm` / `chunkwise_rlm`:
    *   `fit <- fast_rlm(X_proj, Y_proj, ...)`
    *   `beta_mat <- fit$beta`
    *   `row_w <- fit$weights`
    *   `sigma_hat_sq_global = fit$sigma^2` (assuming `sigma` is returned from `cpp_fast_rlm`).
*   **Covariance `Var(β̂) = σ̂² (XᵀWX)⁻¹`:**
    *   `σ̂² = (1/(n-k)) * sum(w_i * r_i_sq_scalar)` where `r_i_sq_scalar` would be `sum_voxels( (Y_iv - X_i β_v)^2 ) / n_voxels` for time point `i`.
    *   Alternatively, a common robust variance estimate is based on the final `sigma` from the IRLS: `σ̂²_global = final_sigma_from_IRLS^2`.
    *   Then `Var(β̂_v)` (per voxel) is `σ̂²_global * solve(XtWX_final)`. `XtWX_final` is common for all voxels. This is a huge advantage.

**5. Extra Polish:**
*   **Leverage-based IRLS skip:** `if max(abs(XᵀY_v)) < ε` then `β_v = 0`. This is a good heuristic for very low signal voxels. Could be done after OLS warm-start for a voxel mask.
*   **GPU GEMM:** Future.
*   **Row-wise temporary censoring (`w_i = 0`):** Interesting idea, makes it behave like scrubbing for truly bad frames. Could be an option.
*   **Configurable `ψ`:** For advanced users, allowing other weight functions (Tukey) is a good extension.
*   **Save weight history:** Useful for diagnostics.

**Conclusion of Review:**
This is an excellent, well-thought-out, and highly practical proposal. The C++/Armadillo sketch is solid and directly implementable. The approach correctly identifies the main performance bottleneck in fMRI robust regression and provides an elegant and efficient solution. The "Extra Polish" ideas are all valuable additions for a production-grade tool.

---

**Proposal: High-Performance Row-Weighted Robust GLM (`fmri_fast_rlm`)**

**1. Objective:**

To implement a high-performance, robust General Linear Model (GLM) estimation function for fMRI data, `fmri_fast_rlm`. This function will replace the current voxel-wise `robustbase::lmrob` loop in `fmri_rlm` with a highly optimized C++/Armadillo backend that leverages row-wise Iteratively Reweighted Least Squares (IRLS) with Huber weights. This approach targets frame-wise outliers common in fMRI, achieving robustness with computational efficiency comparable to standard OLS GLM.

**2. Core Algorithm (`cpp_fast_rlm` - C++/Armadillo Backend):**

*   **Inputs:**
    *   `X_proj`: `n x k` dense design matrix (confounds projected out, run-specific).
    *   `Y_proj`: `n x v` dense BOLD data matrix (confounds projected out, run-specific).
    *   `huber_k_const`: Huber tuning constant (default: 1.345).
    *   `max_iterations`: Max IRLS iterations (default: 20, often converges in < 6-10).
    *   `tolerance`: Convergence tolerance for `β` (default: 1e-5).
*   **Outputs (List/Struct):**
    *   `beta_estimates`: `k x v` matrix of robustly estimated coefficients.
    *   `final_row_weights`: `n x 1` vector of final Huber weights for time points.
    *   `final_sigma_estimate`: Final robust scale estimate (scalar).
    *   `iterations_performed`: Number of iterations taken to converge.
*   **Steps:**
    1.  **OLS Warm-Start:**
        *   `beta_ols = solve(X_projᵀX_proj + tiny_ridge*I, X_projᵀY_proj)`.
        *   `residuals_ols = Y_proj - X_proj %*% beta_ols`.
        *   `sigma_current = median(abs(vectorise(residuals_ols))) / 0.6745`.
        *   `beta_current = beta_ols`.
    2.  **IRLS Loop (until `max_iterations` or convergence):**
        a.  `beta_previous = beta_current`.
        b.  `residuals_iter = Y_proj - X_proj %*% beta_current`.
        c.  `row_residual_energy_scaled = arma::mean(square(residuals_iter), 1) / (sigma_current^2)`. (This is `ρ_i` from proposal, `n x 1`). Ensure `sigma_current` is not zero.
        d.  `u_scaled_magnitude = sqrt(abs(row_residual_energy_scaled))`.
        e.  `row_weights = clamp(huber_k_const / u_scaled_magnitude, 0.0, 1.0)`.
        f.  `row_weights.replace(NaN, 1.0)` (handles `u_scaled_magnitude = 0`).
        g.  `sqrt_row_weights = sqrt(row_weights)`.
        h.  `WX = X_proj.each_col() % sqrt_row_weights` (or row-wise broadcast: `X_proj_row_i * sqrt_row_weights_i`).
        i.  `WY = Y_proj.each_col() % sqrt_row_weights` (or row-wise broadcast).
        j.  `XtWX = WX.t() * WX`.
        k.  `XtWY = WX.t() * WY`.
        l.  `beta_current = solve(XtWX + tiny_ridge_robust*I, XtWY)`. (Tiny ridge for safety in weighted solve).
        m. `residuals_updated_beta = Y_proj - X_proj %*% beta_current`.
        n.  `sigma_current = median(abs(vectorise(residuals_updated_beta))) / 0.6745`. (Option: stop updating `sigma` after a few iterations or if it stabilizes).
        o.  If `max_abs_diff(beta_current - beta_previous) < tolerance * (1 + max_abs(beta_previous))`, then `break`.
    3.  **Return** `beta_current`, `row_weights`, `sigma_current`, `iterations_performed`.

**3. R Wrapper (`fast_rlm`):**

*   A simple R function that calls `cpp_fast_rlm`.
    `fast_rlm(X, Y, huber_k = 1.345, maxit = 20, tol = 1e-5) { cpp_fast_rlm(X, Y, huber_k, maxit, tol) }`

**4. Integration into `fmri_rlm` (User-Facing Function):**

*   The existing `fmri_rlm(formula, block, dataset, ..., robust=TRUE)` function will be modified.
*   If `robust=TRUE`:
    1.  Standard `fmri_lm` preprocessing: build `X_model` (per run/block), get `Y_data`.
    2.  Project out any global confounds from `X_model` and `Y_data` to get `X_proj`, `Y_proj`.
    3.  Call `fit_results <- fast_rlm(X_proj, Y_proj, huber_k, maxit, tol)`.
    4.  Store `fit_results$beta_estimates`.
    5.  Calculate coefficient standard errors and t-statistics:
        *   `XtWX_final` can be reconstructed using `fit_results$final_row_weights` and `X_proj`.
        *   `Var_beta_hat = fit_results$final_sigma_estimate^2 * solve(XtWX_final)`. (This covariance is common for all voxels).
    6.  Proceed with contrast calculations and tidying outputs, similar to `fmri_lm`.
*   If `robust=FALSE`, it falls back to standard `fmri_lm` OLS.

**5. Sprint Plan:**

**Sprint Goal:** Implement and integrate the high-performance `cpp_fast_rlm` backend into `fmri_rlm`, providing a fast and robust alternative to per-voxel `lmrob`.

**Tickets:**

**Epic 1: Core C++/Armadillo Robust GLM Engine**

1.  **`RLM-CPP-01`: Implement `cpp_fast_rlm` Backend**
    *   **Task:** Create the C++/Armadillo function `cpp_fast_rlm` as detailed in "Core Algorithm."
    *   Implement OLS warm-start, IRLS loop with Huber row-weights, scale estimation, and convergence checks.
    *   Ensure efficient matrix operations using Armadillo.
    *   Return `beta`, `final_row_weights`, `final_sigma_estimate`, `iterations_performed`.
    *   **DoD:** C++ function compiles and passes basic tests with small synthetic `X` and `Y`.

2.  **`RLM-CPP-02`: Unit Testing for `cpp_fast_rlm`**
    *   **Task:** Create R/C++ unit tests (e.g., using `testthat` calling Rcpp function):
        *   Test against known results for a simple 1-voxel case (can compare to `robustbase::lmrob` with row weights if such a mode exists, or a manual IRLS in R).
        *   Test convergence properties, edge cases (e.g., all weights = 1, some weights = 0).
        *   Test handling of `huber_k`, `maxit`, `tol`.
    *   **DoD:** Comprehensive unit tests pass.

**Epic 2: R Integration and `fmri_rlm` Update**

3.  **`RLM-R-01`: Create `fast_rlm` R Wrapper**
    *   **Task:** Implement the simple R wrapper function `fast_rlm` that calls `cpp_fast_rlm`.
    *   **DoD:** R wrapper function created and callable from R.

4.  **`RLM-FMRI-02`: Integrate `fast_rlm` into `fmri_rlm`**
    *   **Task:** Modify the existing `fmri_rlm` function:
        *   When `robust=TRUE`, after standard preprocessing (model matrix, confound projection), call `fast_rlm`.
        *   Extract `beta_estimates`.
        *   Compute `XtWX_final` using `final_row_weights` from `fast_rlm` output.
        *   Compute `Var_beta_hat = final_sigma_estimate^2 * solve(XtWX_final)`.
        *   Adapt downstream contrast calculation and result formatting to use these robust estimates.
    *   **DoD:** `fmri_rlm(..., robust=TRUE)` runs using the new backend and produces plausible beta estimates and standard errors.

5.  **`RLM-FMRI-03`: Testing `fmri_rlm` with Robust Backend**
    *   **Task:** Test the updated `fmri_rlm` on at least one real fMRI dataset:
        *   Compare results (beta magnitudes, significant voxels) qualitatively with `robust=FALSE` (OLS) and, if possible, with the old `lmrob`-based `fmri_rlm` on a small subset for sanity check.
        *   Verify performance (speed) improvement.
        *   Check for correctness of contrast calculations and output structure.
    *   **DoD:** `fmri_rlm` with new robust backend is significantly faster and produces stable, interpretable results.

**Epic 3: Documentation & Polish**

6.  **`RLM-DOC-01`: Update `fmri_rlm` Documentation**
    *   **Task:** Update documentation for `fmri_rlm` to explain the new robust estimation method, its assumptions (row-wise outliers), and parameters (`huber_k`, `maxit`, `tol` if exposed).
    *   Explain the benefits (speed, robustness to frame-wise spikes).
    *   **DoD:** User documentation is updated and clear.

7.  **`RLM-POLISH-02`: (Optional Polish) Implement "Extra Polish" Ideas**
    *   **Task (if time permits):**
        *   Leverage-based IRLS skip heuristic in `cpp_fast_rlm` (after OLS warm-start, identify voxels with very low `XᵀY_v` signal and set their `β=0`, skipping IRLS for them).
        *   Option to return weight history from `cpp_fast_rlm` for diagnostics.
    *   **DoD:** Selected polish features implemented and tested.

**Sprint Review Focus:**
*   Demonstrable massive speedup of `fmri_rlm(..., robust=TRUE)`.
*   Correctness of beta estimates and their standard errors from the new backend.
*   Robustness of the implementation (handling of edge cases, numerical stability).
*   Clear documentation for users.

This plan provides a clear path to a significantly improved robust GLM capability. The key is the efficient C++/Armadillo backend.