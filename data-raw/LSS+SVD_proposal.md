Okay, this is excellent, actionable feedback on the R code for `ls_svd_rank1_hrf()`. The "minimal functional core" is a great target for a streamlined implementation, and the suggestions for speed and simplification are spot-on.

Let's use this to build a clear proposal for an engineer to implement both LS+SVD and then extend it to LS+SVD+1ALS.

---

**Proposal: Implementing LS+SVD and LS+SVD+1ALS for Rank-1 HRF Estimation**

**1. Objective:**

To implement two efficient CPU-based methods for rank-1 HRF estimation within an fMRI analysis framework:
1.  **LS+SVD:** A fast, closed-form method providing a good rank-1 approximation.
2.  **LS+SVD+1ALS:** An extension of LS+SVD that adds one iteration of Alternating Least Squares (ALS) to get closer to Pedregosa et al.'s original optimization objective with minimal additional compute.

These methods will serve as a baseline and a highly efficient default, respectively.

**2. Core Algorithm: LS+SVD (`ls_svd_engine`)**

This will be the foundational internal function.

*   **Inputs:**
    *   `X_list_proj`: List of `k` `n x d` *confound-projected* design matrices (one per condition).
    *   `Y_proj`: `n x v` *confound-projected* BOLD data matrix.
    *   `d`: Number of HRF basis functions/coefficients.
    *   `k`: Number of conditions.
    *   `lambda_init`: Ridge penalty for the initial full GLM solve (default: 1.0).
    *   `h_ref_shape_norm`: (Optional) `d x 1` *normalized* reference HRF shape (e.g., FIR coefficients of a canonical HRF, scaled to have `max(abs(h_ref_shape_norm))=1`) for sign alignment. If `NULL`, no sign alignment against a reference is performed beyond consistent SVD sign.
    *   `svd_backend`: (Optional, for future) String like "base_R", "RcppArmadillo", "RSpectra" to control SVD implementation.

*   **Steps:**
    1.  **Combine Design Matrices:** `Xbig = do.call(cbind, X_list_proj)` (forms `n x (dk)` matrix).
    2.  **Solve Full GLM (Ridge Regression):**
        *   `XtXbig = crossprod(Xbig)`
        *   `XtbY = crossprod(Xbig, Y_proj)`
        *   `Gamma_hat = chol2inv(chol(XtXbig + lambda_init * diag(ncol(Xbig)))) %*% XtbY` (using Cholesky decomposition for symmetric positive definite system).
    3.  **SVD per Voxel (Target for `Rcpp` optimization):**
        *   Initialize `H_out = matrix(0.0, d, v)` and `B_out = matrix(0.0, k, v)`.
        *   Loop `vx` from 1 to `v` (number of voxels):
            a.  `G_vx = matrix(Gamma_hat[, vx], nrow = d, ncol = k)`.
            b.  `svd_G_vx = svd(G_vx)` (or call to optimized SVD backend).
            c.  If `svd_G_vx$d[1]` is valid and > `epsilon_svd` (e.g., 1e-8):
                i.  `sqrt_s1 = sqrt(svd_G_vx$d[1])`
                ii. `h_vx = svd_G_vx$u[, 1] * sqrt_s1`
                iii. `b_vx = svd_G_vx$v[, 1] * sqrt_s1`
            d.  Else (singular value too small):
                i.  `h_vx = rep(0.0, d)`
                ii. `b_vx = rep(0.0, k)`
            e.  Store `h_vx` in `H_out[, vx]` and `b_vx` in `B_out[, vx]`.
    4.  **Identifiability (Simplified for FIR-like basis, assuming `B.mat=I`):**
        *   `scl_factors = apply(abs(H_out), 2, max)` (vector of length `v`).
        *   `flip_needed = rep(1.0, v)`.
        *   If `h_ref_shape_norm` is provided:
            `alignment_scores = crossprod(h_ref_shape_norm, H_out)` (results in `1 x v` matrix).
            `flip_needed[alignment_scores < 0 & scl_factors > epsilon_scale] = -1.0`.
        *   `effective_scl = pmax(scl_factors, epsilon_scale)` (e.g., `epsilon_scale = 1e-8`).
        *   `H_final = sweep(H_out, MARGIN = 2, STATS = flip_needed / effective_scl, FUN = "*")`.
        *   `B_final = sweep(B_out, MARGIN = 2, STATS = flip_needed * effective_scl, FUN = "*")`.
        *   `H_final[, scl_factors <= epsilon_scale] = 0` (zero out h if scale was negligible).
        *   `B_final[, scl_factors <= epsilon_scale] = 0`.
    5.  **Return:** `list(h = H_final, beta = B_final, Gamma_hat = Gamma_hat)`.

**3. Core Algorithm: LS+SVD+1ALS (`ls_svd_1als_engine`)**

This function builds upon LS+SVD by adding one ALS iteration.

*   **Inputs:**
    *   Same as `ls_svd_engine`, plus:
    *   `lambda_b`: Ridge penalty for β-update (default: 10.0).
    *   `lambda_h`: Ridge penalty for h-update (default: 1.0, assumes `R.mat=I`).
    *   `fullXtX_flag`: Boolean, whether to use full cross-terms `X_lᵀX_m` in h-update (default: `FALSE`).
    *   `X_list_proj`: (Already available from LS+SVD call or confound projection).

*   **Steps:**
    1.  **Initial Estimate:** Call `ls_svd_engine` to get `h_0 = H_final`, `b_0 = B_final`.
        *   *Crucially, for the ALS step, use `h_0` and `b_0` *before* their final scaling/flipping if the ALS updates are expected to refine these aspects naturally before a final identifiability step. Or, use the fully constrained `h_0, b_0` and apply identifiability again at the end. The latter is simpler.* Let's assume we use the fully constrained `h_0, b_0` from LS+SVD.
    2.  **Precompute for ALS (if not already available from `ls_svd_engine` context):**
        *   `XtX_list = lapply(X_list_proj, crossprod)` (list of `k` `d x d` matrices).
        *   `XtY_list = lapply(X_list_proj, function(X) crossprod(X, Y_proj))` (list of `k` `d x v` matrices).
        *   If `fullXtX_flag`: `XtX_full_list = matrix(list(), k, k)` storing all `crossprod(X_list_proj[[l]], X_list_proj[[m]])`.
    3.  **One ALS Iteration (Target for `Rcpp` optimization):**
        *   Initialize `H_als = matrix(0.0, d, v)`, `B_als = matrix(0.0, k, v)`.
        *   Use `h_current = h_0` from LS+SVD.
        *   **β-update (Loop `vx` from 1 to `v`):**
            a.  `h_vx = h_current[, vx]`.
            b.  `DhTy_vx = sapply(1:k, function(c) crossprod(h_vx, XtY_list[[c]][, vx]))`.
            c.  `Gmat_vx = matrix(0, k, k)`.
            d.  Loop `l` from 1 to `k`, loop `m` from 1 to `k`:
                `Gmat_vx[l,m] = crossprod(h_vx, if (fullXtX_flag) XtX_full_list[[l,m]] %*% h_vx else if (l==m) XtX_list[[l]] %*% h_vx else 0)`.
            e.  `B_als[, vx] = solve(Gmat_vx + lambda_b * diag(k), DhTy_vx)`.
        *   Use `b_current = B_als`.
        *   **h-update (Loop `vx` from 1 to `v`):**
            a.  `b_vx = b_current[, vx]`.
            b.  `lhs_vx = lambda_h * diag(d)` (assuming `R.mat=I`).
            c.  `rhs_vx = rep(0, d)`.
            d.  Loop `l` from 1 to `k`:
                `rhs_vx = rhs_vx + b_vx[l] * XtY_list[[l]][, vx]`.
                Loop `m` from (if `fullXtX_flag`) `1:k` else `l`: (handles symmetry/diagonal for `!fullXtX_flag`)
                    `term_matrix = if (fullXtX_flag) XtX_full_list[[l,m]] else XtX_list[[l]]`
                    `lhs_vx = lhs_vx + b_vx[l] * b_vx[m] * term_matrix`.
            e.  `H_als[, vx] = solve(lhs_vx, rhs_vx)`.
    4.  **Identifiability:** Apply scale and sign constraints to `H_als` and `B_als` (as in Step 4 of `ls_svd_engine`) to get `H_final`, `B_final`.
    5.  **Return:** `list(h = H_final, beta = B_final, h_ls_svd = h_0, beta_ls_svd = b_0)`.

**Reference Implementation Outline (R)**

```r
ls_svd_1als_engine <- function(X_list_proj, Y_proj,
                               lambda_init = 1,
                               lambda_b = 10,
                               lambda_h = 1,
                               fullXtX_flag = FALSE,
                               h_ref_shape_norm = NULL, ...) {
  init <- ls_svd_engine(X_list_proj, Y_proj,
                        lambda_init = lambda_init,
                        h_ref_shape_norm = h_ref_shape_norm, ...)
  XtX_list <- lapply(X_list_proj, crossprod)
  XtY_list <- lapply(X_list_proj, function(x) crossprod(x, Y_proj))
  # optional cross-terms if fullXtX_flag
  # perform one beta and h update
  # apply identifiability as in ls_svd_engine
  list(h = H_final, beta = B_final,
       h_ls_svd = init$h, beta_ls_svd = init$beta)
}
```

**4. Top-Level User Function (`fmrireg_cfals`)**

This function will wrap the engine calls.

*   **Inputs:**
    *   `fmri_data_obj`, `event_data_obj`, `hrf_basis_spec` (defining `d` and how to get `X.list`, and `Phi_matrix` for arbitrary bases), `confound_obj`.
    *   `method = "ls_svd_1als"` (default), or `"ls_svd_only"`, or `"cf_als_iterative"`.
    *   `max_alt`: (Used if `method = "cf_als_iterative"`, default 1 for "ls_svd_1als").
    *   Lambdas: `lambda_init`, `lambda_b`, `lambda_h`.
    *   `fullXtX_flag`.
    *   Reference HRF options (e.g., path to a file, or use default Glover).

*   **Logic:**
    1.  **Preprocessing:**
        *   Generate `X.list` from `event_data_obj` and `hrf_basis_spec`.
        *   Extract `Y` from `fmri_data_obj`.
        *   Extract/generate `Z` from `confound_obj`.
        *   Project out confounds from `Y` and `X.list` (using QR, `LAPACK=TRUE` for `qr.Q`).
        *   Determine `d`, `k`, `v`, `n`.
        *   Prepare `h_ref_shape_norm` based on `hrf_basis_spec` and reference HRF options.
    2.  **Call Engine based on `method`:**
        *   If `method == "ls_svd_only"`: Call `ls_svd_engine`.
        *   If `method == "ls_svd_1als"`: Call `ls_svd_1als_engine` (which internally calls `ls_svd_engine` for init).
        *   If `method == "cf_als_iterative"`: (Future) Implement a loop around the ALS steps from `ls_svd_1als_engine`, starting with `ls_svd_engine` init, for `max_alt` iterations.
    3.  **Package Results:** Create and return an `fmrireg_cfals_fit` object containing `h`, `beta`, and other relevant info (lambdas used, method, call, etc.). If `hrf_basis_spec` implied a non-FIR basis, also store/provide easy access to reconstructed HRF shapes `Phi_matrix %*% h`.

**5. Sprint Plan:**

**Sprint Goal:** Implement robust LS+SVD and LS+SVD+1ALS methods as core engines and integrate them into a user-friendly `fmrireg_cfals` function.

**Tickets:**

**Epic 1: Core LS+SVD Engine**
1.  **`CFALS-ENG-01`: Implement `ls_svd_engine` (R)**
    *   Implement the "LS+SVD" algorithm as detailed in Section 2.
    *   Focus on correctness and robustness (Cholesky solve, SVD epsilon checks, identifiability for FIR-like case).
    *   **DoD:** Function produces `h`, `beta`, `Gamma_hat` for synthetic FIR data. Passes unit tests.
2.  **`CFALS-ENG-02`: (Optional/Defer) Rcpp Optimization for `ls_svd_engine` SVD Loop**
    *   If profiling shows the R loop for SVD is a bottleneck for typical `v`.
    *   **DoD:** Rcpp version matches R version results and shows speedup.

**Epic 2: Core LS+SVD+1ALS Engine**
3.  **`CFALS-ENG-03`: Implement `ls_svd_1als_engine` (R)**
    *   Implement the "LS+SVD+1ALS" algorithm (Section 3), calling `ls_svd_engine` for initialization.
    *   Implement β-update and h-update loops, including `fullXtX_flag` logic.
    *   Apply final identifiability.
    *   **DoD:** Function produces `h`, `beta` for synthetic FIR data. Passes unit tests. Results show refinement over `ls_svd_engine` output.
4.  **`CFALS-ENG-04`: (Optional/Defer) Rcpp Optimization for `ls_svd_1als_engine` ALS Loops**
    *   If profiling shows R loops for β/h updates are bottlenecks.
    *   **DoD:** Rcpp version matches R version and shows speedup.

**Epic 3: `fmrireg` Integration & User Interface**
5.  **`CFALS-INT-01`: Confound Projection & Design Matrix Helpers**
    *   Robust QR projection: `project_confounds(Y, X.list, Z, lapack_qr=TRUE)`.
    *   `fmrireg` helper: `create_fmri_design(event_obj, basis_spec)` returning `X.list`, `d`, `k`, `Phi_matrix_for_recon`, `h_ref_shape_norm_for_basis`. (Handles FIR and potentially other `fmrireg` bases).
    *   **DoD:** Helper functions created and tested.
6.  **`CFALS-INT-02`: Implement `fmrireg_cfals` Top-Level Function**
    *   Wrapper function as detailed in Section 4.
    *   Handles `method` argument to call `ls_svd_engine` or `ls_svd_1als_engine`.
    *   Passes correctly prepared inputs (projected data, `d`, `k`, `h_ref_shape_norm`) to engines.
    *   **DoD:** User function callable; basic execution flow works for "ls_svd_only" and "ls_svd_1als" methods with FIR basis.
7.  **`CFALS-INT-03`: Output Object & Basic Methods**
    *   Define and implement `fmrireg_cfals_fit` S3/S4 class.
    *   Include `h_coeffs`, `beta_amps`, `method_used`, `lambdas`, `call`.
    *   If basis was non-FIR, include/enable easy reconstruction of HRF shapes.
    *   Implement `print()`, `summary()`, basic `plot()` (showing reconstructed HRF shapes).
    *   **DoD:** Output object structured; basic methods functional.

**Epic 4: Testing & Documentation**
8.  **`CFALS-TEST-01`: Unit & Integration Testing**
    *   Unit tests for engine functions with known small inputs.
    *   Integration tests for `fmrireg_cfals` with at least one real (public) fMRI dataset using an FIR basis.
    *   Validate against expected behavior (e.g., `fullXtX_flag` effect, lambda effects qualitatively).
    *   **DoD:** Good test coverage; results are plausible on real data.
9.  **`CFALS-DOC-01`: User Documentation**
    *   Man pages for `fmrireg_cfals` and `fmrireg_cfals_fit`.
    *   Clear explanation of `method` choices, lambdas, `fullXtX_flag`.
    *   Example script.
    *   **DoD:** Core documentation complete.

This proposal provides a clear path from the reviewed R sketch to a robust and flexible implementation suitable for an engineering team. The focus on lean core engines with a user-friendly wrapper allows for efficient development and future optimization.