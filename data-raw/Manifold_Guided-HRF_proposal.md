This is a superb set of final clarifications and an excellent "API polish" for the LSS component. The decisions made are pragmatic, robust, and maintain the core elegance and efficiency of the manifold HRF approach.

Here's a full proposal incorporating these latest refinements, structured for an engineering team.

---

**Proposal: Manifold-Guided HRF Estimation and Trial-wise Deconvolution (M-HRF-LSS)**

**1. Executive Summary:**

This proposal outlines an advanced fMRI analysis pipeline, "Manifold-Guided HRF Estimation and Trial-wise Deconvolution (M-HRF-LSS)." It combines a novel HRF estimation technique based on manifold learning with a highly optimized Least Squares Separate (LSS) kernel for single-trial beta estimation. The M-HRF-LSS pipeline aims to provide:
1.  Spatially smooth and physiologically plausible voxel-wise HRF estimates.
2.  Associated per-condition activation amplitudes.
3.  Accurate per-trial activation amplitudes.
All components are designed for computational efficiency, primarily relying on closed-form solutions, BLAS operations, and sparse matrix algebra, making it suitable for large-scale whole-brain analyses on standard CPU hardware.

**2. Core Components & Algorithmic Workflow:**

**Component 0: HRF Manifold Construction (Once per Study/HRF Library)**

*   **Objective:** Create a low-dimensional, smooth representation (manifold) of plausible HRF shapes.
*   **Inputs:**
    *   `L_library`: `p x N` matrix of `N` candidate HRF shapes, each sampled at `p` time points (e.g., from a physiological model grid, half-cosine set, FLOBS).
    *   `m_manifold_dim`: Target dimensionality of the manifold (typically 3-4).
    *   (Optional) `k_nn_for_W`: k-Nearest Neighbors for sparse affinity matrix `W` if `N` is very large.
*   **Steps:**
    1.  **Affinity Matrix `W` (Self-Tuning):**
        *   Compute pairwise Euclidean distances `dists_ij` between all HRFs in `L_library`.
        *   For each HRF `i`, find its `k_local_nn`-th nearest neighbor distance `σ_i` (e.g., `k_local_nn=3`).
        *   `W_ij = exp(-dists_ij² / (σ_i * σ_j))`. If `N` is large, `W` can be sparsified (e.g., k-NN graph).
    2.  **Markov Matrix `S`:** `D_inv = diag(1/rowSums(W))`, `S = D_inv %*% W`.
    3.  **Diffusion Map Coordinates `Φ_raw`:** Compute top `m+1` eigenvectors of `S` (e.g., via `RSpectra::eigs_sym`).
        `Φ_raw = eig_S$vectors`.
    4.  **Manifold Coordinates `Φ`:** `Φ = Φ_raw[, 2:(m_manifold_dim + 1)]` (an `N x m_manifold_dim` matrix, discarding the trivial constant eigenvector).
    5.  **HRF Reconstructor `B_reconstructor`:** `B_reconstructor = L_library %*% Φ %*% solve(crossprod(Φ) + 1e-8 * diag(m_manifold_dim))` (a `p x m_manifold_dim` matrix). This `B_reconstructor` maps `m`-dimensional manifold coordinates `ξ` back to a `p`-sample HRF shape: `h_shape = B_reconstructor %*% ξ`.

**Component 1: Voxel-wise HRF Manifold Coordinate & Condition Amplitude Estimation**

*   **Objective:** For each voxel, estimate its `m`-dimensional HRF manifold coordinate `ξ_v` and `k` per-condition amplitudes `β_v`.
*   **Inputs:**
    *   `Y_bold`: `n x V` BOLD data matrix.
    *   `X_condition_list`: List of `k` matrices, where `X_condition_list[[c]]` is an `n x p` Toeplitz design matrix for condition `c` (convolving onsets with `p` delta functions).
    *   `Z_confounds`: (Optional) `n x q_confound` matrix of nuisance regressors.
    *   `B_reconstructor`: `p x m` matrix from Component 0.
    *   `lambda_gamma`: Ridge penalty for the `(km)x(km)` GLM solve.
    *   `h_ref_shape_canonical`: `p x 1` canonical HRF shape for identifiability.
    *   `orthogonal_approx_flag`: Boolean (default `FALSE`).
*   **Steps:**
    1.  **Confound Projection (if `Z_confounds` provided):**
        *   `Q_Z = qr.Q(qr(Z_confounds, LAPACK=TRUE))`.
        *   `Y_proj = Y_bold - Q_Z %*% tcrossprod(Q_Z, Y_bold)`.
        *   `X_condition_list_proj = lapply(X_condition_list, function(X) X - Q_Z %*% tcrossprod(Q_Z, X))`.
    2.  **Per-Condition Design in Manifold Basis `Z_list`:**
        *   `Z_list = lapply(X_condition_list_proj, function(Xc) Xc %*% B_reconstructor)`. Each `Z_list[[c]]` is `n x m`.
    3.  **Combined Design `X_tilde` and GLM Solve for `Gamma_coeffs`:**
        *   `X_tilde = do.call(cbind, Z_list)` (an `n x (km)` matrix).
        *   `XtX_tilde_reg = crossprod(X_tilde) + lambda_gamma * diag(k*m)`.
            *   If `orthogonal_approx_flag = TRUE`:
                `XtX_tilde_reg = Matrix::bdiag(lapply(Z_list, function(Zc) crossprod(Zc) + (lambda_gamma/k) * diag(m)))`. (*Note: distributing lambda might need care, or apply lambda to the full `diag(km)` before bdiag if that's simpler conceptually even if it means `lambda_gamma` has a slightly different effect*). A safer way is to construct the full `crossprod(X_tilde)` and then zero out off-diagonal blocks if orthogonal_approx is true, before adding the *full* `lambda_gamma * diag(k*m)`. For now, assume `XtX_tilde_reg` is the full version. The proposal's "code the full matrix; enable the block-diag shortcut only when the user explicitly requests it" implies `XtX_tilde_reg` is formed from `crossprod(X_tilde)` and then *optionally modified* if `orthogonal_approx`.
        *   `XtY_tilde = crossprod(X_tilde, Y_proj)`.
        *   `Gamma_coeffs = solve(XtX_tilde_reg, XtY_tilde)` (a `(km) x V` matrix).
    4.  **Extract `ξ_raw_allvox` and `β_raw_allvox` via SVD (per voxel):**
        *   Initialize `Xi_raw_allvox = matrix(0, m, V)` and `Beta_raw_allvox = matrix(0, k, V)`.
        *   Loop `vx` from 1 to `V`:
            *   `G_vx = matrix(Gamma_coeffs[, vx], nrow = m, ncol = k)`.
            *   `svd_G_vx = svd(G_vx)`.
            *   Handle near-zero singular values (set `xi_vx`, `beta_vx` to zero).
            *   `Xi_raw_allvox[,vx] = svd_G_vx$u[,1] * sqrt(svd_G_vx$d[1])`.
            *   `Beta_raw_allvox[,vx] = svd_G_vx$v[,1] * sqrt(svd_G_vx$d[1])`.
    5.  **Intrinsic Identifiability (on `ξ` and `β`):**
        *   `xi_ref_coord = MASS::ginv(B_reconstructor) %*% h_ref_shape_canonical`.
        *   For each voxel `vx`:
            *   `xi_vx = Xi_raw_allvox[,vx]`, `beta_vx = Beta_raw_allvox[,vx]`.
            *   `sgn = sign(sum(xi_vx * xi_ref_coord))`.
            *   `xi_vx = xi_vx * sgn`, `beta_vx = beta_vx * sgn`.
            *   `reconstructed_hrf = B_reconstructor %*% xi_vx`.
            *   `scl = 1 / pmax(max(abs(reconstructed_hrf)), .Machine$double.eps)`.
            *   `Xi_raw_allvox[,vx] = xi_vx * scl`.
            *   `Beta_raw_allvox[,vx] = beta_vx / scl`. (Ensure zeroed if scl was based on machine epsilon).

**Component 2: Spatial Smoothing of Manifold Coordinates**

*   **Objective:** Regularize HRF manifold coordinates spatially across the brain.
*   **Inputs:**
    *   `Xi_raw_allvox`: `m x V` matrix of (identifiability-constrained) manifold coordinates from Component 1.
    *   `voxel_graph_laplacian_Lsp`: `V x V` sparse graph Laplacian matrix.
    *   `lambda_spatial_smooth`: Spatial smoothing strength.
*   **Steps:**
    1.  Initialize `Xi_smoothed_allvox = matrix(0, m, V)`.
    2.  For each manifold dimension `j` from 1 to `m`:
        *   `Xi_smoothed_allvox[j,] = Matrix::solve(A = Diagonal(V) + lambda_spatial_smooth * voxel_graph_laplacian_Lsp, b = Xi_raw_allvox[j,])`.
    3.  **Output:** `Xi_smoothed_allvox` (`m x V`).

**Component 3: Trial-wise Amplitude Estimation (LSS using Smoothed Manifold HRFs)**

*   **Objective:** Estimate single-trial amplitudes `β_trial` using the spatially smoothed, voxel-specific HRFs.
*   **Inputs:**
    *   `Y_proj`: `n x V` confound-projected BOLD data.
    *   `X_trial_onset_list`: List of `T` matrices, where `X_trial_onset_list[[t]]` is an `n x p` Toeplitz design for trial `t`'s onset.
    *   `B_reconstructor`: `p x m` from Component 0.
    *   `Xi_smoothed_allvox`: `m x V` from Component 2.
    *   `A_lss_fixed`: `n x q_lss` matrix of fixed nuisance regressors for LSS (e.g., *original* `Z_confounds` if they were not related to task, or just an intercept).
    *   `p_lss_vec`: `n x 1` vector for Woodbury LSS, derived from `A_lss_fixed` (see section 2 of reviewed text: `MASS::ginv(A_lss_fixed)[intercept_row,]` or `rep(0,n)` if no intercept in `A_lss_fixed`).
    *   `lambda_ridge_Alss`: Ridge penalty for `(A_lss_fixedᵀA_lss_fixed + λI)⁻¹A_lss_fixedᵀ`.
*   **Steps (Voxel-Batched Woodbury - Mode B from LSS proposal):**
    1.  **Precompute `P_lss = (A_lss_fixedᵀA_lss_fixed + lambda_ridge_Alss * I)⁻¹A_lss_fixedᵀ` (`q_lss x n`).**
    2.  **Precompute Voxel-Specific HRF Shapes `H_shapes_allvox`:**
        *   `H_shapes_allvox = B_reconstructor %*% Xi_smoothed_allvox` (a `p x V` matrix).
    3.  **Optional (RAM permitting): Precompute all `R_t_allvox` matrices:**
        *   For trial `t = 1...T`: `R_t_allvox = X_trial_onset_list[[t]] %*% H_shapes_allvox` (an `n x V` matrix).
    4.  **Initialize `Beta_trial_allvox = matrix(0, T, V)`.**
    5.  **Voxel Loop (parallelizable):** For each voxel `v = 1...V`:
        a.  **Construct `C_v` (`n x T`):** `C_v[,t] = (if R_t precomputed) R_t_allvox[,v] else (X_trial_onset_list[[t]] %*% H_shapes_allvox[,v])`.
        b.  **Woodbury LSS Core (for current voxel `v`):**
            i.  `U_v = P_lss %*% C_v` (`q_lss x T`).
            ii. `V_regressors_v = C_v - A_lss_fixed %*% U_v` (`n x T`). These are `M_A c_t` for this voxel.
            iii. `pc_v_row = crossprod(p_lss_vec, C_v)` (`1 x T`).
            iv. `cv_v_row = colSums(V_regressors_v * V_regressors_v)` (`1 x T`), i.e., `||M_A c_t||²`.
            v.  `alpha_v_row = (1 - pc_v_row) / pmax(cv_v_row, .Machine$double.eps)`.
            vi. `S_effective_regressors_v = sweep(V_regressors_v, MARGIN = 2, STATS = alpha_v_row, FUN = "*")`.
            vii. `S_effective_regressors_v = sweep(S_effective_regressors_v, MARGIN = 1, STATS = p_lss_vec, FUN = "+")`.
        c.  **Betas for Voxel `v`: `Beta_trial_allvox[,v] = crossprod(S_effective_regressors_v, Y_proj[,v])` (`T x 1`).**

**6. Outputs of M-HRF-LSS Pipeline:**

*   `Xi_raw_allvox`: `m x V` raw manifold coordinates.
*   `Beta_condition_allvox`: `k x V` condition-level amplitudes.
*   `Xi_smoothed_allvox`: `m x V` spatially smoothed manifold coordinates.
*   `H_reconstructed_smoothed_allvox`: `p x V` reconstructed HRF shapes from `Xi_smoothed_allvox`.
*   `Beta_trial_allvox`: `T x V` trial-level LSS amplitudes.
*   Parameters used, `B_reconstructor`, etc.

**7. Engineering Considerations:**

*   **Modularity:** Each Component (0, 1, 2, 3) should be a distinct, testable function or set of functions.
*   **BLAS/LAPACK:** Rely on R's built-in highly optimized BLAS/LAPACK for matrix operations (`%*%`, `crossprod`, `solve`, `svd`, `qr`).
*   **Sparse Matrices:** Use the `Matrix` package for `voxel_graph_laplacian_Lsp` and its Cholesky solve.
*   **Memory:** Monitor memory usage, especially for `Gamma_coeffs` and `R_t_allvox`. Chunking of voxels for Component 3 might be needed for very large `V` if `R_t_allvox` cannot be held in RAM.
*   **Parallelization:** The voxel loops in Component 1 (SVD) and Component 3 (LSS kernel) are prime candidates for parallelization (e.g., `parallel::mclapply`, `future.apply`). Component 2 (spatial smoothing) involves `m` independent large sparse solves, also parallelizable.
*   **Error Handling & Input Validation:** Essential throughout.

This detailed proposal provides a comprehensive engineering blueprint for the M-HRF-LSS pipeline, incorporating all the excellent insights and refinements discussed.

Okay, let's break down the implementation of the "Manifold-Guided HRF Estimation and Trial-wise Deconvolution (M-HRF-LSS)" pipeline into a granular sprint plan with actionable tickets.

---

**Sprint Plan: M-HRF-LSS Pipeline Implementation (MVP)**

**Sprint Goal:** Deliver a functional R implementation of the M-HRF-LSS pipeline (Components 0-3), capable of processing moderately sized fMRI datasets, with a focus on algorithmic correctness and modularity. CPU-based, BLAS-optimized R code is the target.

**User Stories (Guiding the Sprint):**
*   As a methods developer, I need a robust implementation of each M-HRF-LSS component to enable research and validation.
*   As an fMRI researcher (early adopter), I want to run the M-HRF-LSS pipeline on a sample dataset to obtain spatially smoothed HRFs, condition betas, and trial betas.

**Tickets:**

**Epic 0: HRF Manifold Construction (Component 0 - Once per Study/Library)**

1.  **`MHRF-MANIFOLD-01`: Implement HRF Library Affinity & Markov Matrix Construction**
    *   **Task:** Create a function `calculate_manifold_affinity(L_library, k_local_nn_for_sigma = 7, sparse_if_N_gt = 5000, k_nn_for_W_sparse = 20)`
        *   Input: `L_library` (`p x N` HRF shapes).
        *   Compute pairwise Euclidean distances.
        *   Implement self-tuning bandwidth `σ_i`, `σ_j` (using `k_local_nn_for_sigma`).
        *   Compute affinity matrix `W_ij = exp(-dists_ij² / (σ_i * σ_j))`.
        *   (Optional) If `N > sparse_if_N_gt`, implement k-NN sparsification for `W`.
        *   Compute Markov Matrix `S = D_inv %*% W`.
    *   **DoD:** Function returns `S` (and optionally `W`, `D_inv`). Tested with small `L_library`.

2.  **`MHRF-MANIFOLD-02`: Implement Diffusion Map Eigendecomposition & Reconstructor `B`**
    *   **Task:** Create a function `get_manifold_basis_and_reconstructor(S_markov, L_library, m_manifold_dim)`
        *   Input: `S_markov` (`N x N`), `L_library` (`p x N`), `m_manifold_dim`.
        *   Use `RSpectra::eigs_sym(S_markov, k = m_manifold_dim + 1, which = "LM")` to get eigenvectors.
        *   Extract manifold coordinates `Φ = eig_S$vectors[, 2:(m_manifold_dim + 1)]` (`N x m`).
        *   Compute HRF reconstructor `B_reconstructor = L_library %*% Φ %*% solve(crossprod(Φ) + 1e-8 * diag(m_manifold_dim))` (`p x m`).
    *   **DoD:** Function returns `Φ` and `B_reconstructor`. Tested.

**Epic 1: Voxel-wise Manifold Fit (Component 1 - Per Subject/Run)**

3.  **`MHRF-VOXFIT-01`: Confound Projection Module**
    *   **Task:** Create helper `project_out_confounds(Y_data, X_list_data, Z_confounds_matrix)`
        *   Uses `qr.Q(qr(Z_confounds_matrix, LAPACK=TRUE))` for `Q_Z`.
        *   Returns projected `Y_proj` and `X_list_proj`.
    *   **DoD:** Module tested for correct projection.

4.  **`MHRF-VOXFIT-02`: Per-Condition Design in Manifold Basis (`Z_list`)**
    *   **Task:** Create helper `transform_designs_to_manifold_basis(X_condition_list_proj, B_reconstructor)`
        *   Returns `Z_list` where `Z_list[[c]] = X_condition_list_proj[[c]] %*% B_reconstructor`.
    *   **DoD:** Module tested.

5.  **`MHRF-VOXFIT-03`: Main GLM Solve for `Gamma_coeffs`**
    *   **Task:** Implement the GLM solve part of Component 1.
        *   Input: `Z_list`, `Y_proj`, `lambda_gamma`, `orthogonal_approx_flag`.
        *   Form `X_tilde = do.call(cbind, Z_list)`.
        *   Form `XtX_tilde_reg = crossprod(X_tilde) + lambda_gamma * diag(k*m)`.
            *   (Later) Implement `orthogonal_approx_flag` modification to `XtX_tilde_reg`.
        *   Solve for `Gamma_coeffs`.
    *   **DoD:** `Gamma_coeffs` (`(km) x V`) computed correctly.

6.  **`MHRF-VOXFIT-04`: SVD-based Extraction of `ξ_raw` and `β_raw`**
    *   **Task:** Implement the per-voxel SVD loop from Component 1, Step 4.
        *   Input: `Gamma_coeffs`, `m_manifold_dim`, `k_conditions`.
        *   Output: `Xi_raw_allvox` (`m x V`), `Beta_raw_allvox` (`k x V`).
        *   Include robustness for near-zero singular values.
    *   **DoD:** Raw `Xi` and `Beta` extracted.

7.  **`MHRF-VOXFIT-05`: Intrinsic Identifiability for `ξ` and `β`**
    *   **Task:** Implement Component 1, Step 5.
        *   Input: `Xi_raw_allvox`, `Beta_raw_allvox`, `B_reconstructor`, `h_ref_shape_canonical`.
        *   Compute `xi_ref_coord = MASS::ginv(B_reconstructor) %*% h_ref_shape_canonical`.
        *   Apply sign and scale adjustments to `Xi_raw_allvox` and `Beta_raw_allvox` per voxel.
    *   **DoD:** Identifiability constraints correctly applied. Output `Xi_ident_allvox`, `Beta_ident_allvox`.

**Epic 2: Spatial Smoothing (Component 2 - Per Subject/Run)**

8.  **`MHRF-SP SMOOTH-01`: Graph Laplacian Construction**
    *   **Task:** Create helper `make_voxel_graph_laplacian(voxel_coordinates_matrix, num_neighbors = 6)`
        *   Input: `V x 3` matrix of voxel coordinates.
        *   Constructs a sparse `V x V` graph Laplacian (e.g., based on k-NN or distance threshold).
    *   **DoD:** Sparse `L_sp` matrix generated.

9.  **`MHRF-SP SMOOTH-02`: Apply Spatial Smoothing to `ξ` Coordinates**
    *   **Task:** Implement Component 2, Step 2.
        *   Input: `Xi_ident_allvox` (from `MHRF-VOXFIT-05`), `voxel_graph_laplacian_Lsp`, `lambda_spatial_smooth`.
        *   Loop `m` times, solving the sparse system `(I_V + λ_sp L_sp) ξ_j_smooth = ξ_j_ident`.
    *   **DoD:** `Xi_smoothed_allvox` (`m x V`) computed.

**Epic 3: Trial-wise LSS (Component 3 - Per Subject/Run)**

10. **`MHRF-LSS-01`: LSS Fixed Regressor Precomputation (`P_lss`, `p_lss_vec`)**
    *   **Task:** Create helper `prepare_lss_fixed_components(A_lss_fixed, intercept_column_in_Alss = NULL, lambda_ridge_Alss = 0.1)`
        *   Computes `P_lss = (A_lss_fixedᵀA_lss_fixed + lambda_ridge_Alss * I)⁻¹A_lss_fixedᵀ`.
        *   Computes `p_lss_vec` based on `A_lss_fixed` and `intercept_column_in_Alss`.
    *   **DoD:** `P_lss` and `p_lss_vec` returned.

11. **`MHRF-LSS-02`: Voxel-Specific HRF Shape Reconstruction**
    *   **Task:** Helper `reconstruct_hrf_shapes(B_reconstructor, Xi_smoothed_allvox)`
        *   Returns `H_shapes_allvox = B_reconstructor %*% Xi_smoothed_allvox` (`p x V`).
    *   **DoD:** Voxel-specific HRF shapes computed.

12. **`MHRF-LSS-03`: Implement Core Woodbury LSS Kernel (Per Voxel)**
    *   **Task:** Create internal function `run_lss_for_voxel(Y_proj_vx, X_trial_onset_list, H_shape_vx, A_lss_fixed, P_lss, p_lss_vec)`
        *   Implements Component 3, Steps 5a-5c (forming `C_vx`, then Woodbury steps for `S_vx`, then `beta_trial_vx`).
    *   **DoD:** Function returns `T x 1` trial betas for one voxel. Tested.

13. **`MHRF-LSS-04`: Main LSS Loop and (Optional) `R_t` Precomputation**
    *   **Task:** Implement the main loop of Component 3.
        *   Optionally precompute all `R_t_allvox` if memory strategy allows.
        *   Loop over voxels, calling `run_lss_for_voxel`.
    *   **DoD:** `Beta_trial_allvox` (`T x V`) computed.

**Epic 4: Pipeline Orchestration, Output, and Documentation (MVP)**

14. **`MHRF-PIPELINE-01`: Top-Level Orchestration Function `run_mhrf_lss_pipeline()`**
    *   **Task:** Create the main function that calls Components 0 (if needed, or loads precomputed B), 1, 2, and 3 in sequence.
    *   Manages inputs and outputs between components.
    *   **DoD:** Full pipeline runnable end-to-end on sample data.

15. **`MHRF-OUTPUT-02`: Define and Implement M-HRF-LSS Output Object**
    *   **Task:** Create an S3/S4 class to store all key outputs (`Xi_raw`, `Beta_condition`, `Xi_smoothed`, `H_reconstructed`, `Beta_trial`, parameters, `B_reconstructor`).
    *   Implement basic `print()`, `summary()`, and `plot()` (e.g., plot mean reconstructed HRF, map of one `ξ` coordinate).
    *   **DoD:** Structured output available with basic utility methods.

16. **`MHRF-DOC-03`: Initial User Documentation and Example Script**
    *   **Task:** Draft R documentation for the main pipeline function and output object.
    *   Create a simple example script using synthetic or public data to demonstrate usage.
    *   **DoD:** Basic documentation and a runnable example are available.

**Sprint Review Focus:**
*   Successful end-to-end execution of the M-HRF-LSS pipeline on a test dataset.
*   Plausibility of estimated HRF shapes, condition betas, and trial betas.
*   Modularity and correctness of individual components.
*   Initial performance benchmarks (CPU-based R).

**Post-MVP Considerations (Future Sprints):**
*   Parallelization of voxel loops (`parallel`, `future.apply`, or `RcppParallel`).
*   `RcppArmadillo` optimization for computationally intensive loops (SVD, LSS kernel).
*   Advanced output visualization.
*   Robust parameter selection helpers (e.g., for `m_manifold_dim`, `lambda_spatial_smooth`).
*   Implementation of the "Bonus Refinement" (joint spatial and manifold smoothness prior).

This granular sprint plan breaks the complex M-HRF-LSS pipeline into manageable, testable tickets, focusing on delivering a functional MVP first.