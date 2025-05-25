Okay, let's flesh this out into a formal proposal for a "Next-Generation CPU-Based LSS Kernel," focusing on BLAS-optimized operations and the two identified operating modes for HRF handling.

---

**Proposal: High-Performance CPU-Based Least Squares Separate (LSS) Kernel**

**1. Introduction & Motivation**

Least Squares Separate (LSS) is a widely used method for estimating fMRI BOLD responses to individual trials by fitting a separate GLM for each trial where the target trial is modeled uniquely and other trials are treated as nuisance. While effective, traditional LSS implementations can be computationally intensive for large datasets (many trials, many voxels). AFNI's `3dLSS` is a notable optimized C implementation. This proposal outlines an even faster CPU-based LSS kernel that achieves significant speedups by:
1.  Eliminating the explicit formation or application of `n x n` projection matrices.
2.  Leveraging the Woodbury matrix identity for efficient residualization.
3.  Maximizing the use of highly optimized BLAS (Basic Linear Algebra Subprograms) level-2 and level-3 operations, vectorizing computations across trials.

This kernel is designed to be integrated with fMRI analysis pipelines, particularly those using data-driven Hemodynamic Response Functions (HRFs) like those from CF-ALS.

**2. Core Algebraic Optimizations**

For a standard LSS setup, for each target trial `c` (an `n x 1` regressor), we want to estimate its beta coefficient in a model that also includes a set of common or nuisance regressors `A` (an `n x m` matrix, where `m` is typically small, e.g., other conditions, motion, baseline).

The traditional approach involves projecting `c` and the data `y` onto the space orthogonal to `A` using the projector `M_A = I - A A⁺`, where `A⁺ = (AᵀA)⁻¹Aᵀ` is the pseudoinverse.
*   `c_perp = M_A c`
*   `y_perp = M_A y`
*   `β_c = (c_perpᵀ c_perp)⁻¹ c_perpᵀ y_perp`

The key optimizations are:

*   **Woodbury for Residualization:** Instead of forming `M_A`, compute `M_A c = c - A(A⁺c)`. This reduces the per-trial cost for residualizing `c` from `O(n²)` or `O(nm)` (if `M_A` is applied naively) to `O(mn)` more efficiently.
*   **Vectorization Across Trials:** If `C` is an `n x T` matrix of all trial-specific regressors:
    *   `U = A⁺C` (all `A⁺c_t` in one go)
    *   `V = C - AU` (all `M_Ac_t` in one go, stored in `V`). `V` is `n x T`.
*   **Efficient Beta Calculation:** The LSS beta for trial `t` can be expressed as `β_t = s_tᵀy`, where `s_t` is an "effective regressor" or "influence vector" for trial `t`.
    *   Based on common LSS derivations (e.g., similar to AFNI's approach), `s_t` can be constructed efficiently from `V_t = (M_A c_t)`, `c_t`, and specific rows of `A⁺` (denoted `p` below).
    *   The denominator `c_tᵀM_A c_t` can be stably computed as `‖M_A c_t‖² = V_tᵀV_t`.
    *   The specific formula for `s_t` used in this proposal is:
        `s_t = p_vec + α_t V_t`, where `V_t` is the `t`-th column of `V`, `p_vec` is an `n x 1` vector derived from `A⁺` (e.g., the row of `A⁺` corresponding to a global "all-other trials" regressor, if used, or another relevant vector for adjusting specific common effects), and `α_t = (1 - p_vecᵀc_t) / (V_tᵀV_t)`.

**3. Algorithmic Steps (CPU Optimized)**

Let:
*   `Y`: `n x Vvox` BOLD data matrix.
*   `A`: `n x m` matrix of common/nuisance regressors (fixed across trials and voxels).
*   `C`: `n x T` matrix of trial-specific regressors (can be common for all voxels or voxel-specific).
*   `p_vec`: An `n x 1` vector derived from `A` or `A⁺` as per the specific LSS formulation (e.g., a specific column of `A` or related to a specific row of `A⁺`). This adjusts the baseline or other common effects.
*   `lambda_ridge_A`: Ridge penalty for `(AᵀA + λI)⁻¹Aᵀ`.

**Precomputation (once per analysis):**
1.  `AtA_reg = AᵀA + lambda_ridge_A * I_m`
2.  `P = (AtA_reg)⁻¹Aᵀ` (This is `A⁺` effectively, `m x n`).
3.  (Optional) If `p_vec` is derived from `P` (e.g., a specific row of `P` transposed), extract it.

**LSS Kernel (Two Operating Modes):**

**Mode A: Shared Trial Regressors `C` (e.g., using a global/ROI-average HRF)**

1.  **Construct `C` (`n x T`):** Each column `c_t` is the regressor for trial `t` (e.g., `X_onset_t @ h_shared`).
2.  **Compute `U = PC`:** `(m x n) %*% (n x T) -> m x T`. (BLAS `gemm`)
3.  **Compute `V = C - AU`:** `(n x T) - (n x m) %*% (m x T) -> n x T`. (BLAS `gemm` for `AU`, then matrix subtraction). `V` contains all `M_A c_t`.
4.  **Compute `pc_row = p_vecᵀC`:** `(1 x n) %*% (n x T) -> 1 x T`. (BLAS `gemv` or `gemm`). This contains `p_vecᵀc_t` for each trial.
5.  **Compute `cv_row = colSums(V * V)`:** `(1 x T)`. This contains `‖M_A c_t‖²` for each trial. Numerically stable.
6.  **Compute `alpha_row = (1 - pc_row) / cv_row`:** `(1 x T)`. Element-wise. Handle potential division by zero if `cv_row` is zero (trial regressor perfectly collinear with `A`).
7.  **Construct `S` (`n x T`):** For each trial `t` (column `t`): `S_t = p_vec + alpha_row[t] * V_t`. This can be done efficiently with `sweep` or `axpy`-like operations per column.
8.  **Compute All Betas `B = SᵀY`:** `(T x n) %*% (n x Vvox) -> T x Vvox`. (BLAS `gemm`).

**Mode B: Voxel-Specific Trial Regressors `C_v` (e.g., using voxel-specific `h_v` from CF-ALS)**

1.  **Precomputation for Trial Onset Bases:** Let `X_onset_list` be a list of `T` matrices, where `X_onset_list[[t]]` is `n x d` (e.g., FIR basis for trial `t`'s onset).
2.  **Precompute All Trial Regressors for All Voxels (Optional, if RAM permits):**
    *   For each trial `t`: `R_t_allvox = X_onset_list[[t]] @ H_allvoxels` where `H_allvoxels` is `d x Vvox`. `R_t_allvox` is `n x Vvox`. Store `T` such matrices.
    *   This step involves `T` BLAS `gemm` operations.
3.  **Voxel Loop (can be parallelized, e.g., OpenMP in C++ or `mclapply` in R):**
    For each voxel `v = 1...Vvox`:
    a.  **Construct `C_v` (`n x T`):**
        *   If `R_t_allvox` are precomputed: `C_v[,t] = R_t_allvox[,v]`. (Memory gather)
        *   Else (compute on-the-fly): `C_v[,t] = X_onset_list[[t]] @ H_allvoxels[,v]`. (BLAS `gemv`)
    b.  **Compute `U_v = PC_v`:** `(m x n) %*% (n x T) -> m x T`. (BLAS `gemm`)
    c.  **Compute `V_v = C_v - AU_v`:** `(n x T) - (n x m) %*% (m x T) -> n x T`. (BLAS `gemm` for `AU_v`, then subtraction).
    d.  **Compute `pc_v_row = p_vecᵀC_v`:** `(1 x n) %*% (n x T) -> 1 x T`. (BLAS `gemv` or `gemm`).
    e.  **Compute `cv_v_row = colSums(V_v * V_v)`:** `(1 x T)`.
    f.  **Compute `alpha_v_row = (1 - pc_v_row) / cv_v_row`:** `(1 x T)`.
    g.  **Construct `S_v` (`n x T`):** `S_v_t = p_vec + alpha_v_row[t] * V_v_t`.
    h.  **Compute Betas for Voxel `v`: `B_v_col = S_vᵀY[,v]`:** `(T x n) %*% (n x 1) -> T x 1`. (BLAS `gemv`). Store `B_v_col` as `B_allvoxels[,v]`.

**4. Expected Performance & Complexity (CPU-Based)**

*   **Mode A (Shared `C`):**
    *   Dominant costs: `PC` (`O(mnT)`), `AU` (`O(nmT)`), `SᵀY` (`O(nTVvox)`).
    *   Total roughly `O(mnT + nTVvox)`. If `m << Vvox`, then `O(nTVvox)`.
*   **Mode B (Voxel-Specific `C_v`):**
    *   Precomputing `R_t_allvox`: `O(TndVvox)`.
    *   Per-voxel kernel (Steps 3a-3h): `O(nmT + nT)`. Total `Vvox * (nmT)`.
    *   Overall: `O(TndVvox + Vvox nmT) = O(TVvox n(d+m))`.
*   Both modes significantly reduce computations compared to naive LSS (`O(Vvox * T * n²)` or `O(Vvox * T * nm)` if `M_A` is applied efficiently but still per trial-GLM).

**5. Integration with CF-ALS Data-Driven HRFs**

*   The CF-ALS method provides `H_allvoxels` (`d x Vvox`), the matrix of HRF basis coefficients for all voxels.
*   This `H_allvoxels` is used directly in Mode B (Step 2 or 3a) to construct voxel-specific trial regressors `C_v`.
*   For Mode A, an `h_shared` (e.g., average of `H_allvoxels` over an ROI or whole brain) would be used to construct a single `C`.

**6. Implementation Considerations (R with BLAS / C++ with BLAS)**

*   **R:** Base R matrix operations (`%*%`, `crossprod`, `solve`, `sweep`, `colSums`) are already linked to BLAS/LAPACK. The structure above is directly translatable. Parallelization of the voxel loop in Mode B can be done with `parallel::mclapply` or `foreach`.
*   **C++ with BLAS (e.g., via RcppArmadillo or standalone):**
    *   Direct calls to BLAS functions (e.g., `sgemm`, `sgemv`, `saxpy`, `sdot`) for maximal performance.
    *   OpenMP for parallelizing the voxel loop in Mode B.
    *   Armadillo library simplifies C++ matrix algebra and links to BLAS.

**7. Prewhitening and Other Extensions**

*   **Prewhitening:** If a whitening matrix `W` (e.g., from `chol(Sigma_inv)`) is available, apply it to `Y`, `A`, and each `c_t` (or `X_onset_list[[t]]`) *before* starting the LSS kernel. The algebra remains the same for the whitened data.
*   **Ridge Penalties:** As noted, `P_λ = (AᵀA + λI)⁻¹Aᵀ` can be used. Additional ridge penalties could be incorporated into the LSS GLM itself if needed, but the current proposal focuses on speeding up standard LSS.

---

**Conclusion:**

This proposal outlines a highly efficient CPU-based LSS kernel leveraging established linear algebra optimizations and BLAS routines. It offers two modes for practical application, especially when combined with voxel-specific HRFs from methods like CF-ALS. Mode A provides maximum throughput when a shared HRF is appropriate, while Mode B offers substantial speedups for fully voxel-specific HRF LSS analyses. This approach promises to make LSS significantly more tractable for large-scale fMRI datasets on standard CPU hardware.