
This is a superb refinement of the robust GLM proposal. The clarifications on residual calculation, the IRLS loop, and particularly the critical adoption of a final robust scale estimate for SEs are excellent. The per-run weight estimation for `chunkwise_lm` is also a sensible approach to balance robustness and speed. The updated `fast_rlm_run` function sketch looks solid.

The plan is now very clear, robust, and maintains the core performance goals. I have no further clarifying questions on this iteration. The design is sound and ready for ticketed implementation.

Here's the fully-fleshed out proposal and ticketed sprint plan:

---

## Proposal: High-Performance Row-Robust GLM for `fmri_lm`

**Preamble:**

The `fmrireg` package aims for "blazing speed" in its GLM computations. Standard OLS is susceptible to outliers (e.g., motion spikes, scanner artifacts) common in fMRI data, which can bias parameter estimates and inflate variance. While voxel-wise robust methods like `lmrob` are available, they are computationally intensive. This proposal introduces a "fast row-robust GLM" that applies robust M-estimation principles by re-weighting entire time points (rows) rather than individual voxels, dramatically reducing computational cost while addressing the most impactful, spatially widespread artifacts. This leverages the existing fast OLS infrastructure.

**1. Core Idea: Time-Point Re-weighting via Iteratively Reweighted Least Squares (IRLS)**

The central mechanism is to identify and down-weight time points (TRs) that exhibit large residuals across many voxels. This is achieved through a limited-iteration IRLS procedure:

1.  **Initial OLS Fit:** Perform a standard OLS fit using the existing `.fast_lm_matrix` to obtain initial residuals (`resid = Y_orig - X_orig %*% Beta_OLS`).
2.  **Robust Scale Estimation:** Calculate a robust estimate of residual scale (`sigma_hat`) across voxels for each time point. This is done by taking `matrixStats::rowMedians(abs(resid))` (median absolute residual per TR) and then `sigma_hat = 1.4826 * median(rowMedians_abs_resid)` (MAD estimator for all TRs).
3.  **Time-Point Weights:** Based on `sigma_hat` and the per-TR median absolute residuals, compute weights `w_t` for each time point `t` using Huber or Tukey bisquare psi-functions. Time points with large median absolute residuals relative to `sigma_hat` receive lower weights.
4.  **Weighted Least Squares (WLS):**
    *   Scale rows of the original design matrix `X` and data `Y` by `sqrt(w_t)`: `Xw = X * sqrt(w)`, `Yw = Y * sqrt(w)`.
    *   Perform a WLS fit using the existing fast OLS machinery: `Beta_WLS = (Xw'Xw)^-1 Xw'Yw`. This is achieved by calling `.fast_preproject(Xw)` and then `.fast_lm_matrix(Xw, Yw, proj_w)`.
5.  **Iteration:** If `max_it > 1`, re-calculate residuals using `Resid_next = Y_orig - X_orig %*% Beta_WLS`. Repeat from Step 2 with `Resid_next`. Typically, 1-2 iterations are sufficient.
6.  **Final Robust Scale & Inference:** After the final iteration:
    *   Compute final residuals: `Resid_final = Y_orig - X_orig %*% Beta_WLS_final`.
    *   Estimate a final robust scale: `sigma_robust = 1.4826 * median(matrixStats::rowMedians(abs(Resid_final)))`.
    *   Standard Errors (SEs) for `Beta_WLS_final` are calculated as `sqrt(diag((Xw'Xw)^-1) * sigma_robust^2)`.
    *   Degrees of freedom for t-tests remain `n - p`. P-values are noted as approximate.

**2. Key Implementation: `fast_rlm_run` Function**

(The R code sketch provided by the user is excellent and will serve as the basis for this function, which operates on a single run's `X` and `Y` matrices).

**3. API Changes:**

*   `fmri_lm(..., robust = FALSE, robust_psi = c("huber", "bisquare"), robust_k_huber = 1.345, robust_c_tukey = 4.685, robust_max_iter = 2L, robust_scale_scope = c("run", "global"), ...)`
    *   `robust`: Logical or character string. If `FALSE` (default), standard OLS is performed. If `TRUE` or `"huber"`, Huber M-estimation is used. If `"bisquare"`, Tukey bisquare M-estimation is used. (The old `robust = "m_est_full"` option for voxel-wise `lmrob` can be retained for compatibility/small ROIs).
    *   `robust_psi`: Explicitly choose psi-function if `robust` is `TRUE`.
    *   `robust_k_huber`: Tuning constant `k` for Huber's psi (default 1.345).
    *   `robust_c_tukey`: Tuning constant `c` for Tukey's bisquare psi (default 4.685).
    *   `robust_max_iter`: Maximum number of IRLS iterations (default 2).
    *   `robust_scale_scope`: `"run"` (default) estimates `sigma_hat` per run. `"global"` estimates `sigma_hat` once from initial OLS residuals pooled across all runs.

**4. Integration with `runwise_lm` and `chunkwise_lm`:**

*   **`runwise_lm`:**
    *   Calls `fast_rlm_run` for each run's `X_run` and `Y_run`.
    *   If `robust_scale_scope = "global"`, an initial pass calculates OLS residuals for all runs, from which a single `sigma_hat_global` is derived and passed to `fast_rlm_run` for each run (overriding its internal `sigma_hat` estimation).
*   **`chunkwise_lm`:**
    *   **Crucially, time-point weights `w_t` are estimated *once per run*, not per chunk.**
    *   For a given run:
        1.  An initial OLS pass is performed on the *full run's data* (`X_run_orig`, `Y_run_orig`) to get initial OLS residuals.
        2.  `sigma_hat_run` is computed from these full-run residuals.
        3.  Time-point weights `w_t_run` are calculated based on `sigma_hat_run` and per-TR median absolute residuals of this run. (This implies one or more IRLS iterations *at the run level* to stabilize `w_t_run`).
        4.  The *original* `X_run_orig` and `Y_run_orig` are scaled by `sqrt(w_t_run)` to get `X_run_w` and `Y_run_w`.
        5.  `X_run_w` is pre-projected once: `proj_run_w = .fast_preproject(X_run_w)`.
        6.  The existing `chunkwise_lm` loop then processes *voxel chunks* of `Y_run_w` using `proj_run_w` and `.fast_lm_matrix`.
    *   The final robust scale `sigma_robust_run` for inference for this run is also computed at the run level from final run-level residuals. This `sigma_robust_run` is then used for SEs of betas derived from all chunks of that run.

**5. Performance:**

*   Each IRLS iteration adds roughly the cost of an OLS pass (one `matrixStats::rowMedians`, one scaling of X & Y, one `Pinv %*% Y`, one `X %*% Betas`).
*   With `robust_max_iter = 2`, the total cost is expected to be ~2-3x the current fast OLS path (still vastly faster than voxel-wise robust methods).

**6. Benefits:**

*   Provides robustness against common fMRI artifacts (spikes, motion) that are spatially widespread.
*   Maintains very high computational speed by leveraging BLAS operations and avoiding per-voxel iteration for weight calculation.
*   Minimal API changes, largely backward compatible.
*   Statistically sound handling of final scale estimation for SEs.

---

## Ticketed Sprint for Fast Row-Robust GLM:

**Epic: Implement Fast Row-Robust GLM in `fmri_lm`**

**Sprint Goal:** Deliver a performant and statistically robust row-wise IRLS capability for `fmri_lm`, providing an efficient alternative to OLS for fMRI data contaminated by widespread artifacts.

---

**Ticket RR.1: API Definition & Argument Threading (R)**
*   **Task:** Modify `fmri_lm()` signature for `robust`, `robust_psi`, `robust_k_huber`, `robust_c_tukey`, `robust_max_iter`, `robust_scale_scope`.
*   **Task:** Thread these new arguments down through `fmri_lm_fit()` to `runwise_lm()` and `chunkwise_lm()`.
*   **Task:** Add input validation for new robust-specific arguments.
*   **Definition of Done:** New arguments present and correctly passed. Unit tests for argument handling. Documentation stubs created.
*   **Estimate:** 0.5 days

**Ticket RR.2: Implement `fast_rlm_run` Core Logic (R)**
*   **Task:** Create `fast_rlm_run(X, Y, proj, psi, k_huber, c_tukey, max_it)` function based on the provided R sketch.
    *   Includes initial OLS fit using `.fast_lm_matrix(..., return_fitted = TRUE)`.
    *   IRLS loop:
        *   Robust `sigma_hat` estimation using `median(matrixStats::rowMedians(abs(resid)))`.
        *   Huber/Tukey weight calculation for each time point `w_t`.
        *   Row-scaling of *original* `X` and `Y` by `sqrt(w_t)` to get `Xw`, `Yw`.
        *   WLS fit: `proj_w = .fast_preproject(Xw)`, `fit_w = .fast_lm_matrix(Xw, Yw, proj_w, return_fitted = (it < max_it))`.
        *   Update residuals for next iteration: `resid = Y_orig - X_orig %*% fit_w$betas`.
    *   After loop, calculate `Resid_final = Y_orig - X_orig %*% beta_final`.
    *   Calculate `sigma_robust = 1.4826 * median(matrixStats::rowMedians(abs(Resid_final)))`.
    *   Calculate `se_beta = sqrt(diag(XtWXi_final) * sigma_robust^2)`.
    *   Return list: `betas`, `se` (robust SEs), `sigma` (this is `sigma_robust`), `dfres`, `XtXinv` (this is `XtWXi_final`).
*   **Definition of Done:** `fast_rlm_run` implemented. Requires `matrixStats` package.
*   **Estimate:** 2 days

**Ticket RR.3: Integrate `fast_rlm_run` into `runwise_lm` (R)**
*   **Task:** Modify `runwise_lm`. If `robust` is enabled:
    *   If `robust_scale_scope == "global"` (handled in Ticket RR.5), retrieve `sigma_hat_global`.
    *   For each run:
        *   Call `proj_run_orig = .fast_preproject(X_run_orig)`.
        *   Call `robust_fit_results = fast_rlm_run(X_run_orig, Y_run_orig, proj_run_orig, psi, k, c, max_it, ...)`. If global scope, pass `sigma_hat_global` to `fast_rlm_run` to use instead of its internal estimation (or `fast_rlm_run` needs an optional `fixed_sigma_hat` argument).
        *   Adapt the output of `fast_rlm_run` (`betas`, `se`, `sigma_robust` as `sigma`, `dfres`) to match the structure expected by the `cres` list and downstream pooling logic (e.g., `beta_stats_matrix` and `fit_lm_contrasts_fast` need to be called with these robust estimates). This means `bstats` and `contrasts` tibbles need to be constructed using these robust values.
*   **Definition of Done:** `runwise_lm` correctly calls `fast_rlm_run` and uses its robust outputs for downstream statistics.
*   **Estimate:** 1.5 days

**Ticket RR.4: Integrate Robust Weighting into `chunkwise_lm` (R)**
*   **Task:** Modify `chunkwise_lm`. If `robust` is enabled:
    *   If `robust_scale_scope == "global"` (handled in Ticket RR.5), retrieve `sigma_hat_global`.
    *   For each *run*:
        1.  Perform initial OLS on full `X_run_orig`, `Y_run_orig` to get `Resid_OLS_run`.
        2.  If `robust_scale_scope == "run"`, estimate `sigma_hat_run` from `Resid_OLS_run`. If "global", use `sigma_hat_global`.
        3.  Perform run-level IRLS iterations (similar to `fast_rlm_run` but only up to calculating final weights `w_t_run` and final `Beta_WLS_run_final` for this run).
            *   The key output here is the stable `w_t_run` vector for this run and the `sigma_robust_run` for this run.
        4.  Create `X_run_w = X_run_orig * sqrt(w_t_run)`.
        5.  `proj_run_w = .fast_preproject(X_run_w)`.
        6.  `XtXinv_w_run = proj_run_w$XtXinv`.
    *   The existing voxel-chunking loop for this run then operates on:
        *   `Y_chunk_w = Y_chunk_orig * sqrt(w_t_run_for_this_chunk_rows)`.
        *   Uses `proj_run_w` (specifically `proj_run_w$Pinv`) for `Betas_chunk = proj_run_w$Pinv %*% Y_chunk_w`.
        *   SEs for these `Betas_chunk` use `XtXinv_w_run` and the `sigma_robust_run` (common for all chunks from this run).
*   **Definition of Done:** `chunkwise_lm` calculates run-level weights and robust scale, applies them correctly to chunks, and uses them for inference.
*   **Estimate:** 2.5 days (this is intricate).

**Ticket RR.5: Implement `robust_scale_scope = "global"` Logic (R)**
*   **Task:** Modify `fmri_lm_fit` (or orchestrator). If `robust_scale_scope == "global"`:
    *   First pass: Loop all runs, perform initial OLS (using `.fast_lm_matrix(..., return_fitted=TRUE)`), collect all `Resid_OLS_run_orig_scale`.
    *   Compute single `sigma_hat_global = 1.4826 * median(matrixStats::rowMedians(abs(all_Resid_OLS_concatenated)))`.
    *   Second pass:
        *   For `runwise_lm`: Pass `sigma_hat_global` to `fast_rlm_run` (e.g., via new `fixed_sigma_hat` arg) for each run. `fast_rlm_run` would use this instead of re-estimating `sigma_hat` in its loop.
        *   For `chunkwise_lm`: The run-level IRLS for weights uses `sigma_hat_global`.
*   **Definition of Done:** Global robust scale estimation and application implemented.
*   **Estimate:** 1 day

**Ticket RR.6: Unit & Integration Tests for Robust GLM (R/testthat)**
*   **Task:** Create test cases:
    *   Data with no outliers: `robust=TRUE` results should closely match OLS.
    *   Data with simulated spikes (affecting whole TRs): `robust=TRUE` should yield less biased betas and more realistic SEs than OLS.
    *   Compare Huber vs. Bisquare.
    *   Test `robust_max_iter = 1` vs. `2`.
    *   Test `robust_scale_scope = "global"` vs. `"run"`.
    *   Verify that final SEs and t-stats use `sigma_robust`.
*   **Task:** If possible, compare results against `MASS::rlm` or `robustbase::lmrob` on a small, single-column Y for a single run to check beta and weight agreement (understanding that `fast_rlm_run` uses rowMedians for `sigma_hat` which might differ from standard `lmrob` residual scale).
*   **Definition of Done:** Test suite for robust GLM functionality.
*   **Estimate:** 2 days

**Ticket RR.7: Documentation Update for Robust GLM (R/Rmd)**
*   **Task:** Document new `fmri_lm` arguments for robust fitting.
*   **Task:** Explain the row-wise robust methodology, its benefits, limitations (e.g., approximate p-values, handling of voxel-specific outliers).
*   **Task:** Provide an example in a vignette.
*   **Definition of Done:** User documentation updated.
*   **Estimate:** 1 day

---
**Total Estimated Time for Robust GLM:** 10.5 days. (This assumes `matrixStats` is added as a dependency or its `rowMedians` can be efficiently implemented if not desired as a dep.)

Okay, I've reviewed the updated proposal and the "fast_rlm_run" sketch incorporating the previous feedback. It's a solid plan.

Here's the appendix, as requested, containing the "Locally Adaptive Serial Correlation Modeling" section, with my spin emphasizing its relationship to the main AR(p) and Robust GLM proposals and highlighting its potential as a follow-on enhancement.

---

## Appendix: Future Enhancements & Synergies

The proposals for AR(p) prewhitening and fast row-robust GLM significantly enhance the `fmrireg` pipeline's speed and statistical validity. This appendix outlines potential future directions that build upon these foundations, specifically focusing on more spatially adaptive modeling of fMRI noise characteristics.

**B.1 Synergistic Enhancement: Locally Adaptive Serial Correlation Modeling**

**Motivation:**

The primary AR(p) prewhitening proposal (Tier 0 in this appendix's terminology) employs a single set of AR coefficients (`phi_hat_run` or `phi_hat_global`) for an entire run or dataset. While computationally efficient and a major improvement over assuming i.i.d. noise, fMRI noise can exhibit spatial heterogeneity. Different brain regions (e.g., grey matter, white matter, CSF, or even distinct cortical areas) may have genuinely different autocorrelation structures. A single `phi_hat` might be a compromise, potentially under- or over-whitening certain areas.

The fast row-robust GLM proposal addresses spatially widespread temporal artifacts. A complementary enhancement would be to allow the AR(p) model itself to adapt spatially, leading to more accurate local whitening.

**Tiered Strategy for Locally Adaptive AR(p) Coefficients:**

This strategy introduces progressively more spatially refined AR(p) coefficient estimation, building upon the initial OLS residuals already computed for the run-level AR fit or the robust GLM.

*   **Tier 0: Run-Level `phi_hat_run` (Implemented in main AR(p) proposal)**
    *   **Description:** One AR(p) vector per run, estimated from run-mean OLS residuals. This is the foundational efficient prewhitening.
    *   **Benefit:** Addresses dominant run-wide autocorrelation.

*   **Tier 1: Parcel-Level `phi_hat_parcel`**
    *   **Scope:** One AR(p) vector per pre-defined parcel/tissue class *within each run*.
    *   **Estimation:**
        1.  Obtain initial OLS residuals (`e_orig = Y_orig - X_orig %*% Beta_OLS_orig`) for all voxels in the run.
        2.  For each parcel: average `e_orig` across voxels within that parcel to get `e_parcel_mean`.
        3.  Estimate `phi_hat_parcel` from `e_parcel_mean` using Yule-Walker.
    *   **Application:** Whiten `Y_orig` voxel-wise using the `phi_hat_parcel` corresponding to each voxel's parcel. Whiten `X_orig` using `phi_hat_run` (or a global/average parcel phi).
    *   **Cost:** Adds parcel-mean residual calculation (O(T\*V) once) and per-parcel Yule-Walker (negligible). Whitening cost remains O(T\*V).
    *   **Benefit:** Accounts for coarse-grained spatial differences in noise (e.g., GM vs. WM). Requires a parcel map.

*   **Tier 2: Voxel-Level `phi_tilde_vox` (Shrunk AR Estimates)**
    *   **Scope:** One AR(p) vector per voxel, regularized by shrinking towards a more stable estimate (e.g., `phi_hat_parcel` or `phi_hat_run`).
    *   **Estimation:**
        1.  Obtain initial OLS residuals `e_orig`.
        2.  **Raw per-voxel estimates (`phi_hat_vox`):** For each voxel, estimate AR(p) coefficients from its `e_orig` (e.g., fast C++ Yule-Walker for AR(1) as sketched previously, or full Yule-Walker for AR(p)).
        3.  **Shrinkage Target (`phi_bar_ref`):** Use `phi_hat_parcel` (if available) or `phi_hat_run`.
        4.  **Empirical Bayes/James-Stein Shrinkage:** Compute `phi_tilde_vox = w_vox * phi_bar_ref_for_vox + (1-w_vox) * phi_hat_vox`. The shrinkage intensity `w_vox` is data-driven, balancing the stability of `phi_bar_ref` with the local detail in `phi_hat_vox`.
    *   **Application:** Whiten `Y_orig` voxel-wise using `phi_tilde_vox`. Whiten `X_orig` using `phi_hat_run`.
    *   **Cost:** Adds per-voxel AR estimation (O(p\*T\*V) for AR(p)) and shrinkage calculations (O(V)). Whitening remains O(T\*V).
    *   **Benefit:** Most spatially precise whitening. Crucial for analyses demanding high local fidelity where noise structure might vary at a fine scale.

**API Integration:**

This could be integrated via an `ar_scope` argument in `fmri_lm`:
`fmri_lm(..., cor_struct = "arp", ar_p = 1, ar_scope = c("run", "parcel", "voxel"), parcels = NULL, ...)`

*   `ar_scope = "run"`: Implements Tier 0 (main AR(p) proposal).
*   `ar_scope = "parcel"`: Implements Tier 1. `parcels` would be a user-supplied `NeuroVol` mask or path to one. If `NULL`, a simple internal tissue segmentation could be attempted (e.g., k-means on mean EPI).
*   `ar_scope = "voxel"`: Implements Tier 2.

**Performance Implications:**

The computational cost remains dominated by operations linear in the number of voxels (V) and timepoints (T).
*   **Tier 1:** Marginal increase over Tier 0.
*   **Tier 2:** Most expensive of the adaptive tiers, but significantly faster than full per-voxel REML/GLS approaches. The per-voxel AR(1) estimation is comparable in cost to calculating temporal variance per voxel.

**Synergy with Row-Robust GLM:**

The locally adaptive AR modeling and the row-robust GLM are complementary:
1.  **Order of Operations:** AR prewhitening (potentially locally adaptive) would generally be applied first to handle serial correlations. The row-robust IRLS procedure would then operate on these (approximately) whitened residuals to address remaining widespread temporal artifacts.
2.  **Residuals for Robust Weighting:** The initial OLS residuals used by the row-robust GLM (`Y_orig - X_orig %*% Beta_OLS`) would be computed *before* any AR whitening. The subsequent IRLS iterations and WLS fits would then operate on data that has *not* been AR-whitened by this adaptive scheme if both are enabled.
    *   *Alternatively, and perhaps more coherently,* the row-robust GLM could operate on the *AR-whitened* data. This means `Y` and `X` fed into `fast_rlm_run` would be `Y_whitened_adaptive` and `X_whitened_run`. The `sigma_hat` for robust weights would then be estimated from residuals of OLS on this *already AR-whitened* data. This seems more statistically principled, as it attempts to make residuals i.i.d. *before* robustly down-weighting outlier timepoints. This needs careful consideration regarding the interpretation of `sigma_robust`.

**Conclusion for Appendix:**

Locally adaptive AR prewhitening offers a pathway to further refine fMRI noise modeling within the `fmrireg` fast GLM framework. By building upon the efficient run-level AR estimation and C++ filtering routines, spatial variations in autocorrelation can be accommodated with manageable increases in computation, providing a valuable extension for advanced use cases. Its interaction with the proposed row-robust GLM, particularly the stage at which robust weighting is applied (pre- or post-AR-whitening), warrants careful design.

---

This appendix situates the locally adaptive AR idea as a natural, performance-conscious extension, connecting it to the main proposals and highlighting potential interactions and design choices.