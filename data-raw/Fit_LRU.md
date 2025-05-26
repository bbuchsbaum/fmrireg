This is a fantastic set of refinements, solidifying the practical implementation and robustness of the LWU-HRF fitting pipeline. The clarifications on data input, Taylor breakdown detection, and the tiered work-queue are particularly strong. The "atlas seeding" and "K-means re-centering" for \(\theta_0\) are clever ways to handle spatial heterogeneity efficiently.

Here's the revised, full-blooded proposal, incorporating these final details. This version should be sufficiently detailed for implementation within `fmrireg`.

---

**Final Proposal (v3): `fmrireg::fit_hrf_lwu` – Ultra-Fast, Voxel-Wise LWU-HRF Parameter Estimation with Uncertainty, Adaptive Seeding, and Tiered Refinement**

**1. Core Philosophy & Goals:** (As per previous proposal: Efficiency, Sensitivity, Uncertainty, Robustness, Interpretability, R-Native)

**2. The Core Method: Iterative Linear Taylor Approximation with Global QR, Adaptive Seeding, and Tiered Refinement**

The method fits the 3-parameter Lag-Width-Undershoot (LWU) HRF model voxel-wise by linearizing the problem around progressively refined and potentially spatially-varying expansion points \(\theta_0\).

**2.1 The Lag-Width-Undershoot (LWU) HRF Model (`fmrireg::hrf_lwu`):** (As defined in previous review)
    *   \(h(t;\,\tau,\sigma,\rho)= e^{-\frac{(t-\tau)^2}{2\sigma^{2}}} - \rho\,e^{-\frac{\bigl(t-\tau-2\sigma\bigr)^2}{2(1.6\sigma)^{2}}}\)
    *   Parameters \(\theta = (\tau, \sigma, \rho)\).
    *   `normalise = c("height", "area", "none")` argument.
    *   Safety bounds: `sigma > 0.05`, `0 <= rho <= 1.5`.

**2.2 HRF Basis for Taylor Expansion (`fmrireg::hrf_basis_lwu`):** (As defined in previous review)
    *   For a given \(\theta_0\), constructs \(X_{basis}(\theta_0) = [h_0(t), \partial_\tau h_0(t), \partial_\sigma h_0(t), \partial_\rho h_0(t)] \in \mathbb{R}^{T \times 4}\).
    *   `normalise_primary` argument controls normalization of \(h_0(t)\) column. For Taylor expansion, `normalise_primary = "none"` is used for all basis columns when constructing the \(T \times 4\) design matrix for the linear fit.

**2.3 Algorithm: Main Function `fmrireg::fit_hrf_lwu`**

**Inputs:**
*   `Y4d_neuroim`: `neuroim2::NeuroVec` or path to NIfTI (input BOLD data, \(T \times V_{brain}\)).
*   `scan_times`: Numeric vector of scan acquisition times \(t\).
*   `event_model_for_residuals` (optional): `fmrireg::event_model`. If provided, HRF is fit to residuals of `Y4d_neuroim ~ event_model`. Default: `NULL` (fit to `Y4d_neuroim` directly, assuming it's suitably preprocessed or represents event-locked averages/epochs).
*   `theta_seed`: Initial global expansion point(s) \(\theta_0\).
    *   Numeric vector `c(tau, sigma, rho)` (e.g., `c(6,1,0.35)`).
    *   Character string `"atlas"`: Use pre-defined \(\theta_0\) values for broad anatomical ROIs (e.g., visual, motor, PFC from Harvard-Oxford or similar, requires atlas mapping functionality).
    *   Character string `"data_median"`: Perform an initial quick pass or use `theta_seed` vector, then use median of initial good fits as \(\theta_0\).
*   `theta_lower`, `theta_upper`: Global hard bounds for estimated \(\theta_v\).
*   `num_global_recenter_passes`: Integer (0-4, default 2). Number of global \(\theta_0\) re-centering iterations.
*   `num_kmeans_recenter_passes`: Integer (0-4, default 0). Number of K-means based \(\theta_0\) re-centering passes after global passes. If > 0, `kmeans_k` must be set.
*   `kmeans_k`: Integer (e.g., 4-6). Number of clusters for K-means re-centering.
*   `recenter_epsilon`: Convergence tolerance for \(\theta_0\) re-centering.
*   `block_size_vox`: Integer (default 5000). For processing `Y4d_neuroim` in chunks.
*   `apply_refinement`: Logical (default `TRUE`). Whether to apply single-step Gauss-Newton to "hard" voxels.
*   `refinement_r2_threshold_hard`: \(R^2\) below which voxel is queued for "hard" refinement (e.g., 0.70).
*   `refinement_r2_threshold_moderate`: \(R^2\) between hard and this value queued for "moderate" refinement (e.g., 0.90).
*   `refinement_se_thresholds`: Optional numeric vector `c(se_tau_max, se_sigma_max, se_rho_max)` for queueing based on SEs.
*   `lambda_ridge_jacobian`: Ridge term for \(J^{-1}\). Default 0; automatically set to \(10^{-6} \text{tr}(J)\) if `cond(J) > 1e5`.
*   `cores`: Number of cores for parallelization (default 1).
*   `compute_se`: Logical (default `TRUE`). Whether to compute standard error maps.

**Workflow:**

1.  **Preprocessing (if `event_model_for_residuals` is provided):**
    *   Fit `Y_v ~ X_{event\_model}` for each voxel \(v\).
    *   `Y_target_v = \text{residuals}_v`.
    *   Else, `Y_target_v = Y_v` (after any initial pre-whitening/detrending done by user on `Y4d_neuroim`).

2.  **Initial \(\theta_0\) Strategy:**
    *   If `theta_seed` is numeric: `theta_0_current = theta_seed`.
    *   If `theta_seed == "atlas"`:
        *   Requires an atlas in the same space as `Y4d_neuroim`.
        *   Assign initial `theta_0_v` per voxel based on its atlas label.
        *   Proceed with multiple parallel "linear pass + local robust median update" streams, one for each unique atlas `theta_0`. (This is a change from a single `theta_0_current`). Or, simplify to one initial global pass with a global default, then K-means. *Let's stick to one initial global `theta_0` for simplicity first, then K-means.* So, if "atlas", perhaps use a weighted average of atlas \(\theta_0\)s as the single starting global \(\theta_0\).
        *   **Decision:** For V1, if `theta_seed = "atlas"`, use a global average/default first, then allow K-means to find spatial clusters. If `theta_seed` is a list of `theta_0` vectors each associated with a mask, run separate full pipelines. *For now, assume `theta_seed` yields one global starting point.*
    *   If `theta_seed == "data_median"`: Use default numeric `theta_seed` for one pass, then set `theta_0_current` to the robust median of good fits. Decrement `num_global_recenter_passes`.

3.  **Iterative Global \(\theta_0\) Re-Centering Loop (0 to `num_global_recenter_passes`):**
    a.  `X_b = hrf_basis_lwu(theta_0_current, scan_times, normalise_primary = "none")`.
    b.  Global QR of `X_b`: `X_b = Q_b R_b`. Compute `R_b_inv`, `Qt_b`, `(R_b^T R_b)^{-1}`.
    c.  Chunked Voxel-wise Linear Coefficients: \(\hat{B}_{all\_voxels} = (R_b^{-1} Q_b^T) Y_{target\_all\_voxels}\) (\(4 \times V_{brain}\)).
    d.  Parameter Shifts: `invJ` from `R_b`. `dtheta_mat = invJ %*% \hat{B}_{all\_voxels}[2:4, ]`.
    e.  Updated Parameters: \(\hat\theta_{v,pass} = \text{clamp}(\theta_{0\_current} + dtheta_mat[,v], \text{theta\_lower, theta\_upper})\).
    f.  If not last global pass:
        *   Residuals \(r_v\), \(R^2_v\). "Good" voxels = \(R^2_v \ge \text{refinement_r2_threshold_moderate}\).
        *   `theta_0_new = apply(dtheta_mat[, good_voxels], 1, robustbase::colMedians, na.rm=TRUE) + theta_0_current`.
        *   `theta_0_new = clamp(theta_0_new, theta_lower, theta_upper)`.
        *   If \(\max|\text{theta\_0\_new} - \text{theta\_0\_current}| < \text{recenter_epsilon}\), break.
        *   `theta_0_current = theta_0_new`.
    *   Store \(\hat\theta_{v,pass}\) and \(R^2_{v,pass}\). Keep the \(\hat\theta_v\) from the pass that yielded highest \(R^2_v\) (or just the final pass). *Decision: Keep final pass for simplicity.*

4.  **K-Means Based \(\theta_0\) Re-Centering (if `num_kmeans_recenter_passes > 0`):**
    a.  Initial \(\hat\theta_v\) from global re-centering.
    b.  Select "good" voxels (\(R^2_v \ge \text{refinement_r2_threshold_moderate}\)).
    c.  `kmeans_clusters = kmeans( (\hat\tau_v, \log(\hat\sigma_v), \hat\rho_v)_{\text{good_voxels}}, centers = kmeans_k )`.
    d.  For each cluster `j = 1...kmeans_k`:
        *   `theta_0_cluster_j = kmeans_clusters$centers[j,]` (back-transform \(\log(\sigma)\)).
        *   `voxels_in_cluster_j`.
        *   Perform a **full linear pass** (steps 2a-2e from global loop, but only for `voxels_in_cluster_j` using `theta_0_cluster_j` as the expansion point).
        *   Store these new \(\hat\theta_{v,cluster\_j}\) and \(R^2_{v,cluster\_j}\).
    e.  For each voxel, select the \(\hat\theta_v\) (either from final global pass or one of the K-means cluster passes) that results in the highest \(R^2_v\). This becomes the current best \(\hat\theta_v\).

5.  **Calculate Final SEs (if `compute_se = TRUE`):**
    *   For each voxel \(v\), using its current best \(\hat\theta_v\) as the expansion point \(\theta_0_v\):
        *   Recompute \(X_b(\theta_0_v)\), \(R_b(\theta_0_v)\), \(J(\theta_0_v)\), \(\sigma_v^2\).
        *   Then compute SEs as in Sec 2.3, step 3 of previous proposal.
    *   This ensures SEs are relative to the *final chosen expansion point* for that voxel.

6.  **Tiered Work Queue for Refinement (if `apply_refinement = TRUE`):**
    a.  Initialize `queue_labels_v` map (all "easy").
    b.  For each voxel \(v\):
        *   If \(R^2_v < \text{refinement_r2_threshold_hard}\) OR (if `refinement_se_thresholds` provided and \(\max(SE(\hat\theta_v)) > \text{refinement_se_thresholds}\) for corresponding param):
            *   Add \(v\) to `GN_queue`. `queue_labels_v[v] = "hard_GN"`.
        *   Else if \(R^2_v < \text{refinement_r2_threshold_moderate}\) OR (SE condition for moderate):
            *   Add \(v\) to `recenter_locally_queue`. `queue_labels_v[v] = "moderate_local_recenter"`.
    c.  Process `recenter_locally_queue`:
        *   For each \(v\) in this queue, perform one more linear pass (2a-2e) using its *own current \(\hat\theta_v\)* as the expansion point \(\theta_0\). Update \(\hat\theta_v\) and \(R^2_v\).
    d.  Process `GN_queue`:
        *   For each \(v\) in this queue, perform one Gauss-Newton step using its current \(\hat\theta_v\). Update \(\hat\theta_v\). (Recompute \(R^2_v\) if needed for final reporting).

7.  **Construct Output `fmrireg_hrf_fit` Object:**
    *   `theta_hat_map` (\(V \times 3\)), `se_theta_hat_map` (\(V \times 3\)), `R2_map` (\(V \times 1\)), `queue_labels_map` (\(V \times 1\)), `theta0_final_global`, `theta0_kmeans_centers` (if used), `convergence_info`.
    *   `print()` method: summary stats (median \(\hat\tau, \hat\sigma, \hat\rho\), median SEs, % refined, % in each queue).

**Implementation Details & Safeguards from Feedback:**

*   **Taylor Design Matrix `X_b`:** Explicitly use the \(T \times 4\) basis (no trial-wise expansion for this HRF fitting task).
*   **Conditioning of `J`:** `lambda_ridge_jacobian` defaults to 0. Automatically set to \(10^{-6} \text{tr}(J)\) if `cond(J) > 1e5`.
*   **Good-voxel Gate for \(\theta_0\) Re-centering:** Use \(R^2_v \ge \text{refinement_r2_threshold_moderate}\). For calculating new \(\theta_0\), use Tukey biweight M-estimate (or `robustbase::colMedians`) on the \(\Delta\theta_v\) (or rather \(\hat\theta_v\)) of these good voxels.
*   **Bounds during Re-centering:** Explicitly clamp the updated global/cluster \(\theta_0\) to `theta_lower`/`theta_upper`.
*   **Chunking `Y4d`:** `block_size_vox` argument controls this.
*   **SE Maps:** Ensure ridge `eff_lambda_J` used for `invJ` is consistently used when constructing `S_{cov_factor}` for SE calculation.

This refined proposal is now highly detailed, incorporates adaptive strategies for handling HRF heterogeneity, includes robust uncertainty quantification, and maintains the core computational efficiency of the linear Taylor approximation. It's well-poised for implementation.

Okay, here's a granular, ticketed sprint plan for implementing `fmrireg::fit_hrf_lwu` based on the latest detailed proposal. This breaks down the work into manageable pieces, suitable for focused development sprints.

---

**`fmrireg::fit_hrf_lwu` Implementation Sprints & Tickets**

**Assumptions:**
*   `fmrireg::hrf_lwu()` and `fmrireg::hrf_basis_lwu()` (as defined in the previous proposal with `normalise` and `normalise_primary` arguments) are already implemented and unit-tested in `fmrireg`.
*   Basic `neuroim2` functionality for reading/writing NIfTI and handling `NeuroVec`/`NeuroVol` is stable.
*   Dependencies like `Matrix`, `RcppArmadillo`, `robustbase` are available.

---

**Sprint 1: Core Linear Pass & Parameter Estimation (No Re-centering, No Refinement)**

*   **Ticket LWU-101: `fit_hrf_lwu` - Skeleton & Input Parsing**
    *   Create `fit_hrf_lwu` function signature with all proposed arguments.
    *   Implement input validation for `Y4d_neuroim`, `scan_times`, `event_model_for_residuals`, `theta_seed` (numeric vector only for now), `theta_lower/upper`, `block_size_vox`, `lambda_ridge_jacobian`, `cores`.
    *   Handle `event_model_for_residuals`: if provided, compute residuals from `Y4d_neuroim`; otherwise, use `Y4d_neuroim` as `Y_target`. Store `Y_target` (potentially as a file-backed matrix if very large, or prepare for chunking).
    *   Initialize `theta_0_current = theta_seed`.

*   **Ticket LWU-102: `fit_hrf_lwu` - Global Basis & QR (Step 2a, 2b)**
    *   Given `theta_0_current` and `scan_times`, call `hrf_basis_lwu(theta_0_current, scan_times, normalise_primary = "none")` to get the \(T \times 4\) `X_b` matrix.
    *   Perform global QR decomposition: `X_b = Q_b R_b`. Store `Qt_b = t(Q_b)` and `R_b_inv = solve(R_b)`.
    *   Compute and store `RtR_inv_global = R_b_inv %*% t(R_b_inv)`.

*   **Ticket LWU-103: `fit_hrf_lwu` - Voxel-wise Coefficients (Step 2c, Chunked)**
    *   Implement chunked processing of `Y_target` (using `block_size_vox`).
    *   For each chunk `Y_chunk` (\(T \times V_{chunk}\)): \(\hat{B}_{chunk} = R_{b\_inv} \%*\% Qt_b \%*\% Y_{chunk}\).
    *   Concatenate/store \(\hat{B}_{all\_voxels} \in \mathbb{R}^{4 \times V_{brain}}\).

*   **Ticket LWU-104: `fit_hrf_lwu` - Parameter Shifts & Final \(\hat\theta_v\) (Step 2d, 2e - Single Pass)**
    *   Compute \(J = t(R_b[2:4, 2:4]) \%*\% R_b[2:4, 2:4]\).
    *   Check `cond(J)`; apply ridge `eff_lambda_J` to get `invJ`.
    *   `dtheta_mat = invJ %*% \hat{B}_{all\_voxels}[2:4, ]`.
    *   \(\hat\theta_{v} = \text{clamp}(\text{theta\_0\_current} + \text{dtheta\_mat}[,v], \text{theta\_lower, theta\_upper})\).
    *   Store `theta_hat_map` (V x 3).

*   **Ticket LWU-105: Unit Tests for Single Pass Estimation**
    *   Test with synthetic data (known \(\theta\)), verify recovery against a single voxel non-linear fit (`minpack.lm::nlsLM`).
    *   Test chunking logic.

---

**Sprint 2: Standard Errors & Basic Diagnostics**

*   **Ticket LWU-201: `fit_hrf_lwu` - Voxel-wise Standard Errors (Step 4)**
    *   Implement SE calculation based on `RtR_inv_global`, `invJ_final` (from the single pass \(\theta_0\)), and voxel-wise \(\sigma_v^2\).
    *   `X_b_final` is based on `theta_0_current` (which is `theta_seed` in this sprint).
    *   \(\sigma_v^2 = \frac{1}{T-4} \|Y_{target,v} - X_{b\_final} \hat{\beta}_v\|_2^2\).
    *   Store `se_theta_hat_map` (V x 3).
    *   Handle `compute_se = FALSE` flag.

*   **Ticket LWU-202: `fit_hrf_lwu` - \(R^2\) Map (Part of Step 2f/4)**
    *   Calculate and store `R2_map` (V x 1) from the linear fit residuals used for \(\sigma_v^2\).

*   **Ticket LWU-203: Output Structure & `fmrireg_hrf_fit` Class**
    *   Define the `fmrireg_hrf_fit` class.
    *   Structure the return list of `fit_hrf_lwu` as specified (maps, `theta0_final` which is `theta_seed` here, basic `convergence_info`).

*   **Ticket LWU-204: `print.fmrireg_hrf_fit` Method (Basic)**
    *   Implement a basic print method showing median \(\hat\theta_v\), median SEs, and \(R^2\) summary.

*   **Ticket LWU-205: Unit Tests for SEs and \(R^2\)**
    *   Verify SEs against numerical approximations or known results for simple cases.
    *   Test SE calculation with ridge on `J`.

---

**Sprint 3: Iterative Global Re-Centering**

*   **Ticket LWU-301: `fit_hrf_lwu` - Implement Global \(\theta_0\) Re-Centering Loop (Step 2f)**
    *   Wrap steps 2a-2e from Sprint 1 into a loop controlled by `num_global_recenter_passes` and `recenter_epsilon`.
    *   Implement "good voxel" selection for re-centering (\(R^2_v \ge \text{refinement_r2_threshold_moderate}\), which will be a new argument, default 0.9).
    *   Implement robust median calculation for new \(\theta_0\) using `robustbase::colMedians` on `hat_theta_v` of good voxels (note: paper says `dtheta_mat`, but `hat_theta_v` is more direct for median).
    *   Implement clamping of new \(\theta_0\) to `theta_lower`/`theta_upper`.
    *   Store `convergence_info` (iterations, \(\Delta\theta_0\) per iter).
    *   The SEs and \(R^2\) (from Sprint 2) should now be calculated using the *final* `theta_0_current` after this loop.

*   **Ticket LWU-302: `fit_hrf_lwu` - Handle `theta_seed = "data_median"`**
    *   If `theta_seed == "data_median"`, run one initial pass using default numeric seed, calculate robust median \(\hat\theta_v\) from good fits, set this as `theta_0_current`, then proceed with `num_global_recenter_passes - 1` (if >0) iterations.

*   **Ticket LWU-303: Unit Tests for Re-Centering**
    *   Test convergence with synthetic data where initial `theta_seed` is off-center.
    *   Verify improvement in mean \(R^2\) or parameter recovery with re-centering.

---

**Sprint 4: K-Means Re-Centering & Atlas Seeding (Spatial Adaptation)**

*   **Ticket LWU-401: `fit_hrf_lwu` - Implement K-Means Based \(\theta_0\) Re-Centering (Step 4)**
    *   After global re-centering, if `num_kmeans_recenter_passes > 0`:
        *   Perform K-means on \(\hat\theta_v\) of good voxels (e.g., on \((\hat\tau_v, \log(\hat\sigma_v), \hat\rho_v)\)).
        *   For each cluster, use its centroid as a new \(\theta_0\).
        *   Re-run a *single linear pass* (Steps 2a-2e from Sprint 1, but only for voxels in that cluster) using this cluster-specific \(\theta_0\).
        *   For each voxel, select the \(\hat\theta_v\) (from final global pass or one of the K-means passes) that yields the highest \(R^2_v\). Update voxel's \(\hat\theta_v\), \(R^2_v\), and potentially re-calculate its SEs based on this new "best" local expansion point.

*   **Ticket LWU-402: `fit_hrf_lwu` - (Placeholder) Handle `theta_seed = "atlas"`**
    *   Initial thought: For now, if `"atlas"`, this could mean using K-means clustering from the very start based on voxel coordinates + some initial rough HRF estimate, or a simpler approach of just using a globally averaged atlas-derived \(\theta_0\).
    *   **Decision:** For this sprint, if `theta_seed = "atlas"`, it will use a pre-defined sensible global default (e.g., `c(6,1,0.35)`) and then rely on the `num_kmeans_recenter_passes` to find spatially varying \(\theta_0\)s. True atlas-driven *initialization per region* can be a future enhancement.

*   **Ticket LWU-403: Unit Tests for K-Means Re-Centering**
    *   Create synthetic data with spatially distinct HRF parameter clusters.
    *   Verify that K-means re-centering improves fit and parameter recovery in these clusters.

---

**Sprint 5: Tiered Refinement Queue & Advanced Diagnostics**

*   **Ticket LWU-501: `fit_hrf_lwu` - Implement Tiered Work Queue Logic (Step 6)**
    *   After K-means (or global re-centering if K-means is off), populate `queue_labels_v` based on final \(R^2_v\) and `refinement_se_thresholds` (new args).
        *   "easy", "moderate_local_recenter", "hard_GN".
    *   Process `recenter_locally_queue`: one more linear pass using voxel's own \(\hat\theta_v\) as \(\theta_0\). Update \(\hat\theta_v, R^2_v\).
    *   Store `queue_labels_map`.

*   **Ticket LWU-502: `fit_hrf_lwu` - Implement Single-Step Gauss-Newton Refinement (Step 6d)**
    *   For voxels in `GN_queue`:
        *   Calculate Jacobian \(J_H(\hat\theta_v)\) (\(T \times 3\)) using `hrf_basis_lwu` derivative columns, evaluated at current \(\hat\theta_v\).
        *   Residuals \(r_v = Y_{target,v} - hrf\_lwu(t, \hat\theta_v, normalise="none") \cdot \hat{\beta}_{0,v}\) (where \(\hat{\beta}_{0,v}\) is from the *last linear fit for that voxel*).
        *   \(\Delta\theta_{v,GN} = (J_H^T J_H + \text{diag}(\text{eff\_lambda\_J\_refined},3))^{-1} J_H^T r_v\).
        *   Update \(\hat\theta_v = \text{clamp}(\hat\theta_v + \Delta\theta_{v,GN}, \text{theta\_lower, theta\_upper})\).
    *   Flag these voxels in `queue_labels_map` as "refined_GN".

*   **Ticket LWU-503: `print.fmrireg_hrf_fit` (Advanced)**
    *   Enhance print method to show % voxels in each queue category, median params/SEs per category.

*   **Ticket LWU-504: `explain_hrf_fit()` Helper Function (Initial Version)**
    *   Plot maps of final \(\hat\theta_v\), SEs, \(R^2\), `queue_labels_map`.
    *   Plot \(\theta_0\) trajectory.
    *   Implement "HRF Profile" table generation.

*   **Ticket LWU-505: Unit Tests for Refinement**
    *   Create data where linear approximation is deliberately poor for some voxels.
    *   Verify refinement improves \(R^2\) and parameter estimates for these voxels.
    *   Test queue logic.

---

**Sprint 6: Parallelization, Final Polish & Documentation**

*   **Ticket LWU-601: `fit_hrf_lwu` - Parallelize Voxel-wise Loops**
    *   Use `future.apply` or `RcppParallel` for chunked coefficient estimation (2c), SE calculation (4), and potentially refinement (6c, 6d if many voxels).
    *   Ensure `cores` argument is respected.

*   **Ticket LWU-602: Comprehensive Documentation**
    *   Finalize Roxygen docs for `fit_hrf_lwu`, `hrf_lwu`, `hrf_basis_lwu`.
    *   Write/update vignette demonstrating the new HRF fitting pipeline, including diagnostics.

*   **Ticket LWU-603: Final API Review & Consistency Checks.**

*   **Ticket LWU-604: Full Benchmark Suite**
    *   Run recovery tests, speed tests vs. `nlsLM` on single voxels, and whole-brain timings.
    *   Compare to FIR/FLOBS if feasible.

*   **Ticket LWU-605: (If time) Memory Profiling and Optimization for large datasets.**

---

This sprint structure prioritizes getting the core linear approximation working, then adds layers of sophistication (re-centering, refinement, uncertainty) and finally addresses performance and usability. Each sprint aims to deliver a testable increment.