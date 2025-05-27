Okay, this is excellent. The responses clarify all previous ambiguities and the proposed patches/adjustments are sensible and maintain the efficiency goals.

Here's a draft proposal followed by a ticketed sprint plan:

## Proposal: High-Performance AR(p) Prewhitening for `fmri_lm`

**Preamble:**

The `fmrireg` package currently offers a highly optimized "fast path" for fitting fMRI linear models assuming i.i.d. errors. This proposal outlines an extension to incorporate AR(p) serial correlation modeling via prewhitening, maintaining the core speed advantages of the existing pipeline. The method involves estimating run-level AR coefficients from OLS residuals, applying a causal filter to the data (Y) and design matrix (X) for that run, and then re-fitting the GLM on the whitened data. This approach offers a substantial portion of the statistical efficiency of full REML/GLS at a fraction of the computational cost.

**1. High-Level Design:**

The process for a single run (or globally, if specified) will be:

1.  **(If `iter == 1`):** Perform an initial OLS fit on the *original* (unwhitened) `Y_run` and `X_run`.
    *   Obtain OLS residuals: `e_OLS = Y_run_orig - X_run_orig %*% Beta_OLS`.
    *   Estimate AR(p) coefficients (`phi_hat`) from the run-mean `e_OLS` using Yule-Walker.
2.  **Whitening:** Apply an AR(p) causal filter to `Y_run_orig` and `X_run_orig` using `phi_hat` to produce `Y_run_whitened` and `X_run_whitened`. This is done in-place via an efficient C++ routine.
3.  **GLS Fit:** Perform a GLM fit on `Y_run_whitened` and `X_run_whitened` using the existing `.fast_lm_matrix` machinery. The resulting `Betas` and `sigma^2` are the GLS estimates.
4.  **(Optional, if `cor_iter > 1`):** Re-estimate `phi_hat` from the residuals of the current GLS fit (`Y_run_whitened - X_run_whitened %*% Beta_GLS`) and repeat steps 2-3.

**2. Key Mechanisms & API Changes:**

*   **API Extension:**
    *   `fmri_lm(..., cor_struct = c("iid", "ar1", "ar2", "arp"), cor_iter = 1, cor_global = FALSE, ar_p = NULL, ar1_exact_first = FALSE, ...)`
    *   `cor_struct`: Specifies the correlation structure. "arp" allows user-defined AR order `p` via `ar_p` argument. Initially, only AR components will be handled (no MA).
    *   `cor_iter`: Number of iterations for GLS estimation (1 = initial OLS + 1 GLS fit).
    *   `cor_global`: If `TRUE`, a single `phi_hat` is estimated from all runs' initial residuals and applied globally. If `FALSE` (default), `phi_hat` is run-specific.
    *   `ar_p`: Integer, specifies the order for `cor_struct = "arp"`.
    *   `ar1_exact_first`: If `TRUE` and `cor_struct = "ar1"`, scales the first observation by `sqrt(1-phi^2)`. Default `FALSE`.

*   **AR Estimation (`.estimate_ar`):**
    *   R helper function using `stats::acf` and `solve(toeplitz(...))` for Yule-Walker.
    *   Input: vector of residuals (typically run-mean OLS residuals).
    *   Output: vector of AR coefficients `phi_hat`.

*   **Whitening Filter (C++):**
    *   New C++ functions (e.g., `ar_whiten_inplace(arma::mat& Y, arma::mat& X, const arma::vec& phi, bool exact_first_ar1)`) using OpenMP.
    *   Performs causal filtering: `y_new(t) = y_old(t) - sum_{k=1 to p} (phi_k * y_old(t-k))`.
    *   Handles `Y` (data) and `X` (design matrix columns) in-place.
    *   Initial `p` values (for `t < p`) are based on `prev_k = 0.0`.
    *   Optionally scales first sample for AR(1) if `exact_first_ar1` is true.

*   **GLM Engine (`.fast_lm_matrix` Modification):**
    *   Signature amended: `.fast_lm_matrix(X, Y, proj, return_fitted = FALSE)`
    *   If `return_fitted = TRUE`, it computes and returns `Fitted = X %*% Betas`. RSS is then `colSums((Y - Fitted)^2)`.
    *   Otherwise (default), `Fitted` is `NULL`, and `rss` is computed via the faster `yTy - beta_XtX_beta` formula.

*   **Integration into `runwise_lm` and `chunkwise_lm`:**
    *   **Initial OLS & `phi_hat` Estimation (if `iter == 1`):**
        *   Call `.fast_lm_matrix(X_run_orig, Y_run_orig, proj_run_orig, return_fitted = TRUE)` to get `ols$betas` and `ols$fitted`.
        *   Compute `Resid_OLS = Y_run_orig - ols$fitted`.
        *   Estimate `phi_hat` using `.estimate_ar(rowMeans(Resid_OLS), p = ar_order)`.
        *   For `cor_global = TRUE`, `phi_hat` estimation happens after an initial OLS pass over *all* runs, aggregating residuals.
    *   **Whitening:** Call `ar_whiten_inplace(Y_run_orig, X_run_orig, phi_hat, ...)` to modify data and design *for the current run*.
        *   For `chunkwise_lm`, whitening of `Y_run_orig` and `X_run_orig` happens *once per run* before chunking for the GLM step.
    *   **GLS Fit:** Call `.fast_preproject(X_run_whitened)` to get `proj_run_whitened`. Then call `.fast_lm_matrix(X_run_whitened, Y_run_whitened, proj_run_whitened)` to get GLS betas and sigma.
    *   **Iteration (if `cor_iter > 1`):**
        *   Compute residuals from the current GLS fit: `Resid_GLS = Y_run_whitened - X_run_whitened %*% Betas_GLS`.
        *   Re-estimate `phi_hat` using `.estimate_ar(rowMeans(Resid_GLS), p = ar_order)`.
        *   Loop back to **Whitening** using the *original* `Y_run_orig` and `X_run_orig` with the *new* `phi_hat`.

*   **Statistical Inference:**
    *   The GLS estimates `Beta_GLS` and `sigma^2_GLS` from the whitened model are used directly.
    *   `Cov(Beta_GLS) = sigma^2_GLS * (X_whitened^T X_whitened)^-1`. This `(X_whitened^T X_whitened)^-1` is already available from `proj_run_whitened$XtXinv`.

**3. Performance:**

*   The AR estimation step is `O(n_timepoints)` per run (negligible).
*   The C++ whitening filter is `O(n_timepoints * (V_voxels + P_predictors))` per run, efficiently parallelized.
*   The GLM fit on whitened data reuses the existing fast path.
*   Expected total cost for one GLS iteration: < 2x the cost of the i.i.d. fast path.

**4. Benefits:**

*   Handles AR(p) serial correlations, improving statistical validity.
*   Maintains high computational performance, suitable for large datasets.
*   Minimal changes to the existing high-level R API (`fmri_lm`).
*   Leverages established statistical approximations common in fMRI software.

---

## Ticketed Sprint to Completion:

**Epic: Implement AR(p) Prewhitening for `fmri_lm`**

**Sprint Goal:** Deliver a functional and tested AR(p) prewhitening capability within `fmri_lm`, significantly improving model accuracy for serially correlated fMRI data while maintaining high performance.

---

**Ticket 1: API Definition & Argument Threading (R)**
*   **Task:** Modify `fmri_lm()` signature to include `cor_struct`, `cor_iter`, `cor_global`, `ar_p`, `ar1_exact_first`.
*   **Task:** Thread these new arguments down through `fmri_lm_fit()` to `runwise_lm()` and `chunkwise_lm()`.
*   **Task:** Add basic input validation for new arguments (e.g., `ar_p` must be positive integer if `cor_struct == "arp"`).
*   **Definition of Done:** New arguments are present in function signatures and passed correctly. Unit tests for argument passing. Documentation stubs for new arguments.
*   **Estimate:** 0.5 days

**Ticket 2: AR Coefficient Estimation (R)**
*   **Task:** Implement R helper function `.estimate_ar(residuals_vec, p_order)` using `stats::acf` and Yule-Walker equations (`solve(toeplitz(...))`).
*   **Task:** Unit test `.estimate_ar` with known AR(1) and AR(2) series to verify coefficient recovery.
*   **Definition of Done:** `.estimate_ar` function implemented and unit tested.
*   **Estimate:** 0.5 days

**Ticket 3: C++ Whitening Filter Implementation (C++/RcppArmadillo)**
*   **Task:** Implement `ar_whiten_inplace(arma::mat& Y, arma::mat& X, const arma::vec& phi_coeffs, bool exact_first_ar1)` in C++.
    *   Handle general AR(p) based on `phi_coeffs.n_elem`.
    *   Implement `prev_k` rolling buffer for past values.
    *   Include OpenMP pragmas for parallelization over voxels (Y) and predictors (X).
    *   Implement optional `exact_first_ar1` scaling for AR(1) case.
*   **Task:** Expose this function to R via Rcpp.
*   **Task:** Unit test the C++ filter directly with simple matrices and known `phi` values, checking output against manual calculation for AR(1) and AR(2). Test `exact_first_ar1`.
*   **Definition of Done:** C++ filter implemented, Rcpp-exposed, and unit tested for correctness on small examples.
*   **Estimate:** 1.5 days

**Ticket 4: Modify `.fast_lm_matrix` (R)**
*   **Task:** Amend signature: `.fast_lm_matrix(X, Y, proj, return_fitted = FALSE)`.
*   **Task:** Implement conditional calculation: if `return_fitted`, compute `Fitted = X %*% Betas` and `rss = colSums((Y - Fitted)^2)`.
*   **Task:** Ensure `list(..., fitted = if(return_fitted) Fitted else NULL, ...)` is returned.
*   **Definition of Done:** `.fast_lm_matrix` updated and existing unit tests (if any specifically for it) still pass or are adapted.
*   **Estimate:** 0.5 days

**Ticket 5: Integrate Whitening into `runwise_lm` (R)**
*   **Task:** Implement the iterative loop logic within `runwise_lm`:
    *   If `cor_struct != "iid"`:
        *   Initial OLS pass (`iter == 1` or if `phi_hat` is `NULL`):
            *   Call `.fast_lm_matrix(X_run_orig, Y_run_orig, proj_run_orig, return_fitted = TRUE)`.
            *   Compute `Resid_OLS = Y_run_orig - ols$fitted`.
            *   Estimate `phi_hat_run = .estimate_ar(rowMeans(Resid_OLS), p = ar_order)`.
            *   (If `cor_global = TRUE`, this `phi_hat_run` contributes to a global estimate; actual global estimation logic is part of Ticket 7).
        *   Make copies `X_run_iter = X_run_orig`, `Y_run_iter = Y_run_orig`.
        *   Call `ar_whiten_inplace(Y_run_iter, X_run_iter, phi_hat_run, ar1_exact_first)`.
        *   `proj_run_iter = .fast_preproject(X_run_iter)`.
        *   `gls_results = .fast_lm_matrix(X_run_iter, Y_run_iter, proj_run_iter)`.
        *   Update `Betas`, `sigma`, etc. from `gls_results` for this iteration's output.
        *   If `cor_iter > iter`:
            *   Compute `Resid_GLS = Y_run_iter - X_run_iter %*% gls_results$betas`.
            *   Re-estimate `phi_hat_run = .estimate_ar(rowMeans(Resid_GLS), p = ar_order)`.
*   **Task:** Ensure correct `X_run_orig`, `Y_run_orig` (unwhitened versions) are used for the *first* OLS residual calculation and as input to `ar_whiten_inplace` in *each* iteration.
*   **Task:** Ensure `proj_run_orig` (from unwhitened `X_run_orig`) is used for the initial OLS, and `proj_run_iter` (from whitened `X_run_iter`) is used for the GLS step.
*   **Definition of Done:** `runwise_lm` correctly performs OLS (if `iter==1`), estimates `phi_hat_run`, whitens, and performs GLS. Results are passed to existing pooling logic.
*   **Estimate:** 2 days

**Ticket 6: Integrate Whitening into `chunkwise_lm` (R)**
*   **Task:** Adapt `chunkwise_lm` to incorporate run-level AR estimation and whitening:
    *   Before chunking loop for a given run:
        *   Perform initial OLS on the *full run's* data (`Y_run_orig`, `X_run_orig`) to get `Resid_OLS_run`.
        *   Estimate `phi_hat_run` using `.estimate_ar(rowMeans(Resid_OLS_run), p = ar_order)`.
        *   (If `cor_global = TRUE`, this contributes to global estimate, see Ticket 7).
    *   Create whitened versions of the *full run's* data: `Y_run_w = Y_run_orig`, `X_run_w = X_run_orig`. Call `ar_whiten_inplace(Y_run_w, X_run_w, phi_hat_run, ar1_exact_first)`.
    *   The existing chunking loop in `chunkwise_lm` then operates on *slices* of `Y_run_w`.
    *   The design matrix `X_run_w` is pre-projected *once* for the whitened run: `proj_run_w = .fast_preproject(X_run_w)`.
    *   Inside the chunk loop, `.fast_lm_matrix` is called with `Y_chunk_w` (from `Y_run_w`) and `proj_run_w`.
*   **Task:** Iteration (`cor_iter > 1`) for `chunkwise_lm` would be more complex as `phi_hat` would need to be re-estimated from full-run GLS residuals. Initially, `chunkwise_lm` might support `cor_iter = 1` only, or this ticket needs careful planning for iteration. (Clarification: for iter > 1, re-estimate `phi_hat_run` from residuals reconstructed from *all chunks* of that run, then re-whiten the full run data and re-chunk).
*   **Definition of Done:** `chunkwise_lm` correctly uses run-level `phi_hat` to whiten full run data before chunking and GLM fitting. Handles at least `cor_iter = 1`.
*   **Estimate:** 2.5 days (extra 0.5 for managing iteration if included)

**Ticket 7: Implement `cor_global = TRUE` Logic (R)**
*   **Task:** Modify `fmri_lm_fit` (or a new orchestrator function called by it).
    *   If `cor_global = TRUE` and `cor_struct != "iid"`:
        *   First pass: Loop through all runs (using `runwise_lm` or `chunkwise_lm` logic for OLS part *only* up to `Resid_OLS_run`).
        *   Collect/Aggregate `Resid_OLS_run` from all runs (e.g., concatenate or average autocovariances).
        *   Estimate a single `phi_hat_global` using `.estimate_ar`.
        *   Second pass: Loop through all runs/chunks again. This time, *all* runs/chunks use `phi_hat_global` for whitening. The GLS estimation proceeds as in Tickets 5 & 6, but `phi_hat` is fixed to `phi_hat_global` and not re-estimated per run/iteration. `cor_iter` for global would typically be 1.
*   **Definition of Done:** `cor_global = TRUE` correctly estimates and applies a global `phi_hat`.
*   **Estimate:** 1 day

**Ticket 8: Unit & Integration Tests (R/testthat)**
*   **Task:** Create test cases:
    *   Simulated data with no serial correlation (`cor_struct = "iid"` should match `cor_struct = "ar1"` with `phi_hat` near 0).
    *   Simulated data with known AR(1) (e.g., `phi = 0.4`). Verify `phi_hat` recovery and that GLS standard errors are more accurate than OLS.
    *   Simulated data with known AR(2).
    *   Test `cor_global = TRUE` vs. `cor_global = FALSE`.
    *   Test `ar1_exact_first = TRUE` vs. `FALSE`.
    *   Test with `cor_iter > 1` (if implemented for `chunkwise_lm`).
*   **Task:** Compare results (betas, SEs, t-stats) against a known package (e.g., `nlme::gls` or a manual GLS calculation on a small dataset) if feasible for validation.
*   **Definition of Done:** Comprehensive test suite covering different `cor_struct` options, iterations, and global/run-specific estimation.
*   **Estimate:** 2 days

**Ticket 9: Documentation & Vignette Update (R/Rmd)**
*   **Task:** Document new `fmri_lm` arguments (`cor_struct`, `cor_iter`, etc.).
*   **Task:** Explain the AR(p) prewhitening method, its benefits, and usage.
*   **Task:** Add an example to a vignette demonstrating AR(1) correction.
*   **Task:** Clarify "arp" means AR-only in docs.
*   **Definition of Done:** User-facing documentation and vignette updated.
*   **Estimate:** 1 day

---
**Total Estimated Time:** 9 days (adjusting for `chunkwise_lm` iteration complexity)

************

Okay, this is a very interesting and practical extension to the AR(p) prewhitening proposal. It addresses a common concern with global or run-level AR models: their potential inadequacy in brain regions with heterogeneous noise characteristics. The tiered approach is sensible.

Here's the draft appendix, incorporating my spin:

---

## Appendix: Locally Adaptive Serial Correlation Modeling

**A.1 Motivation: Beyond Global AR Models**

While run-level AR(p) prewhitening (as outlined in the main proposal) significantly improves upon i.i.d. assumptions with minimal computational overhead, fMRI noise characteristics can vary spatially. Different tissue types (grey matter, white matter, CSF) and even distinct grey matter regions can exhibit different temporal autocorrelation structures. A single AR(p) coefficient vector (`phi_hat_run`) applied to an entire run might be suboptimal for some voxels.

This appendix proposes a tiered strategy to introduce locally adaptive AR(p) coefficients, allowing the whitening process to be more sensitive to these spatial variations, while still aiming to preserve the "lightning-fast" nature of the GLM pipeline.

**A.2 Tiered Strategy for AR Coefficient Estimation & Application**

The core idea is to refine the `phi_hat` used for whitening based on local information, ranging from a single run-level estimate to voxel-specific, regularized estimates.

*   **Tier 0: Run-Level `phi_hat_run` (Default)**
    *   **Scope:** One set of AR(p) coefficients per run.
    *   **Estimation:** As per the main proposal: `phi_hat_run` is estimated from the OLS residuals of the entire run (e.g., from run-mean residuals or pooled autocovariances).
    *   **Application:** The same `phi_hat_run` is used to whiten all voxels within that run.
    *   **Cost:** Already accounted for in the main proposal's O(T\*V) whitening pass.
    *   **Rationale:** This is the new baseline, offering substantial improvement over i.i.d. with minimal cost. Likely sufficient for most standard analyses.

*   **Tier 1: Parcel-Level `phi_hat_parcel`**
    *   **Scope:** One set of AR(p) coefficients per pre-defined parcel or tissue class within each run.
    *   **Estimation:**
        1.  Obtain OLS residuals for all voxels in the run.
        2.  For each parcel/tissue class:
            *   Average the OLS residuals across all voxels belonging to that parcel within the current run.
            *   Estimate `phi_hat_parcel` from these parcel-mean residuals using Yule-Walker.
    *   **Application:** All voxels within a given parcel are whitened using that parcel's specific `phi_hat_parcel`. The C++ whitening loop would use `phi_coeffs[parcel_id_of_voxel_v]` instead of a single `phi_hat_run`.
    *   **Cost:**
        *   OLS residuals: Already computed for Tier 0.
        *   Averaging residuals per parcel: O(T\*V) once (can be combined with residual calculation).
        *   Yule-Walker per parcel: O(T \* N_parcels), typically negligible as N_parcels << V.
        *   Whitening: Still O(T\*V), just with a lookup for `phi`.
    *   **Rationale:** Addresses gross differences in noise profiles (e.g., GM vs. WM vs. CSF). Requires a parcel/tissue map.

*   **Tier 2: Voxel-Level `phi_tilde_vox` (Shrunk Estimate)**
    *   **Scope:** One set of AR(p) coefficients for *each voxel*, regularized by shrinking towards a more stable estimate (e.g., parcel-mean or run-mean).
    *   **Estimation:**
        1.  Obtain OLS residuals for all voxels (`e_vox`).
        2.  **Noisy per-voxel estimates (`phi_hat_vox`):** For each voxel, compute a raw AR(p) estimate from `e_vox` using Yule-Walker (efficient C++ loop shown in user's proposal for AR(1)).
        3.  **Stable reference estimates (`phi_bar_ref`):** These could be `phi_hat_run` (from Tier 0) or `phi_hat_parcel` (from Tier 1 if parcels are provided).
        4.  **Shrinkage:** Apply James-Stein-like or empirical Bayes shrinkage:
            `phi_tilde_vox = w * phi_bar_ref + (1-w) * phi_hat_vox`
            where `w` is a shrinkage weight. `w` can be estimated based on the variance of `phi_hat_vox` within a region and the sampling variance of AR estimates (as sketched in user's proposal for AR(1) parcel-wise shrinkage).
    *   **Application:** Each voxel `v` is whitened using its specific `phi_tilde_vox[v]`.
    *   **Cost:**
        *   OLS residuals: Already computed.
        *   Raw `phi_hat_vox`: O(T\*V) for AR(1) (as per user sketch); O(p\*T\*V) for AR(p) if Yule-Walker per voxel.
        *   Shrinkage calculations: O(V), negligible.
        *   Whitening: Still O(T\*V).
    *   **Rationale:** Offers the most spatial specificity. Shrinkage is crucial to stabilize noisy per-voxel AR estimates. Best for high-resolution analyses or when fine-grained noise modeling is critical.

**A.3 Implementation Details & API Considerations**

*   **C++ Whitening Filter Adaptability:**
    The C++ `ar_whiten_inplace` function can be readily adapted. Instead of a scalar `phi` (for AR(1)) or a single `arma::vec phi_coeffs` (for AR(p)) applied to all voxels/columns, it would take an `arma::mat phi_matrix` (V x p or P x p) where each row `v` (or `j`) contains the AR coefficients for that voxel (or predictor).

    ```cpp
    // For voxel-specific AR(1) coefficients
    void ar1_whiten_voxel_specific_inplace(arma::mat& Y, arma::mat& X,
                                           const arma::vec& phi_per_voxel, /* ... */) {
        // ... omp loop over V voxels ...
        for (unsigned int v = 0; v < V; ++v) {
            const double a = -phi_per_voxel(v); // Voxel-specific coefficient
            double prev = 0.0;
            for (unsigned int t = 0; t < n_time; ++t) {
                double y_old = Y(t,v);
                Y(t,v) = y_old + a*prev;
                prev   = y_old;
            }
        }
        // ... similar loop for X, potentially using run-level phi for X ...
    }
    ```
    *Self-correction:* Whitening X columns (predictors) should ideally use a `phi` representative of the average noise structure that affects *all* voxels, so `phi_hat_run` (Tier 0) or `phi_hat_global` is likely most appropriate for X. Applying voxel-specific `phi` to X columns doesn't quite make sense, as X is common to all voxels. The design matrix X should be whitened using a single, globally (for the run) estimated AR coefficient set. The Y matrix is whitened using the (potentially) spatially varying AR coefficients.

*   **API Sketch Update:**
    ```R
    fmri_lm(...,
            cor_struct = c("iid", "ar1", "ar2", "arp"),
            ar_p       = NULL, # User-defined p for "arp"
            ar_scope   = c("run", "parcel", "voxel"), # New argument
            parcels    = NULL, # User-provided parcel map (NeuroVol or vector)
                               # or path to NIfTI file.
                               # If ar_scope="parcel" and parcels=NULL, can attempt
                               # internal 3-tissue segmentation (e.g. k-means on mean EPI).
            cor_iter   = 1L,
            cor_global = FALSE,
            ar1_exact_first = FALSE,
            ...)
    ```

*   **Workflow for Tiers in `runwise_lm` / `chunkwise_lm`:**
    1.  **Initial OLS:** Always perform initial OLS on original data to get `Resid_OLS_run`.
    2.  **Estimate `phi_hat_run` (Tier 0):** Always compute this from `rowMeans(Resid_OLS_run)`. This serves as a fallback and reference for shrinkage.
    3.  **If `ar_scope == "run"`:** Use `phi_hat_run` for whitening `Y_run_orig` and `X_run_orig`.
    4.  **If `ar_scope == "parcel"`:**
        *   If `parcels` is provided (or generated), compute `phi_hat_parcel` for each parcel in the run from `Resid_OLS_run` (parcel-mean residuals).
        *   Create `phi_vector_for_Y` (length V) by mapping `phi_hat_parcel` to voxels.
        *   Whiten `Y_run_orig` using `phi_vector_for_Y`.
        *   Whiten `X_run_orig` using `phi_hat_run` (or global mean parcel phi).
    5.  **If `ar_scope == "voxel"`:**
        *   Compute raw `phi_hat_vox` for all voxels in the run from `Resid_OLS_run`.
        *   Determine reference `phi_bar_ref` (either `phi_hat_run` or `phi_hat_parcel` if parcels available).
        *   Compute shrinkage weights `w` and shrunken `phi_tilde_vox`.
        *   Create `phi_vector_for_Y` (length V) using `phi_tilde_vox`.
        *   Whiten `Y_run_orig` using `phi_vector_for_Y`.
        *   Whiten `X_run_orig` using `phi_hat_run` (or global mean shrunken phi).
    6.  **Proceed with GLS fit** on whitened data. Iteration (`cor_iter > 1`) would re-estimate the chosen scope of `phi` from GLS residuals.

**A.4 Performance:**

*   **Tier 1 (Parcel):** Adds cost of parcel-mean residual calculation (O(T\*V) once) and per-parcel Yule-Walker (negligible). Total impact on runtime should be small.
*   **Tier 2 (Voxel):** Adds cost of per-voxel raw AR estimation. For AR(1), this is O(T\*V) (as shown in user sketch, comparable to one `var()` pass). For AR(p), if Yule-Walker is per-voxel, it's O(p\*T\*V). Shrinkage math is O(V). This is the most computationally intensive adaptive tier but still linear in V and T.

**A.5 User Guidance & Defaults:**

*   Default remains `ar_scope = "run"` (Tier 0) for maximum speed and good general performance.
*   `ar_scope = "parcel"` is a good intermediate if strong tissue-specific noise is suspected and a parcel map is available.
*   `ar_scope = "voxel"` is for advanced use cases needing maximal spatial adaptivity, with a moderate performance cost.
*   Documentation should clearly state that `X` (design matrix) is typically whitened with a run-level or global AR estimate, not voxel-specific ones.

**A.6 Ticketed Sprint Additions (for Locally Adaptive AR):**

This would follow the initial AR(p) implementation.

*   **Ticket LA.1: API Extension for `ar_scope` and `parcels` (R)**
    *   Add `ar_scope` and `parcels` arguments to `fmri_lm` and thread them down.
    *   Implement logic for basic 3-tissue segmentation (e.g., k-means on mean EPI) if `ar_scope="parcel"` and `parcels=NULL`. (Could be deferred if too complex initially, requiring user-provided parcels).
    *   Estimate: 0.5 - 1 day (depending on auto-segmentation complexity).

*   **Ticket LA.2: C++ Voxel-Specific AR(1) Estimation & Whitening**
    *   Implement C++ function for fast per-voxel AR(1) `phi_hat_vox` estimation from residuals.
    *   Modify `ar_whiten_inplace` (or create new version) to accept a vector `phi_per_voxel` for whitening Y, while still using a single `phi_run` for X.
    *   Extend to AR(p) if feasible within initial scope (might require per-voxel Yule-Walker or a simpler approximation).
    *   Unit test.
    *   Estimate: 1.5 - 2 days.

*   **Ticket LA.3: Tier 1 (Parcel-Level) AR Logic (R)**
    *   Implement logic in `runwise_lm`/`chunkwise_lm` to:
        *   Calculate parcel-mean OLS residuals.
        *   Estimate `phi_hat_parcel` for each parcel.
        *   Construct `phi_vector_for_Y` mapping parcel phis to voxels.
        *   Call C++ whitener with `phi_vector_for_Y` (for Y) and `phi_hat_run` (for X).
    *   Estimate: 1 day.

*   **Ticket LA.4: Tier 2 (Voxel-Level with Shrinkage) AR Logic (R)**
    *   Implement logic in `runwise_lm`/`chunkwise_lm` to:
        *   Call C++ for `phi_hat_vox`.
        *   Implement James-Stein-like shrinkage towards `phi_hat_run` or `phi_hat_parcel`.
        *   Construct `phi_tilde_vox` vector.
        *   Call C++ whitener with `phi_tilde_vox` (for Y) and `phi_hat_run` (for X).
    *   Estimate: 1.5 days.

*   **Ticket LA.5: Testing and Documentation for Locally Adaptive AR**
    *   Add unit tests for parcel and voxel scopes.
    *   Update documentation and vignettes to explain new scopes and `parcels` argument.
    *   Estimate: 1 day.

This tiered approach offers a flexible and performant way to handle spatially varying serial correlations, building incrementally on the foundational run-level AR prewhitening.