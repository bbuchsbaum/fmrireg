This is an *excellent* set of guidelines for code hygiene, modularity, and maintainability. It's exactly what's needed to ensure `fmrireg` remains robust and developer-friendly as it grows. The "TL;DR for Developers" is a perfect summary.

Let's integrate these principles into a revised, comprehensive proposal and ticketed sprint. The existing "Phase 1, 2, 3" structure will be maintained, but the new ARCH tickets will be prioritized as foundational.

---

**Project: Integrated Robust & AR(p) Modeling with Architectural Refinement (Version 4.0)**

**Goal:** Deliver a robust, efficient, user-friendly, and maintainable implementation of fMRI linear modeling in `fmrireg`. This version focuses on integrating Iteratively Reweighted Least Squares (IRLS) with Autoregressive (AR(p)) modeling, underpinned by a modular and clean codebase.

**Core Design & Architectural Principles:**

1.  **Primary Fitting Sequence ("Whiten then Robustly Weight"):**
    *   Optional: Regress out `extra_nuisance` regressors.
    *   Initial OLS/GLS to estimate AR parameters (`phi_hat`).
    *   AR Pre-whitening of data (`Y`) and design (`X`).
    *   IRLS on the whitened data (`Y_w`, `X_w`).
    *   Optional: Re-estimate `phi_hat` and perform a final weighted GLS.
2.  **Modularity (Slice by Responsibility):**
    *   `R/fmri_lm_config.R`: Configuration object (`fmri_lm_config`) creation and validation.
    *   `R/fmri_lm_context.R`: GLM context object (`glm_context`) definition.
    *   `R/fmri_lm_solver.R`: Core GLM solver (`solve_glm_core`) for OLS/WLS.
    *   `R/fmri_ar_modeling.R`: AR parameter estimation (`estimate_ar_parameters`) and data whitening (`ar_whiten_transform`).
    *   `R/fmri_robust_fitting.R`: IRLS engine (`robust_iterative_fitter`) using `solve_glm_core` and `ar_whiten_transform`.
    *   `R/fmri_lm_orchestrators.R`: `runwise_fitter` and `chunkwise_fitter` orchestrating the steps.
    *   `R/fmrilm.R`: Top-level `fmri_lm` and `fmri_lm_fit` functions.
3.  **Minimized Surface Area:** Use `fmri_lm_config` for options and `glm_context` for data transfer between modules.
4.  **Single Source of Truth for Math:** Centralize core matrix operations in `solve_glm_core`.
5.  **CI Guardrails:** Implement `lintr`, `styler`, and code size checks.
6.  **Progressive Disclosure in Docs:** Clear separation of user API and internal engine documentation.
7.  **Encapsulated Configuration:** `fmri_lm_config` object.
8.  **Semantic Tests:** Small, focused tests per module and integration tests.
9.  **Centralized Error Handling:** Utility functions for common validation/error messages.

---

**API Changes (`fmri_lm`):**

*   `robust`: `c(FALSE, "huber", "bisquare")`. Default `FALSE`.
*   `robust_options`: A list or dedicated S3 object (e.g., `robust_control()`) for `k_huber`, `c_tukey`, `max_iter`, `scale_scope`, `reestimate_phi`. Default `NULL` (uses internal defaults).
*   `ar_options`: A list or dedicated S3 object (e.g., `ar_control()`) for `cor_struct`, `cor_iter` (for non-robust), `cor_global`, `ar_p`, `ar1_exact_first`. Default `NULL`.
*   `extra_nuisance`: `NULL`, matrix, or formula.
*   `keep_extra_nuisance_in_model`: `FALSE`.
*   `ar_voxelwise`: `FALSE`.
*   *Removed:* `robust_psi`, `robust_k_huber`, `robust_c_tukey`, `robust_max_iter`, `robust_scale_scope`, `cor_struct`, `cor_iter`, `cor_global`, `ar_p`, `ar1_exact_first` as top-level args.

---

**Immediate "Must Fix" Items (Pre-Sprint):**

*   **Ticket MUST-FIX-001R: Remove `robust && use_fast_path` Blocker**
    *   **Task:** Delete `if (robust && use_fast_path) robust <- FALSE` in `fmri_lm()` and `fmri_lm_fit()`.
    *   **Reason:** Enables new robust fast path.

*   **Ticket MUST-FIX-002R: Explicit `NA` Checks**
    *   **Task:** Add `if (anyNA(Y) || anyNA(X)) stop(...)` in core whitening and robust fitting engines before matrix operations commence on data expected to be NA-free.
    *   **Reason:** Prevent silent errors.

*   **Ticket MUST-FIX-003R: Correct Pooled/Run-Specific Sigma in Robust Paths**
    *   **Task:** Ensure `sigma_robust` (scalar per run or global) is correctly calculated by `robust_iterative_fitter` (new name for `fast_rlm_run`) and correctly used for SE calculation in `beta_stats_matrix` and `fit_lm_contrasts_fast`. Avoid incorrect averaging when pooling.
    *   **Reason:** Statistical correctness of SEs.

---

**Ticketed Sprint: Integrated Robust & AR(p) Modeling with Architectural Refinement**

**Phase 0: Architectural Foundation (Blockers for subsequent work)**

*   **Ticket ARCH-001: Implement `fmri_lm_config` Object & Factory** *(implemented)*
    *   New `fmri_lm_control()` factory now lives in `R/fmri_lm_config.R` and returns
        an `fmri_lm_config` object containing validated `robust` and `ar` option
        lists.
    *   **Task:** Create `R/fmri_lm_config.R`. Define `fmri_lm_config` S3 class. Implement `fmri_lm_control(robust_options = list(...), ar_options = list(...), ...)` factory function that takes all relevant fitting options, applies defaults, validates, and returns an `fmri_lm_config` object.
    *   **Details:** `robust_options` list to contain `type`, `k_huber`, `c_tukey`, `max_iter`, `scale_scope`, `reestimate_phi`. `ar_options` list to contain `struct`, `p`, `iter_gls`, `global`, `voxelwise`, `exact_first`.
    *   **Acceptance:** Config object created, validated, and holds all fitting parameters.

*   **Ticket ARCH-002: Implement `glm_context` Object**
    *   **Task:** Create `R/fmri_lm_context.R`. Define `glm_context` S3 class (a list) to hold `X`, `Y`, `proj` (from `.fast_preproject`), `phi_hat`, `sigma_robust_scale`, `robust_weights`.
    *   **Acceptance:** Context object defined.

*   **Ticket ARCH-003: Create Core Solver `solve_glm_core`**
    *   **Task:** Create `R/fmri_lm_solver.R`. Implement `solve_glm_core(glm_ctx, return_fitted = FALSE)`.
    *   **Details:**
        *   If `glm_ctx$robust_weights` is `NULL`, performs OLS/GLS using `glm_ctx$X`, `glm_ctx$Y`, and `glm_ctx$proj`.
        *   If `glm_ctx$robust_weights` is not `NULL`, it assumes `X` and `Y` in context are *already weighted* (i.e., `X_w`, `Y_w`). It then performs WLS using `glm_ctx$proj` (which should be `proj_w`).
        *   Replaces current `.fast_lm_matrix()`. Output: list with `betas`, `rss`, `sigma2`, `fitted` (if requested).
    *   **Acceptance:** Single, unified low-level solver for OLS/WLS.

*   **Ticket ARCH-004: Modularize AR Functions**
    *   **Task:** Create `R/fmri_ar_modeling.R`. Move `estimate_ar_parameters()` and `ar_whiten_transform()` here. Ensure they are robust.
    *   **Depends on:** (none in this phase)
    *   **Acceptance:** AR utilities are self-contained.

*   **Ticket ARCH-005: Implement Robust Engine `robust_iterative_fitter`** *(implemented)*
    *   **Task:** Create `R/fmri_robust_fitting.R`. Implement `robust_iterative_fitter(initial_glm_ctx, cfg_robust_options, X_orig_for_resid, sigma_fixed = NULL)`.
    *   **Details:** This function takes an *initial* `glm_context` (typically after OLS on original or whitened data). It performs the IRLS loop:
        1.  Calculate residuals (needs `X_orig_for_resid` which is the `X` corresponding to the `Y` in `initial_glm_ctx` *before* robust weighting).
        2.  Estimate `sigma_robust_scale` (using `sigma_fixed` if provided and `cfg_robust_options$scale_scope` is global, else from current residuals).
        3.  Calculate `robust_weights`.
        4.  Create `glm_ctx_weighted` with `X_w = X_orig_for_resid * sqrt(robust_weights)`, `Y_w = initial_glm_ctx$Y * sqrt(robust_weights)`, and `proj_w = .fast_preproject(X_w)`.
        5.  `fit_wls <- solve_glm_core(glm_ctx_weighted)`. Update `betas`.
        6.  Check convergence. Loop or exit.
    *   **Output:** A list containing final `betas_robust`, `XtWXi_final = proj_w$XtXinv` (from last weighted iteration), `sigma_robust_scale_final`, `robust_weights_final`, `dfres`.
    *   **Depends on:** ARCH-002, ARCH-003.
    *   **Acceptance:** IRLS engine implemented.

*   **Ticket ARCH-006: Refactor `fmri_lm` and `fmri_lm_fit` for Config Object**
    *   **Task:** Modify `fmri_lm` to accept `robust_options`, `ar_options`, etc. Create `cfg <- fmri_lm_control(...)` inside `fmri_lm`. `fmri_lm_fit` now takes `cfg` instead of many individual arguments.
    *   **Depends on:** ARCH-001.
    *   **Acceptance:** API simplified. Config object used internally.

*   **Ticket ARCH-007: Lintr/Styler CI + Size Guard**
    *   **Task:** Set up GitHub Action.
    *   **Acceptance:** CI enforces style and size limits.

**Phase 1: Fast Path OLS/GLS and Robust-Only (No Combined AR+Robust Yet)**

*   **Ticket SPRINT3-01R: Fast Path Standard OLS/GLS (Non-Robust)**
    *   **Task:** Refactor `runwise_lm` and `chunkwise_lm` fast paths for `cfg$robust$type == FALSE`.
    *   **`runwise_lm`:**
        1.  For each run, create `glm_ctx_run_orig` (`X_run`, `Y_run`, `proj_run = .fast_preproject(X_run)`).
        2.  If `cfg$ar$struct != "iid"`:
            *   Iteratively estimate `phi_hat_run` (or use `phi_global`), whiten `X_run`, `Y_run` into `glm_ctx_run_whitened`. Use `solve_glm_core` on this context.
        3.  Else (iid): Use `solve_glm_core` on `glm_ctx_run_orig`.
        4.  Collect results.
    *   **`chunkwise_lm`:**
        1.  Precompute `proj_global_orig = .fast_preproject(X_global_orig)`.
        2.  If `cfg$ar$struct != "iid"`: Precompute `phi_hat_run`s (or `phi_global`). Precompute `X_global_w` and `proj_global_w = .fast_preproject(X_global_w)`.
        3.  For each chunk: Whiten `Y_chunk` using appropriate `phi_hat_run`s to get `Y_chunk_w`. Call `solve_glm_core` with `X_global_w`, `Y_chunk_w`, `proj_global_w`.
        4.  Else (iid): Call `solve_glm_core` with `X_global_orig`, `Y_chunk_orig`, `proj_global_orig`.
    *   **Depends on:** ARCH-001 to ARCH-006.
    *   **Acceptance:** Fast OLS and GLS paths work correctly using new modular components.

*   **Ticket SPRINT3-02R: Fast Path Robust-Only (No AR)**
    *   **Task:** Refactor `runwise_lm` and `chunkwise_lm` fast paths for `cfg$robust$type != FALSE` and `cfg$ar$struct == "iid"`.
    *   **`runwise_lm`:**
        1.  For each run: `glm_ctx_run_orig` (`X_run`, `Y_run`, `proj_run`).
        2.  `robust_fit_run <- robust_iterative_fitter(glm_ctx_run_orig, cfg$robust, X_run, sigma_fixed = sigma_global_if_scope_global)`.
        3.  Collect results from `robust_fit_run`.
    *   **`chunkwise_lm` (using "fully pre-weighted" strategy):**
        1.  **Pass-0 (Per Run):** For each run `r`, create `glm_ctx_run_r_orig`. Call `robust_iterative_fitter` to get `w_robust_run_r` and `sigma_robust_run_r`.
        2.  **Global Weighted Matrices:** Form `X_global_robustly_weighted` and `Y_global_robustly_weighted` by applying `sqrt(w_robust_run_r)` to corresponding run segments of `X_global_orig` and `Y_global_orig`.
        3.  `proj_global_robustly_weighted <- .fast_preproject(X_global_robustly_weighted)`.
        4.  **Chunk Processing:** For each chunk of `Y_global_robustly_weighted` (call it `Y_chunk_rw`), use `solve_glm_core(list(X=X_global_robustly_weighted, Y=Y_chunk_rw, proj=proj_global_robustly_weighted))`.
        5.  Stats use `sigma_robust_run_r` (mapped to voxels) and `proj_global_robustly_weighted$dfres`.
    *   **Depends on:** ARCH-001 to ARCH-006.
    *   **Acceptance:** Fast robust-only fitting works correctly.

**Phase 2: Combined AR + Robust Fast Paths & Advanced AR**

*   **Ticket SPRINT3-03R: Fast Path AR + Robust (`runwise_lm`)**
    *   **Task:** Implement the "Whiten then Robustly Weight" pipeline for `runwise_lm` fast path.
        1.  For each run:
            *   Initial OLS on `X_run_orig`, `Y_run_orig` to get residuals (consider `extra_nuisance`).
            *   Estimate `phi_hat_run` (or use `phi_global` if `cfg$ar$global`).
            *   `X_run_w, Y_run_w <- ar_whiten_transform(X_run_orig, Y_run_orig, phi_hat_run, cfg$ar$exact_first)`.
            *   `proj_run_w <- .fast_preproject(X_run_w)`.
            *   `glm_ctx_run_whitened <- list(X=X_run_w, Y=Y_run_w, proj=proj_run_w)`.
            *   `robust_fit_run <- robust_iterative_fitter(glm_ctx_run_whitened, cfg$robust, X_orig_for_resid = X_run_w, sigma_fixed = sigma_global_if_scope_global)`.
            *   **Optional `robust_reestimate_phi="once"`:**
                *   De-whiten `robust_fit_run$residuals_final_robust_weighted` (using `phi_hat_run`).
                *   `phi_hat_run_updated <- estimate_ar_parameters(...)`.
                *   `X_run_w2, Y_run_w2 <- ar_whiten_transform(X_run_orig, Y_run_orig, phi_hat_run_updated, ...)`.
                *   Create `glm_ctx_final_wls` with `X_run_w2 * sqrt(robust_fit_run$weights_final)`, `Y_run_w2 * sqrt(robust_fit_run$weights_final)`, and its projection.
                *   `final_fit <- solve_glm_core(glm_ctx_final_wls)`. Use this for results.
        2.  Pool results.
    *   **Depends on:** Phase 1 tickets, ARCH-004, ARCH-005.
    *   **Acceptance:** `runwise_lm` handles AR+Robust via fast path, including `robust_reestimate_phi`.

*   **Ticket SPRINT3-04R: Fast Path AR + Robust (`chunkwise_lm`)**
    *   **Task:** Implement "Whiten then Robustly Weight" for `chunkwise_lm` fast path. This is complex.
        1.  **Pass-0 (Per Run Precomputation):**
            *   For each run `r`:
                *   Estimate `phi_hat_run` (or global) from initial OLS on `(X_run_orig, Y_run_full)` (after `extra_nuisance`).
                *   `X_run_w, Y_run_full_w <- ar_whiten_transform(X_run_orig, Y_run_orig_full_run, phi_hat_run, ...)`.
                *   `proj_run_w <- .fast_preproject(X_run_w)`.
                *   `robust_fit_details_run <- robust_iterative_fitter(list(X=X_run_w, Y=Y_run_full_w, proj=proj_run_w), cfg$robust, X_orig_for_resid=X_run_w, sigma_fixed=sigma_global_if_scope_global)`.
                *   Store `phi_hat_run`, `w_robust_run = robust_fit_details_run$weights_final`, `sigma_robust_run = robust_fit_details_run$sigma_robust_scale_final`.
                *   (Optional `robust_reestimate_phi="once"`: if so, update `phi_hat_run` here based on de-whitened robust residuals of this full run fit).
        2.  **Global Transformed Matrices Construction:**
            *   `X_global_final_w <- matrix(...)`, `Y_global_final_w <- matrix(...)`.
            *   For each run segment: apply `phi_hat_run` to get `X_seg_w, Y_seg_w`. Then apply `sqrt(w_robust_run)` to get `X_seg_wr, Y_seg_wr`. Concatenate these.
        3.  `proj_global_final_w <- .fast_preproject(X_global_final_w)`.
        4.  **Chunk Processing:** For each chunk of `Y_global_final_w` (call it `Y_chunk_fwr`), use `solve_glm_core(list(X=X_global_final_w, Y=Y_chunk_fwr, proj=proj_global_final_w))`.
        5.  Stats calculation uses run-specific `sigma_robust_run` (mapped to voxels) and `proj_global_final_w$dfres`.
    *   **Depends on:** Phase 1 tickets, ARCH-004, ARCH-005, SPRINT3-03R (for `robust_reestimate_phi` pattern).
    *   **Acceptance:** `chunkwise_lm` AR+Robust fast path operational. This is the most challenging ticket.

*   **Ticket SPRINT3-05R: Implement `ar_voxelwise` (Slow Path Only)**
    *   **Task:** Implement as previously described in `runwise_lm` (slow path `!use_fast_path`).
    *   **Depends on:** ARCH-004.
    *   **Acceptance:** Voxel-wise AR works in slow path.

**Phase 3: Final Touches, Documentation & Testing**

*   **Ticket SPRINT3-06R: Effective Degrees of Freedom & Sandwich Docs**
    *   **Task:** Implement `calculate_effective_df` and integrate into reporting as previously described.
    *   **Acceptance:** Effective DF available. Documentation clear.

*   **Ticket SPRINT3-07R: Comprehensive Documentation**
    *   **Task:** Update all relevant documentation for `fmri_lm` and new internal functions. Explain new config objects (`fmri_lm_control`).
    *   **Acceptance:** User documentation complete.

*   **Ticket SPRINT3-08R: Extensive Testing & CI**
    *   **Task:**
        *   Add the specific CI test for `fast_path=TRUE` vs `FALSE` for AR+Robust.
        *   Test all new API options, config objects, and their interactions.
        *   Validate `robust_reestimate_phi` and `ar_voxelwise`.
        *   Test error handling and `NA` propagation guards.
    *   **Acceptance:** High test coverage, CI passes.

---

This revised plan heavily emphasizes the architectural changes first (Phase 0), then builds the fast paths incrementally (Robust-only, then AR-only, then combined AR+Robust). The `chunkwise` AR+Robust path remains the most intricate. The API is simplified by grouping options into `fmri_lm_control`. The "must-fix" items are critical prerequisites.

---

**Phase 4: Voxelwise AR Contrast Support**

**Problem Statement:**
The current voxelwise AR implementation in the slow path (SPRINT3-05R) does not compute contrasts. It returns an empty contrast list with a comment "would need proper handling of contrasts". This is a critical gap that prevents users from performing hypothesis testing when using voxelwise AR modeling.

**Technical Challenges:**
1. Each voxel has different AR parameters, leading to different whitening transformations
2. The (X'X)^-1 matrix differs for each voxel after whitening
3. Standard errors must account for voxel-specific whitening
4. Memory efficiency is crucial when storing per-voxel covariance matrices

**Proposed Solution:**

*   **Ticket SPRINT4-01R: Refactor Voxelwise AR to Use Modular Components**
    *   **Task:** Replace the current `lm.fit()` calls with the modular `glm_context` and `solve_glm_core` approach.
    *   **Details:**
        1. For each voxel:
            - Create initial `glm_ctx_voxel` with `X_run` and single-voxel `Y`
            - Estimate voxel-specific `phi_voxel`
            - Create whitened `glm_ctx_voxel_w` after `ar_whiten_transform`
            - Call `solve_glm_core` to get betas
            - Store the whitened `proj$XtXinv` for contrast calculations
        2. This ensures consistency with other paths and enables contrast support
    *   **Depends on:** ARCH-002, ARCH-003, ARCH-004
    *   **Acceptance:** Voxelwise AR uses same solver as other paths

*   **Ticket SPRINT4-03R: Choose and Implement Voxelwise Strategy** *(implemented)*
    *   **Task:** Implement the chosen memory-efficient approach for voxelwise contrasts.
    *   **Details:**
        1. **Chunked Processing Approach** *(chosen)*
           - Minimal memory footprint (< 100 MB)
           - Slightly slower but negligible for most analyses
           - Aligns with fmrireg's existing chunking philosophy
        
        2. **Alternative: QR Storage**
           - If speed is critical and memory allows (~500MB-1GB)
           - Better for iterative/interactive analyses
           - Store QR decompositions instead of covariance matrices
        
        3. **Implementation priorities:**
           - Correctness first
           - Memory efficiency second
           - Speed optimization third
           - Parallel processing fourth
    *   **Implementation:** Voxels are processed in small chunks (default size 100),
        whitening and computing contrasts within each chunk. This keeps memory
        usage low while maintaining accuracy.
    *   **Depends on:** SPRINT4-01R
    *   **Acceptance:** Whole-brain voxelwise analysis completes without memory errors

*   **Ticket SPRINT4-02R: Implement Memory-Efficient Voxelwise Contrast Engine**
    *   **Task:** Create contrast calculation system that avoids storing per-voxel covariance matrices.
    *   **Details:**
        **Approach 1 - Store QR Decompositions (Recommended):**
        ```r
        # For each voxel, store compact QR instead of (X'X)^-1
        # QR storage: n×p instead of p×p (typically 3-5x more efficient)
        
        fit_lm_contrasts_voxelwise_qr <- function(betas_matrix, qr_list, sigma_vec,
                                                   conlist, fconlist, dfres) {
          n_voxels <- ncol(betas_matrix)
          
          # Process each contrast
          contrast_results <- list()
          
          for (con_name in names(conlist)) {
            con_spec <- conlist[[con_name]]
            colind <- attr(con_spec, "colind")
            
            # Allocate result vectors
            est_vec <- numeric(n_voxels)
            se_vec <- numeric(n_voxels)
            
            # Batch process voxels with similar QR structures
            for (v in seq_len(n_voxels)) {
              # Use QR to solve for contrast variance efficiently
              qr_v <- qr_list[[v]]
              l_full <- numeric(ncol(qr_v$qr))
              l_full[colind] <- con_spec
              
              # Efficient computation using QR
              # var(l'beta) = sigma^2 * l'(X'X)^-1*l = sigma^2 * ||Q'l||^2
              Qtl <- qr.qty(qr_v, l_full)
              var_contrast <- sum(Qtl[1:qr_v$rank]^2) * sigma_vec[v]^2
              
              est_vec[v] <- sum(l_full * betas_matrix[, v])
              se_vec[v] <- sqrt(var_contrast)
            }
            
            # Compute statistics
            t_vec <- est_vec / se_vec
            p_vec <- 2 * pt(abs(t_vec), dfres, lower.tail = FALSE)
            
            contrast_results[[con_name]] <- create_contrast_tibble(
              con_name, est_vec, se_vec, t_vec, p_vec, sigma_vec
            )
          }
          return(contrast_results)
        }
        ```
        
        **Approach 2 - Chunked Processing (Memory-Constrained):**
        ```r
        # Process voxels in chunks, recomputing XtXinv as needed
        fit_lm_contrasts_voxelwise_chunked <- function(X_run, Y_run, phi_matrix,
                                                        conlist, fconlist, 
                                                        chunk_size = 100) {
          n_voxels <- ncol(Y_run)
          n_chunks <- ceiling(n_voxels / chunk_size)
          
          # Pre-allocate storage for all contrasts
          contrast_storage <- initialize_contrast_storage(conlist, fconlist, n_voxels)
          
          for (chunk_idx in seq_len(n_chunks)) {
            voxel_idx <- ((chunk_idx-1)*chunk_size + 1):min(chunk_idx*chunk_size, n_voxels)
            
            # Process chunk of voxels
            for (v_local in seq_along(voxel_idx)) {
              v_global <- voxel_idx[v_local]
              
              # Whiten for this voxel
              phi_v <- phi_matrix[, v_global]
              tmp <- ar_whiten_transform(X_run, Y_run[, v_global, drop=FALSE], 
                                         phi_v, exact_first = TRUE)
              X_w <- tmp$X
              Y_w <- tmp$Y
              
              # Compute what we need for contrasts
              qr_w <- qr(X_w)
              beta_w <- qr.coef(qr_w, Y_w)
              sigma_w <- sqrt(sum(qr.resid(qr_w, Y_w)^2) / (nrow(X_w) - qr_w$rank))
              
              # Calculate all contrasts for this voxel
              store_voxel_contrasts(contrast_storage, v_global, 
                                    qr_w, beta_w, sigma_w, conlist, fconlist)
            }
          }
          
          return(format_contrast_results(contrast_storage))
        }
        ```
        
        **Memory comparison:**
        - Original: p×p×V matrices = 40×40×100000×8 bytes = 1.28 GB
        - QR approach: n×p×V = 200×40×100000×8 bytes = 640 MB (2x savings)
        - Chunked: p×p×chunk_size = 40×40×100×8 bytes = 1.28 MB (1000x savings)
        
    *   **Depends on:** SPRINT4-01R
    *   **Acceptance:** Memory usage stays under 1GB for typical whole-brain analysis

*   **Ticket SPRINT4-04R: Integrate Voxelwise Contrasts into runwise_lm**
    *   **Task:** Update the voxelwise AR branch in `runwise_lm` to call the new contrast engine.
    *   **Details:**
        1. Replace empty `ret$contrasts <- list()` with:
           ```r
           conres <- fit_lm_contrasts_voxelwise(
             betas_voxelwise, 
             XtXinv_voxel_list,
             sigma_voxelwise,
             simple_conlist_weights,
             fconlist_weights,
             rdf
           )
           ```
        2. Update `beta_stats_matrix` call to handle voxelwise case
        3. Ensure proper structure for downstream pooling
    *   **Depends on:** SPRINT4-02R
    *   **Acceptance:** Voxelwise AR returns valid contrasts

*   **Ticket SPRINT4-05R: Add Voxelwise AR + Robust Support**
    *   **Task:** Implement full IRLS for voxelwise AR + robust combination.
    *   **Details:**
        1. For each voxel:
            - Initial OLS → estimate AR parameters
            - Whiten data with voxel-specific phi
            - Run IRLS on whitened data
            - Store robust weights and final estimates
        2. Use `robust_iterative_fitter` with single-voxel contexts
        3. Calculate contrasts with both AR and robust adjustments
    *   **Depends on:** SPRINT4-01R, ARCH-005
    *   **Acceptance:** Voxelwise AR + Robust fully functional

*   **Ticket SPRINT4-06R: Performance Optimization**
    *   **Task:** Optimize voxelwise calculations for realistic datasets.
    *   **Details:**
        1. Implement parallel processing over voxels (using future/parallel)
        2. Vectorize where possible (batch similar AR parameters)
        3. Add progress reporting for long-running analyses
        4. Consider C++ implementation for inner loops if needed
    *   **Depends on:** SPRINT4-04R
    *   **Acceptance:** Voxelwise processing time is practical (< 10x slower than non-voxelwise)

*   **Ticket SPRINT4-07R: Comprehensive Testing**
    *   **Task:** Add extensive tests for voxelwise AR with contrasts.
    *   **Details:**
        1. Test contrast accuracy by comparing single-voxel results to standard lm()
        2. Test memory efficiency with large numbers of voxels
        3. Test edge cases (constant voxels, singular designs after whitening)
        4. Test pooling of voxelwise results across runs
        5. Verify effective DF calculations include voxel-specific AR orders
    *   **Depends on:** SPRINT4-04R, SPRINT4-05R
    *   **Acceptance:** All tests pass, edge cases handled

**Alternative Approach (if memory is critical):**

Instead of storing XtXinv for each voxel, we could:
1. Compute a "reference" XtXinv using average AR parameters
2. Store only the deviation of each voxel's XtXinv from reference
3. Use perturbation theory to approximate voxel-specific standard errors

This would trade some accuracy for substantial memory savings.

---

**Phase 5: Code Modularization and Cleanup**

**Problem Statement:**
The `fmrilm.R` file has grown to over 2000 lines with mixed responsibilities, duplicated code, and strategy implementations that are nearly 1000 lines each. This violates the modularity principles and makes the code hard to maintain.

*   **Ticket SPRINT5-01R: Extract Model Utilities Module**
    *   **Task:** Create `R/fmri_model_utils.R` and move model-related utilities.
    *   **Details:**
        - Move: `get_formula.fmri_model`, `term_matrices.fmri_model`
        - Move: `create_fmri_model` (consider renaming to `fmri_model_from_formula`)
        - Keep these functions focused on model manipulation, not fitting
    *   **Acceptance:** Clean separation of model creation from fitting logic

*   **Ticket SPRINT5-02R: Extract Results/Methods Module**  
    *   **Task:** Create `R/fmri_lm_methods.R` for all S3 methods and result extraction.
    *   **Details:**
        - Move: `coef.fmri_lm`, `stats.fmri_lm`, `standard_error.fmri_lm`, `print.fmri_lm`
        - Move: `fitted_hrf.fmri_lm`, `pull_stat`, `pull_stat_revised`
        - Move: `reshape_coef` and other result manipulation helpers
    *   **Acceptance:** All S3 methods in one place, separate from fitting logic

*   **Ticket SPRINT5-03R: Create Strategy Base Module**
    *   **Task:** Create `R/fmri_lm_strategies.R` with shared strategy code.
    *   **Details:**
        ```r
        # Extract common patterns into reusable functions:
        
        # Process a single run (used by both strategies)
        process_run_standard <- function(run_data, model, cfg, phi_fixed = NULL) {
          # Common OLS/GLS logic
        }
        
        process_run_robust <- function(run_data, model, cfg, phi_fixed = NULL, sigma_fixed = NULL) {
          # Common robust fitting logic
        }
        
        process_run_ar_robust <- function(run_data, model, cfg, phi_fixed = NULL, sigma_fixed = NULL) {
          # Common AR+Robust logic
        }
        
        # Pooling utilities
        pool_run_results <- function(run_results) {
          # Common pooling logic
        }
        ```
    *   **Acceptance:** Shared logic extracted, no duplication between strategies

*   **Ticket SPRINT5-04R: Refactor runwise_lm Using Shared Components**
    *   **Task:** Rewrite `runwise_lm` to use the shared strategy components.
    *   **Details:**
        1. Reduce from 500+ lines to ~200 lines
        2. Use `process_run_*` functions from strategy base
        3. Clear separation of slow/fast paths
        4. Move voxelwise AR to separate function `runwise_lm_voxelwise`
    *   **Depends on:** SPRINT5-03R
    *   **Acceptance:** Cleaner, more maintainable runwise implementation

*   **Ticket SPRINT5-05R: Refactor chunkwise_lm Using Shared Components**
    *   **Task:** Rewrite `chunkwise_lm` to use the shared strategy components.
    *   **Details:**
        1. Reduce from 945 lines to ~300 lines
        2. Extract pre-computation phase to `prepare_chunkwise_matrices`
        3. Extract chunk processing to `process_chunk`
        4. Use shared pooling utilities
    *   **Depends on:** SPRINT5-03R
    *   **Acceptance:** Cleaner, more maintainable chunkwise implementation

*   **Ticket SPRINT5-06R: Create Internal Utilities Module**
    *   **Task:** Create `R/fmri_lm_internal.R` for low-level utilities.
    *   **Details:**
        - Move: `.fast_preproject`, `.fast_lm_matrix` (if still needed)
        - Move: `is.formula` and other small helpers
        - Document these as internal with proper roxygen2 tags
        - Consider which should be exported for advanced users
    *   **Acceptance:** Clean separation of internal utilities

*   **Ticket SPRINT5-07R: Update Imports and Exports**
    *   **Task:** Ensure all new modules properly import/export functions.
    *   **Details:**
        1. Update NAMESPACE via roxygen2
        2. Add necessary @importFrom statements  
        3. Document internal functions with @keywords internal
        4. Run R CMD check to verify no missing imports
    *   **Depends on:** SPRINT5-01R through SPRINT5-06R
    *   **Acceptance:** Package builds cleanly with new module structure

*   **Ticket SPRINT5-08R: Integration Testing**
    *   **Task:** Verify that modularization doesn't break existing functionality.
    *   **Details:**
        1. Run full test suite
        2. Add integration tests that verify module boundaries
        3. Check that all examples still work
        4. Benchmark to ensure no performance regression
    *   **Depends on:** SPRINT5-07R
    *   **Acceptance:** All tests pass, no performance degradation

**Expected Outcome:**
- `fmrilm.R` reduced from 2000+ lines to ~400 lines (just core API)
- Clear module structure:
  - `fmri_model_utils.R` (~200 lines)
  - `fmri_lm_methods.R` (~300 lines)
  - `fmri_lm_strategies.R` (~400 lines)
  - `fmri_lm_runwise.R` (~300 lines)
  - `fmri_lm_chunkwise.R` (~400 lines)
  - `fmri_lm_internal.R` (~150 lines)
- Easier to maintain, test, and extend
- Follows single-responsibility principle