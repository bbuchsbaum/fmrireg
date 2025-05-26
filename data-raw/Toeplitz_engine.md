# Plan for Toeplitz/FFT Engine Integration in fmrireg

## 1. Evaluation Summary

*   **Concept:** Using a Toeplitz matrix (or FFT) for fixed-kernel convolution is mathematically sound and leverages the LTI property.
*   **Efficiency:** Significant speed-up (5-25x+) is expected for long time series due to better algorithmic complexity (`O(N log N)` or optimized `O(N*L*K)`) compared to R loops or less optimized filtering.
*   **Tools:** `Matrix::bandSparse` is suitable for efficient sparse representation.
*   **Durations:** The cumulative sum trick (`+amp` at onset, `-amp` at offset using `cumsum(h_tr)*TR`) correctly handles non-zero durations.
*   **Multi-Basis HRFs:** Creating separate Toeplitz matrices (`T_b`) for each basis (`h_b`) and convolving the same impulse data (`S`) is correct.
*   **Integration:** Adding an `engine` argument to `event_model` for internal dispatching is a clean approach, maintaining API stability.
*   **Challenges:**
    *   Handling term-specific HRFs (the engine needs to loop through terms and use the *correct* HRF for each).
    *   Ensuring exact output structure match (column names, `term_spans`, `col_indices` attributes) with the classic engine.
    *   Correctly sampling/resampling HRF kernels to TR resolution.
*   **Feasibility:** The approach is feasible and sound.

## 2. Implementation Plan

1.  **Create New File:** `R/engine_toeplitz.R`
2.  **Implement `toeplitz_engine` function:**
    *   **Signature:** `toeplitz_engine(terms, sampling_frame, precision)`
    *   **Block Iteration:** Loop through blocks defined by `sampling_frame`.
    *   **Inner Block Loop (for block `k`, `N_tr` scans, `TR`):**
        *   Initialize `block_dm_parts = list()`.
        *   **Term Iteration:** Loop through `term` in `terms`:
            *   Get term-specific `hrf_obj = attr(term, "hrfspec")$hrf`.
            *   Get `nb = nbasis(hrf_obj)`.
            *   **Build Impulse Matrix `S`:** (`N_tr` rows, `N_cond` columns) for this `term` within this block.
                *   Use `design_matrix(term, drop.empty=FALSE)` for the *unconvolved* amplitudes per condition, restricted to the current block's events.
                *   Get term's `onsets`, `durations`, `blockids` for the current block.
                *   Map `onsets` to TR indices (`idx_on = floor(global_onset / TR) + 1`).
                *   Handle durations: Calculate offset indices (`idx_off = floor((global_onset + duration) / TR) + 1`). Create sparse matrix `S` using `Matrix::sparseMatrix` with `+amp` at `idx_on` and `-amp` at `idx_off` for each condition column. Clamp indices within `1:N_tr`.
            *   Initialize `term_basis_parts = list()`.
            *   **Basis Loop (`b` from 1 to `nb`):**
                *   Get `h_b`: Evaluate `hrf_obj$basis[[b]]` at fine `precision`.
                *   Resample/downsample `h_b` to TR resolution (`h_b_tr`). Determine kernel length `L`.
                *   Determine effective kernel `k_eff`: Use `cumsum(h_b_tr) * TR` if durations were used, otherwise use `h_b_tr`.
                *   Build sparse Toeplitz `T_b = Matrix::bandSparse(N_tr, N_tr, k = 0:(L-1), diagonals = list_of_diagonals_from(k_eff))`.
                *   Convolve: `X_b = as.matrix(T_b %*% S)`. (`N_tr` x `N_cond`).
                *   Store `X_b` in `term_basis_parts`.
            *   `term_dm_block = do.call(cbind, term_basis_parts)`.
            *   **Generate Column Names:** Create names for `term_dm_block` matching the classic engine (using `conditions(term)`, term name prefix, `_b` if `nb > 1`, `make.names(unique=TRUE)`).
            *   Store the named `term_dm_block` (as tibble) in `block_dm_parts`, using the unique term name as the key.
        *   `block_dm = do.call(cbind, block_dm_parts)`.
    *   **Combine Blocks:** Collect `block_dm` for all blocks in `final_dm_list`. Use `dplyr::bind_rows(final_dm_list)` to create `final_dm`.
    *   **Attach Attributes:** Re-calculate and attach `term_spans` and `col_indices` attributes to `final_dm`, matching the structure based on `attr(terms, "interface")`.
    *   **Return:** `final_dm`.
3.  **Modify `event_model` (`R/event_model.R`):**
    *   Add `engine = c("classic", "toeplitz")` argument (default `"classic"`).
    *   Add `engine <- match.arg(engine)`.
    *   Use `if (engine == "toeplitz") toeplitz_engine(...) else build_event_model_design_matrix(...)`.
4.  **Documentation:** Update `event_model` docs for `engine`. Add internal docs for `toeplitz_engine`.
5.  **Testing:** Add `testthat` tests in `tests/testthat/test_engine_toeplitz.R` comparing outputs and benchmarking.

## 3. TODO Items

*   **[TE-1]** Create `R/engine_toeplitz.R`.
*   **[TE-2]** Implement basic `toeplitz_engine` structure (signature, block loop).
*   **[TE-3]** Implement impulse matrix `S` generation within `toeplitz_engine`, handling term-specific conditions and block subsetting.
*   **[TE-4]** Implement duration handling (`cumsum` trick) for impulse matrix `S`.
*   **[TE-5]** Implement HRF kernel sampling and resampling to TR resolution (`h_b_tr`).
*   **[TE-6]** Implement sparse Toeplitz matrix construction (`T_b`) using `Matrix::bandSparse` and the effective kernel.
*   **[TE-7]** Implement the core convolution `T_b %*% S` and basis combination loop.
*   **[TE-8]** Implement robust column naming logic within `toeplitz_engine` to exactly match the classic engine's output.
*   **[TE-9]** Implement block combination using `dplyr::bind_rows`.
*   **[TE-10]** Implement the calculation and attachment of `term_spans` and `col_indices` attributes.
*   **[TE-11]** Modify `event_model` in `R/event_model.R` to add the `engine` argument and the internal dispatch logic.
*   **[TE-12]** Add comprehensive unit tests in `tests/testthat/test_engine_toeplitz.R`, comparing numerical outputs against the classic engine.
*   **[TE-13]** Add documentation for the new `engine` argument in `event_model`.
*   **[TE-14]** (Optional) Add benchmarking tests.
*   **[TE-15]** (Optional) Explore FFT-based alternative within `toeplitz_engine` for potential further speedup, perhaps activated by a sub-option `engine = "toeplitz-fft"`. 