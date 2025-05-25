# Event Model Refactoring To-Do List (event_model_todo.md)

**Strategy:** Refactor `event_model` construction into a three-stage pipeline (`parse` -> `realise terms` -> `build design matrix`) to unify entry points (`.formula`, `.list`, `create_event_model`), eliminate redundant computations, improve performance, and simplify the codebase. Leverage the refactored `event` and `event_term` structures. Update contrast handling to use `rlang` and align with the new pipeline.

---

**Phase 1: Implement Core Pipeline Components**

*   [X] **Ticket EM-1: Implement `parse_event_model` Helper**
    *   **Scope:** Create the first stage of the pipeline: parsing inputs.
    *   **Deliverable:** Internal function `parse_event_model(formula_or_list, data, block, durations)` in `R/event_model_helpers.R` (marked `@noRd`).
        *   Accepts either a formula or a list of `hrfspec` objects.
        *   If formula, calls `parse_event_formula` (Ticket EM-2).
        *   If list, validates `hrfspec` objects, converts to spec table, checks for `onset` column in `data` (with **clear error message** if missing).
        *   Resolves `onsets` (checking NAs).
        *   Resolves `blockids` (handling formula/vector input, checking NAs, using base `stopifnot()` for non-decreasing check, **warning if non-strictly increasing**).
        *   **Canonicalizes `blockids`** to consecutive integers 1..R using `match(..., unique(...))`. 
        *   Resolves `durations` (recycling using `recycle_or_error`, checking NAs).
        *   Returns a list containing: `spec_tbl`, `onsets`, `durations`, canonicalized `blockids`, `data`, `formula_env`.
    *   **Done when:** `parse_event_model` correctly normalizes inputs and extracts/validates/canonicalizes shared event timing information.

*   [X] **Ticket EM-2: Implement `parse_event_formula` Helper**
    *   **Scope:** Create the formula parsing logic needed by `parse_event_model`.
    *   **Deliverable:** Internal function `parse_event_formula(formula, data)` in `R/event_model_helpers.R` (marked `@noRd`).
        *   Takes the model formula and data.
        *   Uses `rlang` to identify terms on RHS.
        *   Includes helper `find_and_eval_hrf_calls` (marked `@noRd`) that recursively finds `hrf()` etc. calls and **gracefully ignores non-target calls/symbols** (returns `list()`).
        *   Evaluates `hrf()`/`trialwise()` calls *once* using **strict evaluation masking** (e.g., `rlang::eval_tidy(..., data=mask, env=f_env)` where mask uses `env_bury`) in the formula's environment.
        *   Handles potential expressions within `hrf()` calls.
        *   Extracts the LHS variable name and evaluates it for onsets.
        *   Returns `spec_tbl`, `onsets` vector, and `formula_env`.
    *   **Done when:** Function correctly parses formula into `hrfspec` objects, extracts onsets, and returns environment.

*   [X] **Ticket EM-3: Implement `realise_event_terms` Helper**
    *   **Scope:** Create the second stage: turning specs into `event_term` objects.
    *   **Deliverable:** Internal function `realise_event_terms(parsed_spec, sampling_frame, drop_empty, progress = FALSE)` in `R/event_model_helpers.R` (marked `@noRd`).
        *   Takes output from `parse_event_model`.
        *   Iterates through `spec_tbl` (e.g., `lapply` or `purrr::map`).
        *   **Optionally displays progress** using `cli::cli_progress_bar`.
        *   Calls `construct_event_term` (Ticket EM-4) for each spec.
        *   Returns a named list of `event_term` objects.
    *   **Done when:** Function correctly generates a list of `event_term` objects, with optional progress display.

*   [X] **Ticket EM-4: Implement/Refactor `construct_event_term` Helper**
    *   **Scope:** Ensure a helper exists to create an `event_term` from a single `hrfspec` and associated model info.
    *   **Deliverable:** Internal function `construct_event_term(hrfspec, model_spec)` in `R/event_model_helpers.R` (marked `@noRd`).
        *   Takes `hrfspec` and the output list from `parse_event_model` (`model_spec`).
        *   Evaluates variables/expressions from `hrfspec` using **strict masking** (using `model_spec$data` and `model_spec$formula_env`).
        *   **Handles potential `hrfspec`-specific onsets** if the spec allows overriding global onsets (check `hrfspec` definition).
        *   Adds check for factor terms resulting in zero events after subsetting; raise informative error (moved to EM-19).
        *   Calls public `event_term()` with evaluated variables and shared timing info from `model_spec`.
        *   Returns the resulting `event_term` object.
    *   **Done when:** Function correctly creates a single `event_term` with robust evaluation.

*   [X] **Ticket EM-5: Implement `build_design_matrix.event_model` Helper**
    *   **Scope:** Create the third stage: building the full design matrix efficiently.
    *   **Deliverable:** Internal function `build_event_model_design_matrix(terms, sampling_frame, precision, parallel = FALSE)` in `R/event_model_helpers.R` (marked `@noRd`).
        *   Takes list of realised `event_term` objects.
        *   Calls `convolve.event_term` for each term. **Consider `future.apply::future_lapply` if `parallel=TRUE`.**
        *   `cbind`s resulting matrices.
        *   Calculates cumulative column counts and attaches as `term_spans` attribute. Also attach `col_indices` (named list: term name -> numeric index vector).
        *   Returns final design matrix `tibble`.
    *   **Done when:** Function correctly combines term matrices, attaches attributes, with optional parallel hook.

*   [ ] **Ticket EM-18: Implement HRF Basis Cache/Memoization (Optional)**
    *   **(Deferred)**

---

**Phase 2: Refactor Public API and Main S3 Methods**

*   [ ] **Ticket EM-6: Refactor `event_model()` Public Constructor**
    *   **Scope:** Replace old entry points with calls to the new pipeline.
    *   **Deliverable:** Update `event_model()` S3 generic/methods.
        *   Calls `parse_event_model`, `realise_event_terms`, `build_event_model_design_matrix`.
        *   Constructs final `event_model` object (with `$terms`, `$design_matrix` containing `term_spans` & `col_indices` attrs, etc.).
        *   Assigns class `c("event_model", "list")`.
        *   **Deprecates** old `event_model.list`, `event_model.formula`, `create_event_model` using `.Deprecated()` wrappers.
    *   **Done when:** `event_model()` uses new pipeline, returns expected structure, old entries deprecated.

*   [X] **Ticket EM-7: Implement Vectorized `.sanitizeName` Helper**
    *   **Scope:** Provide a robust name sanitization helper.
    *   **Deliverable:** Add vectorized `.sanitizeName` to `R/utils-internal.R` (e.g., `make.names(..., unique=TRUE)`). **Store original names as attribute** for lossless round-tripping.
    *   **Done when:** Vectorized sanitization helper with attribute is available.

*   [X] **Ticket EM-8: Replace `design_matrix.event_term`**
    *   **Scope:** Replace the current `design_matrix.event_term` with the more direct 25-line version proposed in the feedback, removing the `locenv` detour.
    *   **Deliverable:** Updated `design_matrix.event_term` function in `R/event_vector.R` (or `R/event-term.R` if created).
    *   **Done when:** The method uses the direct data frame approach.

*   [X] **Ticket EM-9: Update `hrf()` for Contrasts**
    *   **Scope:** Allow only symbolic contrasts in `hrf()`.
    *   **Deliverable:** Modified `hrf()` in `R/hrf.R` with updated validation. **Emit clear error with hint if numeric contrast is passed.**
    *   **Done when:** `hrf()` only accepts symbolic contrasts with clear error message.

---

**Phase 3: Update Contrast System**

*   [X] **Ticket EM-10: Update Contrast Weights Functions to use `rlang`**
    *   **Scope:** Modernize contrast weight calculation.
    *   **Deliverable:** Modify `contrast_weights.*` in `R/contrast.R`:
        *   Replace `lazyeval` calls with `rlang` equivalents (`is_formula`, `f_rhs`, `eval_tidy`).
        *   **Optional:** Memoize `where` clause evaluation.
    *   **Done when:** `lazyeval` dependency removed.

*   [X] **Ticket EM-11: Update `makeWeights` Helper**
    *   **Scope:** Make internal weight calculation robust.
    *   **Deliverable:** Replace `makeWeights` in `R/contrast.R` with `.make_weights` using logical masks and returning **named numeric vector**. Update callsites.
    *   **Done when:** Robust `.make_weights` implemented and used.

*   [X] **Ticket EM-12: Clean up `Fcontrasts.event_term`**
    *   **Scope:** Remove unnecessary complexity due to guaranteed full-rank design matrix from `design_matrix.event_term`.
    *   **Deliverable:** Modify `Fcontrasts.event_term` to **remove SVD/QR fallback logic**.
    *   **Done when:** SVD/QR fallback removed.

---

**Phase 4: Cleanup and Verification**

*   [ ] **Ticket EM-13: Update Tests for `event_model`**
    *   **Scope:** Ensure unit tests cover the new `event_model` pipeline and verify backward compatibility.
    *   **Deliverable:** Update tests in `tests/testthat/test_event_model.R`:
        *   Test `event_model()` implicitly covering pipeline stages (formula/list inputs, multiple runs, HRFs, interactions, modulators, **expressions**, **contrasts**).
        *   Add/update round-trip sanity checks (**incl. cat × cont interaction**).
        *   Verify `event_model` object structure (**`term_spans`, `col_indices` attrs**).
        *   Verify `contrast_weights.event_model` works with new attributes.
        *   Use snapshots for `design_matrix(event_model)` comparison.
    *   **Done when:** Tests cover new structure and confirm compatibility.

*   [ ] **Ticket EM-14: Remove Old Code / Deprecate**
    *   **Scope:** Delete/deprecate functions/helpers made redundant.
    *   **Deliverable:**
        *   **Deprecate** `event_model.list`, `.formula`, `create_event_model` (handled in EM-6).
        *   Remove internal helpers from `R/event_model.R`: `construct_model` (old), `.pre_eval_hrf_calls`, `extract_terms`, `extract_variables`, `parse_term`, `.extract_expr_variables_into_events`.
    *   **Done when:** Redundant code removed/deprecated.

*   [ ] **Ticket EM-15: Update Documentation**
    *   **Scope:** Document the new structure and update existing documentation.
    *   **Deliverable:**
        *   Update Roxygen for `event_model()`.
        *   Add docs for new internal helpers (as needed).
        *   Add deprecation notes for old constructors.
        *   Update vignettes (`a_04_event_models.Rmd`, `Overview.Rmd`, `a_05_contrasts.Rmd`) re: API changes (symbolic contrasts only) and **add note about cross-term contrasts being future work**.
        *   Document raw numeric contrasts are unsupported.
    *   **Done when:** Documentation is accurate and reflects refactor/deprecations.

*   [ ] **Ticket EM-16: Final Namespace/Import Check**
    *   **Scope:** Verify imports and exports.
    *   **Deliverable:** Ensure necessary imports (e.g., `rlang`, `tibble`, `purrr`, `cli`, **`future.apply` (optional)**) present. Remove unused (`lazyeval`). Check exports/internal markers. Check `R CMD check`.
    *   **Done when:** Namespace is clean.

*   [ ] **Ticket EM-17: Add Performance Benchmark Script**
    *   **Scope:** Add a script to track performance improvements.
    *   **Deliverable:** Create script in `inst/benchmarks/` (or similar) comparing old vs. new `event_model` design matrix build time on a representative synthetic dataset. Potentially add assertion for CI.
    *   **Done when:** Benchmark script exists.

*   [ ] **Ticket EM-19: Add Error Check for Empty Factors**
    *   **(Moved from EM-4)**
    *   **Scope:** Improve error handling.
    *   **Deliverable:** Ensure `construct_event_term` (or earlier) raises informative error if a factor has 0 events after subsetting.
    *   **Done when:** Informative error is raised.

*   [ ] **Ticket EM-20: Implement Cross-Term Contrasts (Post-Refactor)**
    *   **Scope:** Allow users to define contrasts spanning multiple terms within an `event_model`.
    *   **Deliverable:**
        *   Implement a new S3 method, likely `contrast_weights.contrast_formula_spec_event_model` (or potentially adapt `contrast_weights.contrast_formula_spec` to detect context).
        *   This method should handle contrast formulas referencing conditions within specific terms using a `TermName$ConditionName` syntax (e.g., `~ Term1$Face - Term2$Scene`).
        *   It should use the `term_spans` and/or `col_indices` attributes stored on the event model's design matrix to identify the correct columns corresponding to `TermName$ConditionName`.
        *   It will evaluate the contrast formula in an environment where `TermName$ConditionName` is mapped to a logical mask or index vector for the relevant columns in the full design matrix.
        *   Leverage the `.make_weights` helper.
        *   Add corresponding unit tests for cross-term contrast definition and calculation.
        *   Update documentation and vignettes to explain and demonstrate cross-term contrasts.
    *   **Done when:** Users can successfully define and compute weights for contrasts spanning multiple event terms using the symbolic formula interface. 