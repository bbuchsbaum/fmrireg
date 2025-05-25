# Event Vector Refactoring To-Do List

**Strategy Update:** Based on the "Modernisation Road-Map" feedback and subsequent discussion, this refactor focuses on unifying event representation around a single internal constructor (`event()`) and a consistent internal structure (**matrix-based value storage**), simplifying S3 methods while maintaining backward compatibility through **thin wrappers** for existing `event_*` functions. Higher-level components (`event_term`, `design_matrix`, contrasts) are expected to require minimal-to-no changes due to preserved S3 method output signatures.

---

**Phase 1: Unify Event Object Creation & Core Methods**

*   [X] **Ticket EV-1: Implement Internal `event()` Constructor**
    *   **Scope:** Create a single, internal constructor for all event types.
    *   **Deliverable:** `R/event-classes.R` containing:
        *   Internal `event(...)` function calling `.checkEVArgs` first, then applying subsetting.
        *   Stores core data: `$varname`, `$onsets`, `$blockids`, `$durations`, `$continuous` (logical flag).
        *   Stores payload in `$value` (**always** a matrix, numeric type).
            *   **Factors/Characters:** Stored as integer codes; actual levels stored in `$meta$levels`.
            *   **Basis Objects:** Numeric matrix (`basis$y`) stored; full subsetted basis object stored in `$meta$basis`.
            *   **Numeric/Matrix:** Stored as is (subsetted).
        *   Sets class `c("event", "event_seq")`.
    *   **Done when:** `event()` constructor implemented with the specified internal storage and attribute conventions.

*   [X] **Ticket EV-2: Refactor `event_*` Functions as Thin Wrappers**
    *   **Scope:** Turn existing `event_factor`, `event_variable`, `event_matrix`, `event_basis` into thin wrappers around the new `event()` constructor.
    *   **Deliverable:** Modify existing functions in `R/event-classes.R` (or keep in `R/event_vector.R` initially) to:
        *   Call the new internal `event()` constructor, passing arguments appropriately (e.g., `event_factor` sets `durations=0` and ensures value is character).
        *   **No `.Deprecated()` calls initially** to maintain backward compatibility without warnings. These wrappers become the primary user-facing API for now.
    *   **Done when:** Old functions are simple wrappers calling the new `event()` constructor.

*   [X] **Ticket EV-3: Unify Core S3 Methods**
    *   **Scope:** Consolidate essential S3 methods to work on the base `event` class using the new internal structure and attributes.
    *   **Deliverable:**
        *   Implement simple, generic S3 methods in `R/event-classes.R` (or similar):
            *   `levels.event(x)`: Retrieve names. Check `x$meta$basis` -> `levels(basis_obj)`; Check `x$meta$levels` -> factor levels; Default to `colnames(x$value)`.
            *   `columns.event(x)` -> `levels.event(x)` (Alias).
            *   `is_continuous.event(x)` -> `x$continuous`.
            *   `is_categorical.event(x)` -> `!x$continuous`.
            *   `cells.event(x, ...)`: Generate tibble. For factors (if `!is_continuous(x)` and `!is.null(x$meta$levels)`), reconstruct character values from integer codes and `levels` attribute. Ensure `attr(., "count")` is preserved.
        *   Implement `elements.event(x)`.
        *   Remove redundant specific S3 methods for old `event_*` classes.
    *   **Done when:** Generic methods work correctly, handling attributes; redundant methods removed.

---

**Phase 2: Slim Helpers & Utilities**

*   [X] **Ticket EV-4: Review & Replace Helper Functions**
    *   **Scope:** Replace older/custom helper functions used within the event system with base R or standard alternatives where possible.
    *   **Deliverable:**
        *   Review functions like `is.increasing`/`is.strictly.increasing` (replaced with `!is.unsorted(...)`), `.sanitizeName` (kept, moved to utils), `conform_len` (not used).
        *   Refactor/simplify `.checkEVArgs` used by `event()` (kept, moved to utils, minor improvements made).
    *   **Done when:** Legacy helpers used *by the event system* are replaced or confirmed necessary; code simplified.

*   [X] **Ticket EV-5: Centralize Subsetting Logic**
    *   **Scope:** Ensure subsetting logic (via the `subset` argument) is applied consistently and only once within the `event()` constructor path.
    *   **Deliverable:** Modify `event()` constructor (and/or `.checkEVArgs`) to handle the `subset` argument definitively. Remove any lingering subsetting logic from wrapper functions.
    *   **Done when:** Subsetting occurs only within the core constructor path (Confirmed done in EV-1/EV-2/EV-4).

*   [X] **Ticket EV-20: Warn on NA Event Values**
    *   **Scope:** Add robustness to the `event()` constructor.
    *   **Deliverable:** Modify `event()` in `R/event-classes.R` to issue a `warning()` if the input `value` (before processing) contains any `NA` values.
    *   **Done when:** Warning is implemented.
    
*   [X] **Ticket EV-21: Optimize `elements.event(values=FALSE)`**
    *   **Scope:** Improve efficiency and reduce memory usage for `elements.event` when `values=FALSE`.
    *   **Deliverable:** Modify `elements.event` in `R/event-classes.R` so that for continuous/matrix/basis types when `values=FALSE`, it returns a single row matrix/vector of the column names/levels instead of repeating them for every event.
    *   **Done when:** `elements.event` returns non-repeated names/levels for `values=FALSE`.

---

**Phase 3: Verify/Refactor `event_term` and `design_matrix`**

*   [X] **Ticket EV-6: Verify `event_term` Compatibility**
    *   **Scope:** Ensure `event_term` works correctly with the new unified `event` objects created via the wrappers.
    *   **Deliverable:** Verify that the `EV` factory function within `event_term` calls the `event_*` wrappers, which return valid `event` objects, and that subsequent operations within `event_term` (which likely use `levels`, `cells`, etc.) function correctly due to the preserved output signatures of the S3 methods (EV-3). Minimal changes to `event_term` itself are expected.
    *   **Done when:** `event_term` correctly produces terms containing unified `event` objects and downstream methods work (Confirmed compatible).

*   [X] **Ticket EV-7: ~~Simplify~~ Verify `design_matrix.event_term`**
    *   **Scope:** Verify that `design_matrix.event_term` produces the correct output.
    *   **Deliverable:**
        *   Confirm that `design_matrix.event_term` works correctly. **Crucially, verify it checks `is_categorical(ev)` / `!is_continuous(ev)` and uses `attr(ev$value, "levels")` when calling `model.matrix` for factor-like events now stored as integers.** (Confirmed compatible - uses `event_table` built from `elements(..., values=FALSE)` which handles reconstruction).
        *   (Optional simplification remains optional).
    *   **Done when:** `design_matrix.event_term` confirmed to work correctly with integer-coded factors.

*   [ ] **Ticket EV-8: Consolidate `design_matrix` Methods (If Applicable) - DEFERRED**
    *   **Note:** Deferring this for now. Current focus is on core event object refactor.
    *   **Scope:** If other `design_matrix.*` methods exist for specific old event-related classes, consolidate them into `design_matrix.fmri_term` or remove if redundant after EV-3/EV-7 verification.
    *   **Deliverable:** Reduced number of `design_matrix` S3 methods if applicable.
    *   **Done when:** Code duplication across methods is minimized.

---

**Phase 4: Verify Contrast Generation**

*   [X] **Ticket EV-9: ~~Implement `build_contrast` Helper~~ Verify Contrast Logic**
    *   **Scope:** Verify that existing contrast logic works unchanged.
    *   **Deliverable:** Confirm that functions relying on contrast generation (e.g., `Fcontrasts`, `contrast_weights`) work correctly because they use methods like `cells()` which preserve the required `attr(., "count")` (EV-3). (Confirmed compatible).
    *   **Done when:** Existing contrast generation confirmed compatible.

*   [X] **Ticket EV-10: ~~Refactor `Fcontrasts.event_term`~~ Verify `Fcontrasts.event_term`**
    *   **Scope:** Verify that `Fcontrasts.event_term` works correctly without modification.
    *   **Deliverable:** Confirm function produces identical contrast matrices due to preserved input expectations (from `cells`, `levels`, etc.). (Confirmed compatible).
    *   **Done when:** Function is confirmed compatible.

---

**Phase 5: Plotting & Printing (No Change from Original Plan)**

*   [ ] **Ticket EV-11: Implement `autoplot` Methods**
    *   **Scope:** Add ggplot2-based plotting for the unified `event` objects and `event_term`.
    *   **Deliverable:** `autoplot.event`, `autoplot.event_term` methods.
    *   **Done when:** Visualizations can be generated using `autoplot`.

*   [X] **Ticket EV-12: Deprecate Base `plot` Methods**
    *   **Scope:** Add deprecation warnings to old plot methods if they exist for specific `event_*` classes.
    *   **Deliverable:** `.Deprecated()` calls added to `plot.event_*` methods, pointing to `autoplot`.
    *   **Done when:** Plotting warns and points to `autoplot`.

*   [X] **Ticket EV-13: Consolidate `print` Methods**
    *   **Scope:** Create unified print methods using `cli`.
    *   **Deliverable:** Implement `print.event` and `print.event_term` using `cli`. Deprecate/remove specific `print.event_*` methods.
    *   **Done when:** Printing uses new methods; redundancy removed.

---

**Phase 6: Documentation & Tests (Adjusted)**

*   [X] **Ticket EV-14: Update Roxygen Docs & Examples**
    *   **Scope:** Document the refactored system, focusing on the unified `event` object's internal structure and the public-facing wrappers.
    *   **Deliverable:** Updated Roxygen for `event()` (internal), `event_*` wrappers, unified S3 methods (`levels.event`, `cells.event`, etc.), `event_term`. Examples should continue using the `event_*` wrappers primarily. **Add note to `cells.event_term` docs re: potential slowness.**
    *   **Done when:** Documentation is accurate, reflects the internal structure, and builds correctly.

*   [ ] **Ticket EV-15: Update Vignettes**
    *   **Scope:** Revise vignettes using event objects (e.g., Overview, Event Models) to ensure accuracy, even if user-facing code examples don't change significantly.
    *   **Deliverable:** Vignettes updated and build correctly.
    *   **Done when:** Vignettes reflect any subtle changes and are accurate.

*   [X] **Ticket EV-16: Add/Update Unit Tests**
    *   **Scope:** Ensure test coverage for the refactored event system.
    *   **Deliverable:** Update `tests/testthat/test-event_vector.R`. Tests must cover:
        *   `event()` constructor: correct handling of inputs, **correct integer codes and attributes (`levels`, `basis`) in output `$meta`**.
        *   Unified S3 methods (EV-3): **correct output across different event types, especially correct level/column names and reconstructed cell values using attributes.**
        *   Downstream compatibility (EV-6, EV-7, EV-9, EV-10): **Verify identical results for `event_term`, `design_matrix`, contrasts compared to pre-refactor code, confirming handling of new internal format (via snapshots).**
        *   **Add specific "round-trip sanity" test covering wrappers -> event_term -> design_matrix/levels/cells.**
    *   **Done when:** Tests pass with high coverage, confirming internal logic and backward compatibility.

---

**Phase 7: Cleanup & File Organization (Adjusted)**

*   [X] **Ticket EV-17: Remove Orphaned Code**
    *   **Scope:** Delete unused helper functions and **redundant S3 methods specific to old `event_*` classes** (e.g., `levels.event_factor`, `cells.event_matrix`).
    *   **Deliverable:** Removal of identified dead code.
    *   **Done when:** `R CMD check` passes cleanly.

*   [X] **Ticket EV-18: Organize Files**
    *   **Scope:** Structure code into logical files.
    *   **Deliverable:** Code organized into files like `R/event-classes.R` (containing `event()` constructor, wrappers, generic S3 methods), `R/utils-internal.R` (helpers). `event_term` related code remains in `R/event_vector.R` for now (or dedicated file if refactored later).
    *   **Done when:** Files reorganized logically.

*   [X] **Ticket EV-19: Final Namespace/Import Check**
    *   **Scope:** Verify all necessary imports and namespacing.
    *   **Deliverable:** `R CMD check` passes namespace checks.
    *   **Done when:** Namespace is clean. 