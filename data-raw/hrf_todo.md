# HRF Refactoring To-Do List

This plan outlines the incremental steps to refactor the `R/hrf.R` module based on the principles in `hrf_refactor.md`. The goal is a more modular system using `as_hrf`, decorators, and a single `evaluate.HRF`, verified against `tests/testthat/test_hrf.R`.

---

**Phase 1: Core Constructor and Basis Binding**

*   [x] **Ticket 1: Implement `as_hrf` Constructor**
    *   Create the `as_hrf(f, name, nbasis, span, params)` function in a new file (e.g., `R/hrf-core.R` or within `R/hrf.R`).
    *   Ensure it correctly structures the input function `f` with the specified attributes (`class`, `name`, `nbasis`, `span`, `param_names`, `params`).
    *   Add basic unit tests for `as_hrf`.

*   [x] **Ticket 2: Implement `bind_basis` Decorator**
    *   Implement the `bind_basis(...)` function.
    *   It should take multiple HRF objects and combine them into a single multi-basis HRF object using `as_hrf`.
    *   It must correctly calculate the combined `nbasis` and maximum `span`.
    *   The resulting function should `cbind` the outputs of the input HRFs.
    *   Add unit tests for `bind_basis`.

*   [x] **Ticket 3: Refactor Canonical HRF Objects**
    *   Rewrite `HRF_GAMMA`, `HRF_SPMG1`, `HRF_GAUSSIAN` using `as_hrf`.
    *   Rewrite `HRF_SPMG2` and `HRF_SPMG3` using `as_hrf`, the base `hrf_spmg1*` functions, and `bind_basis`.
    *   Rewrite `HRF_BSPLINE` using `as_hrf` and `hrf_bspline`. Consider how `nbasis` and `span` parameters are handled.
    *   **Verification:** Run `tests/testthat/test_hrf.R`. All existing tests related to these objects should still pass, verifying structure and numerical output.

**Phase 2: Implement Decorators**

*   [x] **Ticket 4: Implement `lag_hrf` Decorator**
    *   Implement `lag_hrf(hrf, lag)`.
    *   It should return a new HRF object (via `as_hrf`) whose function evaluates the original `hrf` at `t - lag`.
    *   Update the `name` and `span` attributes correctly.
    *   Add unit tests for `lag_hrf`, potentially adapting existing tests from `gen_hrf_lagged`.

*   [x] **Ticket 5: Implement `block_hrf` Decorator**
    *   Implement `block_hrf(hrf, width, precision, half_life, summate, normalise)`.
    *   It should return a new HRF object (via `as_hrf`) whose function implements the block convolution logic (sampling offsets, applying decay, summing or taking max).
    *   Update the `name` and `span` attributes correctly.
    *   Add unit tests for `block_hrf`, potentially adapting existing tests from `convolve_block` or `gen_hrf_blocked`.

*   [x] **Ticket 6: Implement `normalise_hrf` Decorator**
    *   Implement `normalise_hrf(hrf)`.
    *   It should return a new HRF object (via `as_hrf`) whose function evaluates the original `hrf` and scales the output so the maximum absolute value is 1.
    *   Update the `name` attribute.
    *   Add unit tests for `normalise_hrf`.

**Phase 3: Integrate Decorators and Refactor `gen_hrf`**

*   [x] **Ticket 7: Refactor `gen_hrf`**
    *   Modify the existing `gen_hrf` function.
    *   Instead of its custom logic for `lag` and `width`, it should now call the `lag_hrf` and `block_hrf` decorators conditionally based on the arguments.
    *   The core `hrf` function passed to `gen_hrf` should first be converted using `as_hrf` if it's a raw function.
    *   Handle the `normalize` argument potentially by calling `normalise_hrf`.
    *   **Verification:** Ensure tests involving `gen_hrf` in `tests/testthat/test_hrf.R` still pass.

**Phase 4: Unify Evaluation and Remove Redundancy**

*   [x] **Ticket 8: Implement Unified `evaluate.HRF`**
    *   Replace the existing `evaluate.HRF` method body with the logic proposed in `hrf_refactor.md`.
    *   This method should handle both point evaluation (`duration = 0`) and block evaluation (`duration > 0`), incorporating `amplitude`, `precision`, `summate`, and `normalize` arguments.
    *   **Verification:** Rerun all tests in `tests/testthat/test_hrf.R` involving `evaluate()`.
        Ensure results are consistent with previous versions.

*   [x] **Ticket 9: Remove Redundant Block/Convolution Functions**
    *   Delete the `convolve_block` function.
    *   Delete the `gen_hrf_blocked` function (and its alias `hrf_blocked`).
    *   Search the codebase for any internal uses and replace them with calls to `block_hrf` or direct use of `evaluate.HRF` with a duration.
    *   **Verification:** Rerun all tests.

*   [ ] **Ticket 10: Remove Redundant Lag/Set Generators**
    *   Delete `gen_hrf_lagged` (and alias `hrf_lagged`). Its functionality is covered by `lag_hrf`.
    *   Delete `gen_hrf_set`. Its functionality is covered by `bind_basis`.
    *   Search the codebase for any internal uses and replace them.
    *   **Verification:** Rerun all tests.

**Phase 5: Final Integration and Cleanup**

*   [x] **Ticket 11: Review `hrfspec` and `construct` Methods**
    *   Examine `hrf()`, `hrfspec()`, `construct.hrfspec()`, and related methods (`trialwise`, `afni_hrf`, etc.).
    *   Ensure they correctly handle the creation and use of HRF objects built with `as_hrf` and decorators.
    *   Simplify the logic in `hrf()` (e.g., `make_hrf` might be simplified or removed if `gen_hrf` handles the creation properly).
    *   Ensure the `basis` argument processing is robust.

*   [ ] **Ticket 12: Review Remaining `R/hrf.R` Functions**
    *   Review AFNI functions, `hrf_time`, `hrf_ident`, `hrf_sine`, `hrf_inv_logit`, `hrf_half_cosine`, `hrf_toeplitz`, Daguerre functions, print methods, `list_available_hrfs`.
    *   Refactor `getHRF`: Implement the registry (`.HRF_LIBRARY`) and factory pattern (`hrf_factory`) approach to replace the `switch` statement. Keep the external function name as `getHRF` but use the new internal logic. Ensure aliases (`gam`, `bs`) are handled in the registry.
    *   Ensure print methods work with the new structure.
    *   Update `list_available_hrfs` if necessary.
    *   Refactor any remaining functions to use the new core components where applicable.

*   [ ] **Ticket 13: Final Cleanup and Documentation**
    *   Remove any unused internal helper functions (e.g., potentially parts of `make_hrf`).
    *   Add/update Roxygen documentation (docstrings) for all new and refactored functions (`as_hrf`, decorators, `evaluate.HRF`).
    *   Update examples in the documentation to demonstrate the new composition pipeline (`hrf |> lag_hrf(2) |> block_hrf(5)`).
    *   Ensure consistent coding style.
    *   **Verification:** Final run of all tests in `tests/testthat/test_hrf.R`. 

*   [x] **Ticket 13: Refactor trialwise() as a Thin Wrapper**
    *   [x] Update `trialwise()` in `R/hrf.R` to use `rlang::quo(.trial_factor(length(onsets)))` and pass `id=label`.
    *   [x] Remove `width`/`normalize` args from `trialwise()`. 
    *   [x] Add `$add_sum` flag to the `hrfspec` returned by `trialwise()`.

*   [x] **Ticket 14: Remove Bespoke trialwisespec Code**
    *   [x] Delete `trialwisespec` class and `construct.trialwisespec` method from `R/hrf.R`.
    *   [x] Update `construct.hrfspec` in `R/hrf-formula.R` to handle the `$add_sum` flag, calculating and appending the mean column.

*   [x] **Ticket 15: Add Helper for Trial ID Factor**
    *   [x] Implement internal `.trial_factor()` helper in `R/hrf.R`.

*   [ ] **Ticket 16: Verify and Test Refactored trialwise()**
    *   Add unit tests in `tests/testthat/test_hrf.R` to compare the new `trialwise()` output with explicit `hrf()` calls.
    *   Test the `add_sum=TRUE` functionality and column naming.
    *   Check for backward compatibility and ensure no regressions in `event_model` formulas using `trialwise`.

*   [x] **Related: Refactor hrf() in R/hrf-formula.R**
    *   [x] Remove debug messages.
    *   [x] Add `name=NULL` formal argument as alias for `id`.

---
## Folding the "odds‑and‑ends" into the one‑page HRF core

Below I show how each of the three helpers melts away once we have the canonical
constructor as_hrf() and the composable decorator bind_basis().

⸻

1 Empirical HRF → 3‑line wrapper

# old ------------------------------------------------------
gen_empirical_hrf <- function(t, y, name = "empirical_hrf") {
  f <- approxfun(t, y, yright = 0, yleft = 0)
  HRF(f, name = name, nbasis = 1)
}

# new ------------------------------------------------------
empirical_hrf <- function(t, y, name = "empirical_hrf") {
  as_hrf(approxfun(t, y, yright = 0, yleft = 0),
         name = name, nbasis = 1, span = max(t))
}

Exactly the same public behaviour, just written with the single constructor.
If we want to keep the old name for backward compatibility:

gen_empirical_hrf <- function(...) {
  .Deprecated("empirical_hrf()")
  empirical_hrf(...)
}



⸻

2 Basis‑set combiner → bind_basis()

# old ------------------------------------------------------
gen_hrf_set <- function(..., span = 32, name = "hrf_set") { … }

# new ------------------------------------------------------
hrf_set <- function(..., name = "hrf_set") {
  bind_basis(...) |>
    (\(h) { attr(h, "name") <- name; h })()
}

All the busywork (checking each element is a function, cbind‑eval, computing
nbasis and span) already lives inside bind_basis().
Again, deprecate the alias if you wish:

gen_hrf_set <- function(...) {
  .Deprecated("hrf_set()")
  hrf_set(...)
}



⸻

3 Parameter‑grid library → one tidy‑map call

# old ------------------------------------------------------
gen_hrf_library <- function(fun, pgrid, ...) { … }

# new ------------------------------------------------------
hrf_library <- function(fun, pgrid, ...) {
  purrr::pmap(pgrid, \(...) fun |> purrr::partial(...) |> as_hrf(...)) |>
    (\(fs) do.call(bind_basis, fs))()
}

Differences / improvements

Aspect	Old code	New one‑liner
Handles arbitrary parameter names?	yes	yes (pmap)
Carries each HRF's own attributes?	partly	yes (as_hrf() inside the loop)
Combines into basis set?	ad‑hoc gen_hrf_set()	bind_basis()
Lines of code	~25	5

Deprecation wrapper:

gen_hrf_library <- function(...) {
  .Deprecated("hrf_library()")
  hrf_library(...)
}



⸻

Minimal impact checklist ✔︎

Down‑stream function	Still works?	Why
evaluate()	✔︎	All objects inherit function + "HRF", attributes unchanged.
nbasis(), print.*	✔︎	Attributes set by as_hrf() / bind_basis().
Formula DSL hrf(cond, basis = …)	✔︎	make_hrf() calls as_hrf()‑based constructors.
Legacy code	✔︎ w/ message	Compatibility shims using .Deprecated().



⸻

Take‑away

Those "odds and ends" collapse to ≤ 20 tidy lines once the single‑source
primitives (as_hrf(), bind_basis()) exist.
This keeps the public API intact and enforces the one constructor / one
decorator rule we set for maximal elegance and maintainability. 