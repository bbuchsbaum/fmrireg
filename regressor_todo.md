# Regressor Refactoring To-Do List (v2 - Incorporating Feedback)

Based on the plan in `regressor_refactor.md` and subsequent detailed feedback.

---

**📦 Phase A: Core Constructor & API Setup**

*   [x] **Ticket A-1: Implement Core `Reg` Constructor & Public Facade**
    *   **Scope:** Create internal `Reg()` & make public `regressor()` call it.
    *   **Deliverable:** `R/regressor.R` (or new `R/reg-constructor.R`) containing:
        *   Internal `Reg()` constructor using `recycle_or_error` (B-2), one-pass amplitude filtering (B-3), `make_hrf` call, span logic (`%||%`), negative duration/span checks, and setting `class = c("regressor", "Reg", "list")`.
        *   Public `regressor()` function as a simple wrapper calling `Reg(...)`.
    *   **Done when:** `regressor()` creates objects with class `c("regressor", "Reg", "list")`, passes validation/filtering tests.

*   [x] **Ticket A-2: Verify Helper Constructors** (was API-1)
    *   **Scope:** Ensure `single_trial_regressor` and `null_regressor` function correctly as wrappers.
    *   **Deliverable:** Confirm `single_trial_regressor` and `null_regressor` call `Reg()` correctly.
    *   **Done when:** Wrappers function as intended without deprecation.

*   [x] **Ticket A-3: Export `evaluate` Generic (Verify)** (was A-4)
    *   **Scope:** Ensure `evaluate` S3 generic is defined and exported.
    *   **Deliverable:** Check `R/all_generic.R` (or similar) for `evaluate <- function(x, grid, ...) UseMethod("evaluate")` and `@export` tag. Add if missing.
    *   **Done when:** Generic is confirmed exported.

*   [x] **Ticket API-2: Remove `shift.regressor` Method**
    *   **Scope:** Delete the redundant S3 method for the legacy class.
    *   **Deliverable:** Remove `shift.regressor` definition.
    *   **Done when:** Code compiles; `shift()` calls dispatch correctly to `shift.Reg`.

---

**🧬 B. Internal representation & Helpers**

*   [x] **Ticket B-1: Define `Reg` S3 Class Structure Docs**
    *   **Scope:** Document the fields stored by `Reg()`.
    *   **Deliverable:** Roxygen documentation in `Reg()` outlining list components: `onsets`, `duration`, `amplitude`, `hrf`, `span`, `summate`.
    *   **Done when:** Roxygen docs accurately describe the structure.

*   [x] **Ticket B-2: Implement `recycle_or_error` & `%||%` Helpers** (was B-2)
    *   **Scope:** Create internal utility functions.
    *   **Deliverable:** `R/utils-internal.R` containing the simple base-R `recycle_or_error` function and the `%||%` operator.
    *   **Done when:** Helpers exist. `Reg()` constructor updated to use `recycle_or_error` (as part of A-1).

*   [x] **Ticket B-3: Verify One-Pass Filtering in `Reg`** (was B-3)
    *   **Scope:** Ensure `Reg()` uses the `keep <- ...` logic correctly.
    *   **Deliverable:** Confirmation that the `Reg()` implementation (from A-1) correctly filters based on non-zero/non-NA amplitude *after* recycling.
    *   **Done when:** Filtering logic confirmed correct and tested (F-3).

*   [x] **Ticket B-4: Simplify `shift.Reg` Method**
    *   **Scope:** Update `shift.Reg` based on the simpler `Reg()` constructor.
    *   **Deliverable:** S3 method `shift.Reg(x, shift_amount, ...)` calls `Reg()` directly for reconstruction.
    *   **Done when:** `shift.Reg` works and returns object with class `c("regressor", "Reg", "list")`. Unit tests pass.

*   [x] **Ticket B-3.1: Re-validate/Recycle Duration After Filtering** (New from feedback)
    *   **Scope:** Ensure duration vector length matches filtered onsets.
    *   **Deliverable:** Add call to `recycle_or_error(duration, length(onsets), "duration")` *after* the amplitude filtering logic inside `Reg()`.
    *   **Done when:** Test cases with mismatched initial/filtered lengths are handled correctly.

---

**⚙️ C. Evaluation engine Structure**

*   [x] **Ticket C-1: Implement `evaluate.Reg` Dispatcher & Remove Old Methods**
    *   **Scope:** Create the single `evaluate.Reg` dispatcher and remove evaluate methods for `regressor`, `single_trial_regressor`, `null_regressor`.
    *   **Deliverable:**
        *   `R/reg-evaluate.R` (or similar) containing `evaluate.Reg(x, grid, precision, method, sparse, ...)` that calls `prep_reg_inputs` and dispatches to `eval_*` functions via `switch()`.
        *   Remove `evaluate.regressor`, `evaluate.single_trial_regressor`, `evaluate.null_regressor` method definitions.
    *   **Done when:** `evaluate()` correctly dispatches to `evaluate.Reg`. Old methods removed.

*   [x] **Ticket C-2: Implement `prep_reg_inputs` Helper** (was E-2 part 1)
    *   **Scope:** Create the evaluation pre-processing helper.
    *   **Deliverable:** `R/evaluate-helpers.R` (or `utils-internal.R`) containing `prep_reg_inputs(x, grid, precision)` function handling event filtering (by grid+span) and HRF memoization/sampling (using D-1 helper).
    *   **Done when:** Helper function implemented and used by `evaluate.Reg`.

*   [x] **Ticket C-3: Implement Evaluation Engines (`eval_*`)** (was E-2 part 2)
    *   **Scope:** Create internal functions for each evaluation method.
    *   **Deliverable:** `R/evaluate-helpers.R` containing `eval_fft`, `eval_conv`, `eval_Rconv`, `eval_loop`. Each takes the list from `prep_reg_inputs`.
        *   `eval_fft`, `eval_conv` call the C++ wrapper (C-4).
        *   `eval_Rconv` implements logic using `stats::convolve` (D-2).
        *   `eval_loop` implements the R loop logic.
    *   **Done when:** Engine functions are defined and callable by `evaluate.Reg` dispatcher.

*   [x] **Ticket C-4: Integrate C++ Wrapper** (was C-2)
    *   **Scope:** Ensure C++ code is correctly called.
    *   **Deliverable:** `eval_fft` and `eval_conv` correctly call `evaluate_regressor_cpp`.
    *   **Done when:** C++ code compiles; FFT/Conv methods produce correct results via internal engines.

*   [x] **Ticket C-5: Verify HRF Integration** (was C-4)
    *   **Scope:** Ensure `evaluate.HRF` is used correctly for HRF sampling.
    *   **Deliverable:** Confirmation that `prep_reg_inputs` (using `.memo_hrf`) correctly calls `evaluate.HRF`.
    *   **Done when:** Tests confirm correct HRF evaluation.

---

**🚀 D. Performance & memory Enhancements**

*   [x] **Ticket D-1: Implement Memoized HRF Helper (`.memo_hrf`)** (was D-1)
    *   **Scope:** Create the memoized HRF sampling function.
    *   **Deliverable:** Define `.memo_hrf <- memoise::memoise(...)` in `R/evaluate-helpers.R` (or `utils-internal.R`).
    *   **Done when:** Helper function defined. `prep_reg_inputs` updated to use it.

*   [x] **Ticket D-2: Implement `eval_Rconv` using `stats::convolve`** (was D-2)
    *   **Scope:** Remove `pracma` dependency.
    *   **Deliverable:** `eval_Rconv` engine uses `stats::convolve(..., type="open")`.
    *   **Done when:** "Rconv" method works without `pracma`.

*   [x] **Ticket D-3: Add Optional Sparse Output**
    *   **Scope:** Allow sparse matrix return from `evaluate.Reg`.
    *   **Deliverable:** `evaluate.Reg` has `sparse = FALSE` argument and converts final output using `Matrix::Matrix(..., sparse=TRUE)`.
    *   **Done when:** `evaluate(..., sparse=TRUE)` returns a sparse matrix.

---

**🖼️ E. Visual helpers & Print**

*   [x] **Ticket E-1: Implement `autoplot.Reg` & Deprecate `plot.regressor`**
    *   **Scope:** Create ggplot-based plot method.
    *   **Deliverable:** `autoplot.Reg` function (e.g., in `R/autoplot-methods.R`). Add `.Deprecated("autoplot")` to `plot.regressor`.
    *   **Done when:** `autoplot()` works; `plot()` warns.

*   [x] **Ticket E-2: Implement `print.Reg` & Deprecate `print.regressor`**
    *   **Scope:** Create cli-based print method.
    *   **Deliverable:** `print.Reg` function (e.g., in `R/print-methods.R` or `R/regressor.R`) using `cli`. Update `print.regressor` to call `print.Reg` and issue `.Deprecated("print")`.
    *   **Done when:** `print()` uses new method; old `print.regressor` warns.

---

**📚 F. Documentation & tests**

*   [x] **Ticket F-1: Add/Update Roxygen Docs**
    *   **Scope:** Document all new/refactored public/S3 functions.
    *   **Deliverable:** Updated Roxygen for `regressor`, `evaluate.Reg`, `shift.Reg`, `autoplot.Reg`, `print.Reg`. Add/update deprecation notes.
    *   **Done when:** Docs are accurate and build correctly.

*   [x] **Ticket DOC-1: Update Roxygen Examples**
    *   **Scope:** Ensure examples use only `regressor()`.
    *   **Deliverable:** Review all Roxygen examples, replacing deprecated function calls.
    *   **Done when:** Examples are up-to-date.

*   [x] **Ticket F-2: Update Vignette (`a_02_regressor.Rmd`)**
    *   **Scope:** Revise the regressor vignette.
    *   **Deliverable:** Vignette uses `regressor()`, `evaluate()`, `autoplot()`. Removes refs to deprecated functions.
    *   **Done when:** Vignette builds and reflects current API.

*   [x] **Ticket F-3: Add/Update Unit Tests**
    *   **Scope:** Ensure comprehensive testing.
    *   **Deliverable:** `tests/testthat/test-regressor.R` updated/expanded to cover constructor logic, evaluation parity/methods, shifting, plotting/printing.
    *   **Done when:** All tests pass; code coverage high.

---

**🧹 G. Cleanup & release**

*   [x] **Ticket G-1: Remove Orphaned Helpers & Comments** (was G-1)
    *   **Scope:** Delete unused internal code and comments.
    *   **Deliverable:** Remove `dots`, `conform_len`, old C++ wrappers (`evaluate_regressor_fast`, `evaluate_regressor_convolution` R functions), RANN comments from `eval_loop`.
    *   **Done when:** `R CMD check` clean.

*   [x] **Ticket G-2: Organize Files** (New, from Elegance 3.1)
    *   **Scope:** Group related functions into logical files.
    *   **Deliverable:** Code moved into files like `R/reg-constructor.R`, `R/reg-evaluate.R`, `R/evaluate-helpers.R`, `R/print-shift.R`, `R/utils-internal.R`.
    *   **Done when:** Files reorganized.
    
*   [x] **Ticket G-3: Check Namespacing/Imports** (New, from Correctness 1.3)
    *   **Scope:** Ensure all internal/external function calls are properly namespaced or imported.
    *   **Deliverable:** Review code (esp. calls like `make_hrf`, `evaluate.HRF`) and add `fmrireg::` or `@importFrom` tags as needed.
    *   **Done when:** `R CMD check` passes namespace checks.

*   [x] **Ticket G-4: Update NEWS and Version** (was G-2)
    *   **Scope:** Document changes.
    *   **Deliverable:** Entry in `NEWS.md`; Version bumped in `DESCRIPTION`.
    *   **Done when:** Changes documented.

*   [ ] **Ticket G-5: Final Benchmarks** (was G-3)
    *   **Scope:** Compare performance.
    *   **Deliverable:** Script in `bench/`; report confirms performance >= old.
    *   **Done when:** Performance verified.