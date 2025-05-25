# Regressor Module Refactoring Plan

**Goal:** Refactor the `R/regressor.R` module to create a unified, efficient, and maintainable system for representing and evaluating fMRI regressors, leveraging the refactored HRF pipeline (`as_hrf`, decorators, `evaluate.HRF`).

**Guiding Principles:**

1.  **Single Core Constructor:** Implement a primary internal constructor (e.g., `Reg()`) responsible for creating regressor objects, handling input validation, recycling, and storing data efficiently.
2.  **Backward Compatibility:** Maintain the existing `regressor()` function as the primary public interface. Internally, `regressor()` will call the new core constructor (`Reg()`) and potentially add the `"regressor"` class tag for compatibility with `inherits(x, "regressor")` checks.
3.  **Deprecate Redundancy:** Soft-deprecate `single_trial_regressor()` and `null_regressor()`, pointing users towards `regressor()` (which now handles these cases via `Reg()`).
4.  **Unified Evaluation:** Implement a single `evaluate.Reg()` S3 method that handles different evaluation engines (FFT, direct convolution, R loop) and delegates HRF evaluation to the refactored `evaluate.HRF`.
5.  **Composable HRF:** The regressor system should rely entirely on the public HRF interface (`evaluate.HRF`, `nbasis`, attributes) without knowledge of internal HRF implementation details.
6.  **Performance & Memory:** Optimize evaluation speed and memory usage, potentially using memoization, efficient C++ routines, and optional sparse matrix representations.
7.  **Modernization:** Update plotting (`autoplot`) and printing (`print`) methods using modern packages like `ggplot2`, `cli`, and `glue`.

**High-Level Phases:**

1.  **Core Implementation:** Define the `Reg` class structure and its core constructor. Implement input validation and recycling.
2.  **Public API & Compatibility:** Make `regressor()` a wrapper for `Reg()`, add deprecation shims for other constructors. Define the `evaluate.Reg` generic.
3.  **Evaluation Engine:** Implement the different evaluation methods (`fft`, `conv`, `loop`) within `evaluate.Reg`, leveraging C++ where appropriate and ensuring delegation to `evaluate.HRF`.
4.  **Optimization:** Add performance enhancements like memoization and consider sparse matrix support.
5.  **Helpers & Documentation:** Update plotting, printing, tests, and documentation to reflect the new structure and capabilities.
6.  **Cleanup:** Remove redundant internal functions after verification.

This refactoring will align the regressor system with the elegance and efficiency of the new HRF architecture, providing a cleaner codebase while maintaining user-facing stability. 