Okay, let's scope Sprint 2 with the goal of producing a runnable `fmri_model` from the `fmri_config` object generated in Sprint 1. This involves translating the validated DSL specifications into `fmrireg`'s internal R objects (like `hrfspec`, `event_term`, `baseline_model`, etc.) and assembling them.

## fmrireg DSL Sprint 2: `fmri_config` to Runnable `fmri_model`

**Sprint Goal:** Enable the `build_fmri_model_from_config(fmri_config, subject_id)` function to successfully generate a complete `fmri_model` object for a given subject, ready for use with `fmri_lm()` or other fitting functions. This includes processing transformations, HRFs, terms, baseline, and nuisance regressors.

---

**Ticket Prefix:** `DSL2-` (Distinguishing from Sprint 1)

---

### **Epics:**

*   **EPI-DSL2-DATA:** Subject-Specific Data Loading & Preparation
*   **EPI-DSL2-HRF:** HRF Instantiation from DSL
*   **EPI-DSL2-VARPROC:** Variable Transformation & Basis Expansion (from DSL `transformations` and `terms.*.modulator_basis`)
*   **EPI-DSL2-TERMS:** `hrfspec` & `covariatespec` Generation from DSL `terms`
*   **EPI-DSL2-MODEL:** `event_model`, `baseline_model`, and `fmri_model` Assembly

---

### **Tickets for Sprint 2:**

**Subject-Specific Data Loading & Preparation (EPI-DSL2-DATA)**

1.  **DSL2-101: Implement `load_and_prepare_subject_data(fmri_config, subject_id)` (Refined)**
    *   **Description:** Enhance/refine the existing helper to:
        *   Use `fmri_config$project` (bidser object), `subject_id`, `fmri_config$tasks`, `fmri_config$runs` to load specific event files (`bidser::read_events`).
        *   Apply run filtering based on `fmri_config$runs`.
        *   Load relevant confound files (`bidser::read_confounds`) based on columns needed by `fmri_config$variables` with `role: NuisanceSource` and selected by `fmri_config$confound_groups`.
        *   Resolve `TR` and `run_length` for each run using `fmri_config$dataset$scan_params` (defaults and overrides), falling back to BIDS metadata if available, and erroring if ambiguous.
        *   Combine event data across selected runs for the subject into a single data frame.
        *   Combine selected confound data across runs, ensuring alignment with scan timing.
    *   **AC:**
        *   Returns a list with `events_df` (all events for subject), `confounds_df` (all selected confounds for subject, aligned), `run_lengths` (vector), `TR` (scalar).
        *   Handles cases where some runs/tasks might be missing for a subject gracefully (e.g., warning, skipping those).
        *   Correctly applies `TR`/`run_length` overrides.
    *   **Estimate:** 4 points

**HRF Instantiation from DSL (EPI-DSL2-HRF)**

1.  **DSL2-201: Implement `get_hrf_from_dsl(hrf_name, fmri_config)`**
    *   **Description:** Creates an `fmrireg::HRF` object based on a named HRF definition in `fmri_config$hrfs`.
    *   **AC:**
        *   Handles `SPMCanonical` (maps to `HRF_SPMG1`).
        *   Handles `SPMCanonicalDerivs` (maps to `HRF_SPMG2` or `HRF_SPMG3` based on `derivatives` array).
        *   Handles `GammaFunction`, `Gaussian` by calling `gen_hrf()` with parameters from DSL.
        *   Handles `BSplineBasisHRF`, `TentBasisHRF`, `FourierBasisHRF`, `DaguerreBasisHRF` by calling appropriate `fmrireg` generator functions (e.g., `hrfspline_generator`) with parameters.
        *   Handles `CustomR`: retrieves function name/path from `definition` and attempts to get/source it.
        *   Applies global `lag`, `width`, `summate`, `normalize` if specified in the DSL HRF definition.
        *   Returns a valid `fmrireg::HRF` object.
    *   **Estimate:** 4 points

**Variable Transformation & Basis Expansion (EPI-DSL2-VARPROC)**

1.  **DSL2-301: Process `variables`, `transformations`, and `modulator_basis` into Unified Model Environment**
    *   **Description:** **[DECISION: Option A - Basis as Transform]** Create a unified environment-based approach that holds all "model variables" available for terms. This involves:
        *   Starting with raw variables defined in `fmri_config$variables` (mapped from BIDS columns in `subject_data$events_df` or `subject_data$confounds_df`).
        *   Iteratively applying operations defined in `fmri_config$transformations` to create derived variables.
        *   Processing `modulator_basis` specifications from `fmri_config$terms` as additional transformation steps, creating basis-expanded derived variables (e.g., `rt_poly_deg2` from `Poly(rt, degree=2)`).
        *   Store all variables (raw, transformed, and basis-expanded) in a unified environment with unique names.
    *   **Implementation Details:**
        *   Use `fmrireg`'s existing `ParametricBasis` classes (`Poly`, `BSpline`, `Scale`, etc.) from `basis.R`.
        *   For each term with `modulator_basis`, create a `ParametricBasis` object and extract its `$y` matrix as a derived variable.
        *   Derived variable naming: `{mod_var}_{basis_type}_{params}` (e.g., `rt_poly_deg2`, `difficulty_bspline_deg3`).
    *   **AC:**
        *   All variables listed in `fmri_config$variables` are accessible using their "model variable name" key.
        *   All derived variables from `fmri_config$transformations` are computed and accessible.
        *   All basis-expanded variables from `terms.*.modulator_basis` are computed and accessible as matrix-valued variables.
        *   Handles dependencies (a transformation using the output of a previous one).
        *   Applies simple transform ops: `"center"`, `"scale-sd"`, `"z-score"`, `"log"`, `"exp"`, `"factorize"`.
        *   Implements complex transform ops: `"scale-within-group"`, `"clip"`, `"recode-levels"`.
    *   **Estimate:** 6 points (increased due to basis integration)

**`hrfspec` & `covariatespec` Generation from DSL `terms` (EPI-DSL2-TERMS)**

1.  **DSL2-401: Convert DSL `terms` to `hrfspec` and `covariatespec` Objects**
    *   **Description:** **[DECISION: Use existing covariate() infrastructure for NuisanceRegressors]** For each term defined in `fmri_config$model$terms`:
        *   Retrieve the full term definition from `fmri_config$terms`.
        *   Based on `term$type`:
            *   **EventRelated, ParametricModulation, Trialwise**: Create `fmrireg::hrfspec` objects using the existing `hrf()` function.
            *   **NuisanceRegressors**: Create `fmrireg::covariatespec` objects using the existing `covariate()` function.
    *   **Details for `hrfspec` creation:**
        *   `vars`: Collect quosures of the relevant *model variable names* from the unified environment (DSL2-301). For `ParametricModulation` with `modulator_basis`, reference the basis-expanded derived variable (e.g., `rt_poly_deg2`).
        *   `hrf`: Get the `HRF` object using `get_hrf_from_dsl(term$hrf, fmri_config)`. Apply `term$lag` to this HRF.
        *   `id`: Use the term's key from the DSL.
        *   `subset`: Parse `term$subset` string into an expression.
        *   `contrasts`: **[DECISION: Store contrast names, resolve later]** Store the names of contrasts from `fmri_config$model$contrasts` that reference this term. The existing `contrast_weights.event_model()` system will resolve these names later.
    *   **Details for `covariatespec` creation (NuisanceRegressors):**
        *   Use `fmrireg::covariate()` function with variables from `term$nuisance_source_variables`.
        *   `id`: Term's key.
        *   This leverages the existing, tested infrastructure for non-convolved regressors.
    *   **AC:**
        *   Returns a named list of `hrfspec` and `covariatespec` objects, one for each term in `fmri_config$model$terms`.
        *   `hrfspec$vars` correctly refers to processed variables from the unified environment.
        *   Basis-expanded modulators are correctly referenced as derived variables.
        *   HRF objects are correctly instantiated and lagged.
        *   Subset expressions are captured.
        *   Contrast names are stored for later resolution.
    *   **Estimate:** 4 points (reduced due to clearer decisions)

**`event_model`, `baseline_model`, and `fmri_model` Assembly (EPI-DSL2-MODEL)**

1.  **DSL2-501: Construct `baseline_model` from DSL**
    *   **Description:**
        *   Parse `fmri_config$model$baseline$basis` (e.g., "BSpline(3)" or structured) to determine type and degree.
        *   Resolve `fmri_config$model$baseline$include_confound_groups` using `fmri_config$confound_groups` to identify specific confound columns.
        *   Extract these selected confound columns from `subject_data$confounds_df`.
        *   Call `fmrireg::baseline_model()` with the derived basis, degree, intercept type, `sampling_frame` (from DSL2-101), and the `nuisance_list` (confound matrices split by run).
    *   **AC:**
        *   A valid `fmrireg::baseline_model` object is created.
        *   Nuisance regressors are correctly included if specified.
        *   Shorthand basis strings (e.g., "BSpline(3)") are correctly parsed.
    *   **Estimate:** 3 points

2.  **DSL2-502: Construct `event_model` from Processed Term Specs**
    *   **Description:** **[DECISION: Leverage existing covariate() integration]**
        *   Take the mixed list of `hrfspec` and `covariatespec` objects from DSL2-401.
        *   Call `fmrireg::event_model()` with the combined list, as the existing system already handles mixed term types.
        *   Pass `data` as the unified variable environment from DSL2-301, `block` column from subject data, and `sampling_frame`.
        *   The existing `event_model` constructor will handle both convolved (`hrfspec`) and non-convolved (`covariatespec`) terms appropriately.
    *   **AC:**
        *   A valid `fmrireg::event_model` object is created.
        *   Design matrix includes both convolved event-related regressors and direct nuisance regressors.
        *   Column naming follows the established `fmrireg` conventions.
    *   **Estimate:** 2 points (reduced due to leveraging existing infrastructure)

3.  **DSL2-503: Assemble `fmri_model`**
    *   **Description:** Call `fmrireg::fmri_model(event_model, baseline_model)`.
    *   **AC:**
        *   A valid `fmrireg::fmri_model` object is returned.
    *   **Estimate:** 0.5 points

4.  **DSL2-504: Main Function `build_fmri_model_from_config(fmri_config, subject_id)`**
    *   **Description:** Orchestrate all steps from DSL2-101 to DSL2-503 for a given `subject_id`.
    *   **AC:**
        *   Takes a validated `fmri_config` and `subject_id`.
        *   Returns a runnable `fmri_model` object.
        *   Handles errors gracefully with informative messages.
        *   Provides progress feedback for long-running operations.
    *   **Estimate:** 1 point

---

### **Resolved Key Decisions:**

1.  **Parametric Modulator Basis Expansion:** **Option A (Basis as Transform)** - Process `modulator_basis` in DSL2-301 as transformation steps, creating derived variables that `hrfspec` objects can reference. This leverages `fmrireg`'s existing `ParametricBasis` infrastructure and provides clean separation of concerns.

2.  **Integration of `type: NuisanceRegressors` terms:** **Option 2 (Use existing `covariate()` infrastructure)** - Create `covariatespec` objects using `fmrireg::covariate()` for terms with `type: NuisanceRegressors`. This leverages the existing, tested infrastructure for non-convolved regressors.

3.  **Contrast Handling:** **Store contrast names in `hrfspec`, resolve during `contrast_weights.event_model`** - Store contrast names from the DSL in `hrfspec` objects and let the existing `contrast_weights.event_model()` system resolve them using the global contrast definitions.

---

### **Implementation Architecture:**

**Unified Variable Processing Pipeline (DSL2-301):**
```r
process_variables_and_transformations <- function(fmri_config, subject_data) {
  model_env <- new.env()
  
  # Add raw variables from BIDS columns
  for (var_name in names(fmri_config$variables)) {
    var_def <- fmri_config$variables[[var_name]]
    model_env[[var_name]] <- extract_bids_column(subject_data, var_def$bids_column)
  }
  
  # Process transformations in dependency order
  for (derived_name in names(fmri_config$transformations)) {
    transform_def <- fmri_config$transformations[[derived_name]]
    source_data <- model_env[[transform_def$source_variable]]
    model_env[[derived_name]] <- apply_transform_ops(source_data, transform_def$ops)
  }
  
  # Process modulator_basis from terms (Option A)
  for (term_name in names(fmri_config$terms)) {
    term <- fmri_config$terms[[term_name]]
    if (!is.null(term$modulator_basis)) {
      basis_obj <- create_basis_from_dsl(term$modulator_basis, model_env[[term$mod_var]])
      basis_var_name <- paste0(term$mod_var, "_", term$modulator_basis$type, "_deg", term$modulator_basis$parameters$degree %||% "")
      model_env[[basis_var_name]] <- basis_obj$y  # Store the expanded matrix
    }
  }
  
  return(model_env)
}
```

**HRF Mapping (DSL2-201):**
```r
get_hrf_from_dsl <- function(hrf_name, fmri_config) {
  if (hrf_name == "canonical") {
    return(HRF_SPMG1)  # Default canonical HRF
  }
  
  hrf_def <- fmri_config$hrfs[[hrf_name]]
  
  switch(hrf_def$type,
    "SPMCanonical" = HRF_SPMG1,
    "SPMCanonicalDerivs" = if("Temporal" %in% hrf_def$derivatives) HRF_SPMG2 else HRF_SPMG3,
    "BSplineBasisHRF" = hrfspline_generator(degree = hrf_def$parameters$degree),
    # ... other mappings
  )
}
```

---

**Definition of Done for Sprint 2:**

*   All tickets (DSL2-101 to DSL2-504) are completed and merged.
*   `build_fmri_model_from_config(config, subject_id)` successfully produces an `fmri_model` object that:
    *   Contains correctly convolved event-related regressors based on DSL `terms` of type `EventRelated`, `ParametricModulation`, and `Trialwise`.
    *   Correctly incorporates basis expansions for parametric modulators as derived variables.
    *   Includes regressors from `terms` of type `NuisanceRegressors` *without* HRF convolution using the `covariate()` infrastructure.
    *   Includes a correctly formed `baseline_model` with drift, intercept, and selected `confound_groups`.
    *   The `design_matrix(fmri_model)` output has correctly named columns according to the `fmrireg` naming conventions.
*   Unit tests cover the instantiation of different HRF types, variable transformations, term processing (especially parametric modulation with basis), and baseline construction.
*   The system handles at least one common BIDS dataset example defined via the v2.4 DSL.
*   Documentation for `build_fmri_model_from_config` is updated.
*   Integration tests demonstrate that the resulting `fmri_model` can be successfully used with `fmri_lm()`.

---

### **Stretch Goals (If time permits, or for next sprint):**

*   **DSL2-601: Implement R code generation for `CustomR` HRFs**
    *   **Description:** Safely evaluate or source the `definition` string from `hrfs.*.definition` for `type: CustomR`.
    *   **Implementation:** Use `eval(parse(text=definition))` with proper error handling and sandboxing.
    *   **AC:** `get_hrf_from_dsl()` can handle `CustomR` HRF definitions by executing user-provided R code.
    *   **Estimate:** 2 points
    *   **Note:** This extends DSL2-201 beyond the basic HRF types.

*   **DSL2-602: Implement parser for shorthand `basis` strings (e.g., "BSpline(3)")**
    *   **Description:** Parse shorthand basis specifications like `"BSpline(3)"`, `"Polynomial(2)"` in `baseline.basis` and convert them to structured form.
    *   **Implementation:** Create a regex-based parser that extracts type and parameters from shorthand strings.
    *   **AC:** DSL2-501 can handle both structured basis definitions and shorthand strings interchangeably.
    *   **Estimate:** 1.5 points
    *   **Note:** This enhances DSL2-501 with more user-friendly syntax.

*   **DSL2-603: Advanced `transformations.ops` (e.g., `Clip`, `RecodeLevels`)**
    *   **Description:** Implement additional transformation operations beyond the basic set:
        *   `clip`: Constrain values to a min/max range
        *   `recode-levels`: Map factor levels to new values
        *   `winsorize`: Replace extreme values with percentile cutoffs
        *   `rank-transform`: Convert to ranks
    *   **AC:** DSL2-301 supports the extended set of transformation operations with proper parameter validation.
    *   **Estimate:** 3 points
    *   **Note:** This extends DSL2-301 beyond the basic transformation operations.

*   **DSL2-604: Enhanced Error Handling and Diagnostics**
    *   **Description:** Implement comprehensive error reporting and diagnostic information:
        *   Detailed error messages with suggestions for common mistakes
        *   Validation warnings for potentially problematic configurations
        *   Progress reporting for long-running operations
        *   Debug mode with intermediate object inspection
    *   **AC:** All DSL2 functions provide clear, actionable error messages and optional diagnostic output.
    *   **Estimate:** 2 points

*   **DSL2-605: Performance Optimization for Large Datasets**
    *   **Description:** Optimize the DSL pipeline for large BIDS datasets:
        *   Lazy loading of confound data (only load needed columns)
        *   Parallel processing of transformations where possible
        *   Memory-efficient handling of basis expansions
        *   Caching of expensive computations
    *   **AC:** `build_fmri_model_from_config()` handles large datasets efficiently with optional parallel processing.
    *   **Estimate:** 4 points

*   **DSL2-606: Advanced Contrast Resolution**
    *   **Description:** Implement full contrast resolution from DSL specifications:
        *   Parse contrast formulas and map them to design matrix columns
        *   Handle complex contrast types (polynomial, interaction, etc.)
        *   Validate contrast specifications against available terms
        *   Generate contrast weight matrices automatically
    *   **AC:** The `fmri_model` object includes fully resolved contrasts that can be used directly with `fmri_lm()`.
    *   **Estimate:** 5 points
    *   **Note:** This completes the contrast handling that was deferred in DSL2-401.

---

**Notes on Sprint 1 Stretch Goals Integration:**

*   **DSL-401 (CustomR HRFs):** → **DSL2-601** (New stretch goal)
*   **DSL-402 (Shorthand basis parsing):** → **DSL2-602** (New stretch goal) 
*   **DSL-403 (Advanced transformations):** → **DSL2-603** (New stretch goal)
*   **DSL-404 (build_fmri_model_from_config):** → **DSL2-504** (Main Sprint 2 ticket - already covered)

This sprint leverages the existing, well-tested `fmrireg` infrastructure while providing the clean abstractions needed for the DSL, resulting in a more robust and maintainable implementation.