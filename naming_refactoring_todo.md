# Naming Refactoring To-Do List

This checklist tracks the implementation of the design matrix column naming refactor based on the plan in `naming_refactoring.md`.

## Phase 0: Prep Work - Naming Utilities

- [X] **Create File:** Create `R/naming-utils.R`.
- [X] **Implement `zeropad()`:**
    - `zeropad(i, n_total)`
    - Returns `i` zero-padded to width `ceiling(log10(n_total+1))`.
    - Example: `zeropad(3, 15)` -> `"03"`, `zeropad(12, 150)` -> `"012"`.
    - Internal (`@keywords internal`).
- [X] **Implement `sanitize()`:**
    - `sanitize(x, allow_dot = TRUE)`
    - Calls `make.names(x, unique = FALSE)`. If `allow_dot = FALSE`, replaces `.` with `_`.
    - Handles multiple inputs vectorized.
    - Export (`@export`).
    - *Note:* Maintain backward compatibility for existing internal calls (`fmrireg:::sanitize`) until merge, possibly via re-export if needed temporarily.
- [X] **Implement `basis_suffix()`:**
    - `basis_suffix(j, nb)`
    - Returns `paste0("_b", zeropad(j, nb))`.
    - Vectorised in `j`.
    - Export (`@export`).
- [X] **Implement `make_unique_tags()`:**
    - `make_unique_tags(tags)`
    - Calls `make.unique(tags, sep = "#")`.
    - Internal (`@keywords internal`).
- [X] **Implement `make_term_tag()`:**
    - `make_term_tag(hrfspec, existing_tags = character())`
    - Generates tag from `hrfspec$id` or `paste(vapply(hrfspec$vars, rlang::as_label, ""), collapse = "_")`.
    - Sanitizes tag: `tag <- sanitize(tag, allow_dot = FALSE)`.
    - Makes unique using *sanitized* existing tags: `make_unique_tags(c(sanitize(existing_tags, FALSE), tag))[length(existing_tags) + 1L]`.
    - Internal (`@keywords internal`).
- [X] **Implement `level_token()`:**
    - `level_token(var, lev)`
    - Returns `paste0(sanitize(var, allow_dot = TRUE), ".", sanitize(lev, allow_dot = TRUE))`.
    - Internal (`@keywords internal`).
- [X] **Implement `continuous_token()`:**
    - `continuous_token(colname)`
    - Returns `sanitize(colname, allow_dot = TRUE)`.
    - Internal (`@keywords internal`).
- [X] **Implement `make_cond_tag()`:**
    - `make_cond_tag(tokens)`
    - Returns `paste(tokens, collapse = "_")`.
    - Internal (`@keywords internal`).
- [X] **Implement `add_basis()`:**
    - `add_basis(cond_tags, nb)`
    - If `nb == 1L`, return `cond_tags`.
    - Otherwise, return `as.vector(outer(cond_tags, basis_suffix(seq_len(nb), nb), paste0))`.
    - Internal (`@keywords internal`).
- [X] **Implement `make_column_names()`:**
    - `make_column_names(term_tag, cond_tags, nb)`
    - Add guard: `stopifnot(!grepl("__", term_tag))`.
    - Calls `cond_tags <- add_basis(cond_tags, nb)` internally.
    - Returns `paste0(term_tag, "_", cond_tags)`.
    - Internal (`@keywords internal`).
- [X] **Implement `is_valid_heading()` (Test Helper):**
    - `is_valid_heading(x)`
    - Uses `grepl("^[A-Za-z\\.][A-Za-z0-9\\._#]*$", x)` (allowing `#` for unique tags).
    - Place in `tests/testthat/helper-naming.R` (not exported).
- [X] **Add Imports:** Ensure `rlang` is imported if `as_label` is used.
- [X] **NAMESPACE:** Add `export(sanitize, zeropad, basis_suffix)`. Tag internal functions.
- [X] **Unit Tests:** Add basic unit tests for all helpers in `tests/testthat/test-naming-utils.R`.

## Phase 1: `realise_event_terms()`

- [X] **Locate:** Find the section in `event_model_helpers.R` where unique term names/tags are generated.
- [X] **Generate `uid`:** Add logic to generate sequential UIDs (e.g., `sprintf("t%02d", seq_along(term_list))`).
- [X] **Store `uid`:** Store the generated UID on *both* the `hrfspec` (e.g., `spec$uid <- uid`) and the resulting `event_term` (`attr(term, "uid") <- uid`).
- [X] **Generate `term_tag`:** Replace existing unique naming logic with a call to `make_term_tag(spec, sanitized_existing_tags)` (ensure `existing_tags` are sanitized with `allow_dot=FALSE` before passing).
- [X] **Store `term_tag`:** Store the result on *both* the `hrfspec` (e.g., `spec$term_tag <- term_tag`) and the resulting `event_term` (`attr(term, "term_tag") <- term_tag`).
- [X] **Backward Compatibility:** Do *not* remove the old `$varname` field from the term object.
- [ ] **Update Tests:** Ensure tests relying on old term names/attributes are updated or removed if obsolete.

## Phase 2: `columns.*()` Family

- [X] **`columns.Poly`:**
    - Modify to return `continuous_token(paste0("poly_", x$argname, "_", zeropad(1:x$degree, x$degree)))`. (Use `degree` for padding width).
    - Remove any logic related to `nbasis`.
- [X] **`columns.BSpline` (and similar basis functions like DCT):**
    - Modify to return `continuous_token(paste0("bs_", x$argname, "_", zeropad(1:x$degree, x$degree)))` (using appropriate prefix like `bs`, `dct`, and using `degree` for padding width).
    - Remove any logic related to `nbasis`.
- [X] **`columns.Scale`:**
    - Modify to return `continuous_token(paste0("z_", x$argname))`.
    - Constructor also updated to store this name.
- [X] **`columns.ScaleWithin`:**
    - Modify to return `continuous_token(paste0("z_", x$argname, "_by_", x$grpname))`.
    - Constructor also updated to store this name.
- [X] **`columns.RobustScale`:**
    - Modify to return `continuous_token(paste0("robz_", x$argname))`.
    - Constructor also updated to store this name.
- [X] **`columns.Standardized`:**
    - Ensure it returns `continuous_token(x$name)`.
    - Constructor also updated to store this name.
- [X] **`columns.Ident`:**
    - Ensure it returns one `continuous_token(var)` per variable if multiple inputs (e.g., `Ident(a, b)` -> `c("a", "b")` after tokenization).
- [X] **Update Tests:** Adjust tests checking specific output strings from these methods.
- [ ] *Note:* Remember related S3 methods like `nbasis.*` might need updating if `degree` handling changes internal column counts.

## Phase 3: `conditions.event_term()`

- [X] **Locate:** Find `conditions.event_term` in `conditions.R` (or similar).
- [X] **Implement Caching:** Check for `val <- attr(x, "..conds")` first; return `val` if found and options match. Store result in `attr(x, "..conds") <- result` before returning.
- [X] **Implement Shortcut:** If term is single continuous variable with 1 column (check `length(x$events) == 1` and `ncol(columns(x$events[[1]])) == 1`), return `columns(x$events[[1]])` directly (after potential basis expansion).
- [X] **Rewrite Token Generation:**
    - Remove old logic for assembling condition parts.
    - Use `lapply` over `x$events`.
    - If `is_categorical(ev)` (using existing helper), get `levels(ev$value)` and generate tokens using `level_token(ev$varname, levs)`.
    - If continuous, call `columns(ev)` to get basis tokens (already sanitized via `continuous_token`).
    - Store these component tokens.
- [X] **Combine Tokens:**
    - Use `expand.grid` on the list of component tokens.
    - Use `apply` with `make_cond_tag` to paste tokens with `_` separator.
    - This yields `base_cond_tags`.
- [X] **Handle `drop.empty`:** Removed flawed `drop.empty` logic.
- [X] **Handle `expand_basis`:**
    - If `expand_basis` is TRUE:
        - Get `nb <- nbasis(attr(x, "hrfspec")$hrf)`.
        - Call `final_cond_tags <- add_basis(base_cond_tags, nb)`.
    - Else:
        - `final_cond_tags <- base_cond_tags`.
- [X] **Return:** Return `final_cond_tags`.
- [X] **Remove Old Code:** Delete special case handling that the new logic replaces.
- [X] **Update Tests:** Test output with `expand_basis=FALSE` and `expand_basis=TRUE`, including interactions, caching, and shortcuts.

## Phase 4: `convolve.event_term()`

- [X] **Locate:** Find `convolve.event_term` in `convolve.R` (or similar).
- [X] **Remove Old Naming Logic:** Delete all existing code that manually creates column names, adds basis suffixes, or handles term tags/prefixes.
- [X] **Get Inputs:**
    - `term_tag <- attr(x, "term_tag") %||% stop("term_tag attribute not found on event_term")`.
    - `nb <- nbasis(hrf)`.
    - Use `colnames(design_matrix(...))` instead.
- [X] **Generate Column Names:**
    - `cn <- make_column_names(term_tag, base_cnames, nb)`.
- [X] **Assign Column Names:**
    - `colnames(cmat) <- cn`.
- [X] **Add Debug Validation:**
    - `if (getOption("fmrireg.debug", FALSE)) { stopifnot(all(is_valid_heading(cn))) }` (requires `is_valid_heading` accessible).
- [X] **Update Tests:** Ensure tests check the final column names assigned to the convolved matrix match the new `term_tag_condition_tag[_b##]` format.

## Phase 5: `build_event_model_design_matrix()`

- [X] **Locate:** Find `build_event_model_design_matrix` function.
- [X] **Remove Name Handling:** Delete any lines that modify or assign `colnames(dm_mat)` after `cbind`-ing term matrices.
- [X] **Verify `col_indices`:** Ensure `attr(dm, "col_indices")` still correctly records column ranges based on the incoming (already correctly named) term matrices.
- [X] **Set Final Name Flag:** Add `attr(dm, "colnames_final") <- TRUE` after `cbind`.
- [ ] **Add Optional Debug Check:** Add `if (getOption("fmrireg.debug", FALSE)) { old_names <- colnames(dm); ... later ...; stopifnot(identical(colnames(dm), old_names)) }` around potentially problematic downstream code sections within this function.
- [X] **Update Tests:** Verify the final design matrix has the correct column names inherited from `convolve.event_term`.

## Phase 6: Alias Cleanup

- [X] **Review `columns.Standardized`, `columns.Scaled` etc.:** Check if any contain duplicated logic for string manipulation that is now handled by the primary `columns.*` methods and `continuous_token`. Remove redundancy.
- [X] **Review `levels.*` methods:** Check for similar redundant string formatting logic and remove.

## Phase 6.5: Contrast Engine Update

- [X] **Add `translate_legacy_pattern()` Helper:** 
    - Define internal helper `translate_legacy_pattern(pattern)` in `R/contrast.R`.
    - Should replace `:basis\\[(\d+)\\]` with `_b\1` (or similar robust pattern).
    - Should replace `:` with `_` (for interactions).
    - Consider order of replacements.
- [X] **Modify `.col_index()`:**
    - Remove internal logic checking `nbasis` to decide on `expand_basis`.
    - Call `pattern <- translate_legacy_pattern(pattern)` at the start.
    - Ensure `conditions(..., expand_basis=TRUE)` is used when appropriate (likely already correct).
- [X] **Update Internal `paste(collapse=":")`:**
    - Search `R/contrast.R` for `paste(..., collapse = ":")` used for internal keys/rownames.
    - Change separator to `_`.
- [X] **Remove Obsolete `gsub`:**
    - Delete the `gsub(":", ".", ...)` call within `contrast_weights.contrast_formula_spec`.
- [X] **Remove Manual Basis Broadcasting:**
    - In relevant `contrast_weights.*` methods (e.g., `pair_contrast_spec`), remove code like `rep(weights, each = nbasis)`.
    - Ensure weights are calculated relative to base conditions and applied via name matching to potentially expanded conditions.
- [X] **Update Roxygen Examples:**
    - In `R/contrast.R`, update regex examples in function docs (e.g., `column_contrast`) to reflect new syntax (`_b03`, `Var.Level_Mod`).
- [X] **Add Documentation Section:**
    - Add `@section Column-name grammar...` to `R/contrast.R` roxygen block.
- [X] **Export Helper Alias:**
    - Add `@export translate_legacy_pattern` to its roxygen block.
- [X] **Update Contrast Tests:**
    - Adjust regex patterns and expected names in `tests/testthat/test-contrast-column.R` (or relevant files like `test_event_model.R`).
    - Review and update relevant snapshot tests (e.g., for `conditions()` output if used in contrast tests).
    - (Skipped specific illegal character test case).

## Phase 7: Tests & Snapshots

- [X] **Golden Heading Test:**
    - Create a test model with diverse terms (factors, poly, interactions, scaled vars).
    - Save `colnames(design_matrix(model))` to a snapshot/reference.
    - Ensure test fails if names change unexpectedly.
- [X] **Clash Test:**
    - Use `event_model(on ~ hrf(x, subset = x == 1) + hrf(x, subset = x == 2), ...)`.
    - Assert that resulting column names start with `x_`

## Phase 8: Contrast Engine Simplification (Post-Legacy Deprecation)

*This phase assumes the core naming refactor (Phases 0-7) is complete, stable, and a release cycle with deprecation warnings for legacy contrast patterns (from Phase 6.5) has passed.*

- [ ] **Remove `translate_legacy_pattern()`:** Delete the helper function entirely from `R/contrast.R`.
    - Add a check for reverse dependencies (users directly importing `fmrireg:::translate_legacy_pattern`).
- [ ] **Remove Legacy Checks in `.col_index()`:** Delete the call to `translate_legacy_pattern()` and any associated messaging. Ensure it consistently uses `conditions()` with appropriate `expand_basis`.
    - Change return guarantee: Instead of `integer(0)` on failure/no match, ensure it throws an explicit error if inputs are invalid (though `grep` naturally returns `integer(0)` for no matches, which is fine).
- [ ] **Remove Duplicated `expand_basis` Logic:**
    - Remove internal logic determining `expand_basis` based on HRF specs within `.col_index`, `contrast_weights.column_contrast_spec`, etc. Rely on explicit arguments or the shared `condnames()` helper.
    - *Confirm:* Ensure `contrast_weights.pair_contrast_spec` (or its helpers) is the *only* place checking `nbasis(hrf)` directly for basis broadcasting logic.
- [ ] **Implement Shared `condnames()` Helper:**
    - Create an internal helper `condnames(term, expanded = TRUE/FALSE)` that reliably calls `conditions(term, drop.empty = FALSE, expand_basis = expanded)`. Mark `@keywords internal`.
- [ ] **Implement Shared `mask_to_weights()` Helper:**
    - Create an internal helper `mask_to_weights(names, A_mask, B_mask = NULL)` based on the minimal API example, handling +1/nA and -1/nB logic.
    - Add an assertion/check that `abs(sum(weights)) < tol` (e.g., `tol = 1e-8`) for A-B type weights (where `B_mask` is not NULL and `sum(B_mask) > 0`).
- [ ] **Simplify `contrast_weights.column_contrast_spec`:**
    - Refactor to:
        - Get potentially expanded condition names using `condnames(term, expanded = TRUE)`.
        - Use `grep` directly on these names for `pattern_A` and `pattern_B`.
        - *Retain:* Check for overlapping indices between A and B matches (`idx_A %in% idx_B`) and `stop()` if overlap occurs.
        - Use the shared `mask_to_weights()` helper to calculate final weights.
        - Structure and return the result list.
- [ ] **Simplify `contrast_weights.pair_contrast_spec`:**
    - Refactor to:
        - Get base condition names using `condnames(term, expanded = FALSE)`.
        - Get `cells()` data frame.
        - Apply `where` clause if present to filter `cells_df`.
        - Evaluate `A` and `B` formulas on the (filtered) `cells_df`.
        - Map logical results back to base condition names using `apply(filtered_cells_df, 1, make_cond_tag)` and `%in%`.
        - Use `mask_to_weights()` to calculate base weights (vector named by base condition names).
        - Determine if basis expansion is needed (`nb <- nbasis(...)`).
        - If `nb > 1`:\n            - Get expanded condition names using `condnames(term, expanded = TRUE)`.\n            - Create a zero vector/matrix for expanded names.\n            - Map the calculated base weights to the corresponding expanded names efficiently (e.g., using `match(sub(\\"_b\\\\\\\\d+$\\", \\"\\", expanded_names), names(base_weights))` to find indices and assign `base_weights`).\n            - Add a `warning` if any non-zero `base_weights` do not find a corresponding match in `expanded_names` (signals potential internal mismatch).
        - Else (`nb == 1`):\n            - Use the base weights directly.\n        - Structure and return the result list.
- [ ] **Remove Error Fallback Logic:** Replace `tryCatch` blocks or similar logic (that previously returned defaults on error) with direct calls or `stop()`, *except* around user-provided formula evaluation (`eval_tidy`), where a `tryCatch` converting errors to informative `stop()` messages should be retained.
- [ ] **Review Other `contrast_weights.*` Methods:** Check `oneway_contrast`, `interaction_contrast`, `poly_contrast` etc. for simplification opportunities (removing redundant checks, using `cells()`/`conditions()` directly).
    - *Example:* `oneway_contrast` might reuse `mask_to_weights()` internally for generating each contrast column comparing a level to others.
- [ ] **Update Contrast Tests:**
    - Remove tests specifically targeting legacy pattern translation or behavior.
    - Ensure tests cover the simplified logic and matching against canonical names (both base and expanded where relevant).
    - Add expectations checking `sum(weights) == 0` (within tolerance) for all A-B type contrasts (e.g., `pair_contrast`, `column_contrast` with A and B).
    - Update snapshots if necessary.
- [ ] **Update Documentation:**
    - Remove references to legacy syntax support and translation from `?contrast`, `?column_contrast`, `?pair_contrast`, etc.
    - Ensure all examples use the canonical `term_tag_condition_tag[_b##]` grammar.
    - Verify `?fmrireg_naming` help page exists and is accurate.
    - *Check Vignettes:* Search `vignettes/*.Rmd` files for examples using legacy regex patterns (e.g., `:basis\\[`, `Var[Level]`) and update them.

## Phase 9: Naming Refinements & Bug Fixes (Post-Core Refactor)

*This phase addresses specific naming issues and refines the logic after the main refactor (Phases 0-8) is stable.*

- [ ] **Refine `term_tag` Generation for `Ident()`-only Terms:**
    - **Issue:** `hrf(Ident(var1, var2))` can lead to `term_tag` like `Ident_var1_var2`, and final columns `Ident_var1_var2_var1`, `Ident_var1_var2_var2`.
    - **Goal:** Achieve simpler names, e.g., `var1`, `var2` if `id` is not provided, or `customID_var1`, `customID_var2` if `hrf(Ident(var1,var2), id="customID")`.
    - **Plan:**
        - Modify `make_term_tag()` in `R/naming-utils.R`.
        - If `hrfspec$id` is NULL *and* the *only* variable in `hrfspec$vars` is an `Ident()` call *and* that `Ident()` call itself contains multiple simple variable names:
            - Consider *not* prepending "Ident" to the `term_tag`.
            - The `term_tag` could be an empty string `""` in this specific scenario, or based on the `id` if provided.
            - If the `term_tag` becomes empty, `make_column_names` would effectively produce `condition_tag[_b##]`.
            - Alternatively, if `hrfspec$id` is NULL, and the term is `Ident(V1, V2, ...)`, the `term_tag` could be derived from `V1_V2` (names within Ident), rather than `Ident_V1_V2`.
        - If `hrfspec$id` is provided (e.g., `hrf(Ident(RT1,RT2), id="myterm")`), the `term_tag` should be based on `myterm` (e.g. `myterm`), resulting in `myterm_RT1`, `myterm_RT2`.
        - Ensure `columns.Ident` continues to return the bare variable names (e.g., `"RT1"`, `"RT2"`) as the `condition_tag` components.
        - Ensure uniqueness is still handled correctly if multiple such `Ident()` terms exist.
    - **Tests:** Add tests for `Ident(var1,var2)`, `Ident(var1,var2, id="myid")`, and multiple `Ident()` terms to verify clear and non-redundant naming.

- [ ] **Clarify and Test Naming for Combined Basis Functions:**
    - **Issue:** Ensure clarity in how names are formed when an input variable basis function (e.g., `Poly(RT,2)`) is combined with an HRF basis expansion (e.g., `HRF_SPMG2`).
    - **Goal:** Verify the `term_tag_PolyName.Var_k_b##` structure is correctly and consistently generated.
    - **Plan:**
        - This is mostly a verification and testing task, as the current plan in `naming_refactoring.md` *should* handle this correctly.
        - Review `columns.Poly`, `columns.BSpline` etc. to ensure they only produce the `PolyName.Var_k` part (e.g., `poly_RT_01`).
        - Review `convolve.event_term` and `make_column_names` to confirm they correctly append the `_b##` (from HRF) to the `condition_tag` (which itself might be `poly_RT_01`).
    - **Tests:** Add specific tests for models like `event_model(onset ~ hrf(Poly(RT,2), basis="spmg2"), ...)` and `event_model(onset ~ hrf(BSpline(Age,3), basis="gamma3"), ...)` to assert the expected column names (e.g., `term_poly_RT_01_b01`, `term_poly_RT_01_b02`, `term_poly_RT_02_b01`, `term_poly_RT_02_b02`).