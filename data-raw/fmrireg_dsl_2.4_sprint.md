Okay, here's a sprint plan broken down into tickets for implementing the fmrireg DSL v2.4. This assumes a sprint-based agile workflow.

**Sprint Goal:** Deliver a robust parser and validator for fmrireg DSL v2.4, enabling the construction of `fmri_config` objects that can subsequently be used to build `fmri_model` instances. Initial focus is on parsing, validation, and core object creation, with model building (`build_fmri_model_from_config`) as a stretch or follow-on sprint.

---

## fmrireg DSL v2.4 Implementation Sprint

**Ticket Prefix:** `DSL-`

---

### **Epics:**

*   **EPI-DSL-CORE:** Core DSL Parsing & Validation Infrastructure
*   **EPI-DSL-SEMANTIC:** Semantic Validation & BIDS Integration
*   **EPI-DSL-BUILD:** `fmri_config` Object Construction

---

### **Tickets for the Sprint:**

**Core Parsing & Validation Infrastructure (EPI-DSL-CORE)**

1.  **DSL-101: Setup YAML Parsing & Basic Structure Validation**
    *   **Description:** Implement function `parse_yaml_to_list(yaml_file_path)` that reads the YAML file into an R list. Include basic error handling for malformed YAML.
    *   **Acceptance Criteria (AC):**
        *   Function reads a valid v2.4 YAML file into a nested R list.
        *   Throws an error if YAML is syntactically invalid.
        *   Handles file-not-found errors.
    *   **Estimate:** 1 point

2.  **DSL-102: Implement `ValidationErrors` R6 Class**
    *   **Description:** Create or refine the `ValidationErrors` R6 class as outlined (add\_error, is\_valid, format\_errors, stop\_if\_invalid).
    *   **AC:**
        *   Class can collect multiple errors with paths.
        *   `is_valid()` works correctly.
        *   `format_errors()` produces readable output.
        *   `stop_if_invalid()` stops execution with formatted errors.
    *   **Estimate:** 1 point

3.  **DSL-103: Implement Schema Check Helpers (`check_required`, `check_type`, `check_enum`, `check_pattern`)**
    *   **Description:** Implement the internal helper functions for schema validation based on the v2.4 spec (e.g., `check_type` handles `array`, `items`, `oneOf` for simple cases like `relpath`). Use JSON-Schema keywords where applicable.
    *   **AC:**
        *   Each helper correctly validates its specific constraint against a list node.
        *   Errors are added to the `ValidationErrors` object with correct paths.
        *   `allow_null` (or equivalent logic for optional fields not in `required` arrays) is handled.
    *   **Estimate:** 3 points

4.  **DSL-104: Implement `apply_defaults(data, defaults)` Function**
    *   **Description:** Create a function to recursively merge default values into the parsed YAML list.
    *   **AC:**
        *   Correctly applies defaults defined in the DSL spec (e.g., `events.onset_column = "onset"`).
        *   Handles nested default structures.
        *   User-provided values override defaults.
    *   **Estimate:** 2 points

5.  **DSL-105: Schema Validation - Top-Level Sections & `dataset` Block**
    *   **Description:** Implement schema validation for `dataset`, `events`, `variables`, `terms`, `models` top-level presence. Validate the `dataset` block structure and types, including `scan_params` flattening.
    *   **AC:**
        *   Correctly validates presence of required top-level keys.
        *   Validates `dataset.path`, `dataset.relpath` (string or array).
        *   Validates `dataset.subjects`, `dataset.tasks`, `dataset.runs` (types, patterns).
        *   Validates simplified `dataset.scan_params.TR`, `TR_overrides`, `run_lengths`, `run_length_overrides`.
    *   **Estimate:** 3 points

6.  **DSL-106: Schema Validation - `events`, `hrfs` Blocks**
    *   **Description:** Implement schema validation for the `events` block (column names) and the `hrfs` block (HRF definitions).
    *   **AC:**
        *   `events` block validates `onset_column`, `duration_column`, `block_column`.
        *   `hrfs` block (if present) validates structure of each HRF definition: `type` (enum), optional `parameters` (object), `derivatives` (for SPMCanonicalDerivs), `definition` (for CustomR).
    *   **Estimate:** 2 points

7.  **DSL-107: Schema Validation - `variables`, `transformations` Blocks**
    *   **Description:** Implement schema validation for the new `variables` and `transformations` blocks.
    *   **AC:**
        *   `variables` validates each entry has `bids_column`, `role` (enum).
        *   `transformations` validates each entry has `source_variable`, `ops` (array).
        *   `ops` items correctly validate simple strings (enum) or complex objects (type, params).
    *   **Estimate:** 3 points

8.  **DSL-108: Schema Validation - `confound_groups`, `terms` Blocks**
    *   **Description:** Implement schema validation for `confound_groups` and `terms`.
    *   **AC:**
        *   `confound_groups` validates structure (select\_by\_pattern/column).
        *   `terms` validates structure of each term: `type` (enum), conditional presence of `event_variables` / `selector_vars`+`mod_var` / `nuisance_source_variables`.
        *   Validates `hrf`, `subset`, `lag`.
        *   Validates `modulator_basis` structure (type, parameters, conditional required like `degree` for Polynomial).
    *   **Estimate:** 4 points

9.  **DSL-109: Schema Validation - `contrasts`, `models`, `default_model`, `validation_settings` Blocks**
    *   **Description:** Implement schema validation for `contrasts`, the `models` array, `default_model`, and `validation_settings`.
    *   **AC:**
        *   `contrasts` validates structure based on `type`.
        *   `models` validates as an array of model objects. Each model object: `name`, `baseline` (shorthand or structured `basis`, `intercept`, `include_confound_groups`), `terms` (array of strings), `contrasts` (array of strings).
        *   `default_model` is a string.
        *   `validation_settings` structure is validated.
    *   **Estimate:** 3 points

10. **DSL-110: Integrate Parsing, Defaults, and Schema Validation in `parse_and_validate_config(yaml_file)`**
    *   **Description:** Create the main entry function that calls DSL-101, DSL-104, and then the schema validation sequence (DSL-105 to DSL-109).
    *   **AC:**
        *   Function returns a validated R list (Internal Object Representation - IOR) if all checks pass.
        *   Stops with formatted errors via `ValidationErrors$stop_if_invalid()` if any schema check fails.
        *   Defaults are correctly applied before validation where appropriate.
    *   **Estimate:** 2 points

**Semantic Validation & BIDS Integration (EPI-DSL-SEMANTIC)**
*(These might spill into the next sprint or be done by a separate "BIDS integration" team/person)*

1.  **DSL-201: BIDS Project Loading and Basic Checks**
    *   **Description:** In `build_config_from_ior`, load the BIDS project using `bidser::bids_project(dataset$path)`. Check `dataset$path` existence.
    *   **AC:**
        *   BIDS project is loaded successfully.
        *   Error if `dataset$path` doesn't exist or isn't a valid BIDS dir.
        *   `ValidationErrors` updated for failures.
    *   **Estimate:** 1 point

2.  **DSL-202: Semantic Validation - Subject, Task, Run Availability**
    *   **Description:** Validate `dataset.subjects.include/exclude`, `dataset.tasks`, `dataset.runs` against actual content from `bidser::participants()`, `bidser::tasks()`, `bidser::runs()`.
    *   **AC:**
        *   Correctly identifies selected subjects, tasks, runs.
        *   Errors/warns (per `validation_settings.bids_content_checks`) if specified items are not found.
    *   **Estimate:** 2 points

3.  **DSL-203: Semantic Validation - Event Column Existence**
    *   **Description:** For a representative subject/task/run, load events.tsv and verify that `events.onset_column`, `duration_column`, `block_column` exist.
    *   **AC:**
        *   Correctly checks for essential event column existence.
        *   Errors/warns if columns are missing.
    *   **Estimate:** 2 points

4.  **DSL-204: Semantic Validation - `variables` to BIDS Column Mapping**
    *   **Description:** For representative data, verify that all `variables.*.bids_column` exist in either the events.tsv or available confounds.tsv columns for the relevant files.
    *   **AC:**
        *   Correctly checks mapping against representative BIDS files.
        *   Errors/warns for unmappable `bids_column` entries.
    *   **Estimate:** 3 points

5.  **DSL-205: Semantic Validation - Confound Group Resolution & Baseline Confounds**
    *   **Description:** Resolve `confound_groups.*.select_by_pattern` and `select_by_bids_column` against available confounds from a representative confounds.tsv. Validate that `models.*.baseline.include_confound_groups` refer to defined groups.
    *   **AC:**
        *   Correctly expands confound group patterns.
        *   Errors/warns if specified confound groups or their contents are not found.
    *   **Estimate:** 3 points

6.  **DSL-206: Semantic Validation - Cross-References (HRFs, Terms, Contrasts in Models)**
    *   **Description:** Implement the cross-reference checks:
        *   `terms.*.hrf` must exist in `hrfs` (or be "canonical").
        *   `models.*.terms` must exist in global `terms`.
        *   `models.*.contrasts` must exist in global `contrasts`.
        *   `default_model` must exist in `models.*.name`.
    *   **AC:**
        *   Errors/warns (per `validation_settings.cross_references`) for broken references.
    *   **Estimate:** 2 points

7.  **DSL-207: Semantic Validation - Role and Type Compatibility (Initial Pass)**
    *   **Description:** Implement basic checks for role/type compatibility (e.g., `modulator_basis` applied to a `Numeric` variable; `ParametricModulation` terms have a `Numeric` modulator).
    *   **AC:**
        *   Basic compatibility checks implemented.
        *   Errors/warns for obvious incompatibilities.
    *   **Estimate:** 2 points (can be expanded later)

**`fmri_config` Object Construction (EPI-DSL-BUILD)**

1.  **DSL-301: Define and Construct `fmri_config` S3 Object**
    *   **Description:** Define the structure of the `fmri_config` object. Create `build_config_from_ior(validated_ior)` function that populates this object after all schema and semantic validations pass.
    *   **AC:**
        *   `fmri_config` object contains `spec` (the validated IOR), `project` (bidser object), resolved `subjects`, `tasks`, `runs`, processed `events_info`, `confounds_info` (with resolved patterns), inferred/validated `variable_roles`, and `validated = TRUE` flag.
    *   **Estimate:** 2 points

2.  **DSL-302: Main User Function: `load_fmri_config(yaml_file)`**
    *   **Description:** Implement the top-level user function that orchestrates `parse_and_validate_config` and `build_config_from_ior`.
    *   **AC:**
        *   Takes YAML file path as input.
        *   Returns a fully validated `fmri_config` object or stops with errors.
    *   **Estimate:** 1 point

---

### **Stretch Goals (If time permits, or for next sprint):**

*   **DSL-401: Implement R code generation for `CustomR` HRFs.**
    *   Safely evaluate or source the `definition` string.
*   **DSL-402: Implement parser for shorthand `basis` strings (e.g., "BSpline(3)").**
*   **DSL-403: Advanced `transformations.ops` (e.g., `Clip`, `RecodeLevels`).**
*   **DSL-404: Begin `build_fmri_model_from_config(config, subject_id)` function.**
    *   Focus on loading subject-specific data using resolved paths from `fmri_config`.
    *   Setting up `sampling_frame`.

---

**Definition of Done for the Sprint:**

*   All tickets from EPI-DSL-CORE, EPI-DSL-SEMANTIC, and EPI-DSL-BUILD are completed and merged.
*   `load_fmri_config("my_analysis.yaml")` successfully parses, validates (schema and basic semantics), applies defaults, and returns a usable `fmri_config` object for a representative set of valid v2.4 YAML files.
*   Comprehensive unit tests cover parsing, default application, and validation logic for all major DSL sections.
*   Clear error messages are produced for invalid DSL files.
*   Basic documentation for the DSL v2.4 structure is drafted.

Addendum: working with bidser package:

Okay, here's a concise guide for the engineer working on the BIDS integration parts of Sprint 1 (DSL-201 to DSL-205), focusing on how to leverage the `bidser` package.

**Key `bidser` Functions for DSL Sprint 1 (Semantic Validation & BIDS Integration):**

The engineer will primarily interact with the `fmri_config$project` object, which is a `bids_project` object created by `bidser::bids_project()`.

1.  **Loading the BIDS Project (DSL-201):**
    *   **`bidser::bids_project(path_to_bids_root, fmriprep = FALSE)`:** This is the first step. The `path_to_bids_root` comes from `fmri_config$spec$dataset$path`. `fmriprep = FALSE` is specified in `fmridsl_builder.R::build_config_from_ior`, meaning we are primarily interested in the raw BIDS structure for now, though `bidser` can also handle derivatives. The `fmri_config$project` slot will hold this object.
    *   *Engineer's Task:* Ensure this is called correctly and the resulting object is stored. Check for basic path existence before calling.

2.  **Validating Subject, Task, Run Availability (DSL-202):**
    *   **`bidser::participants(bids_project_object)`:** Returns a character vector of participant IDs (e.g., `"01"`, `"02"`) found in `participants.tsv`. These are *without* the `"sub-"` prefix.
        *   *Engineer's Task:* Compare the `fmri_config$spec$dataset$subjects$include` / `exclude` lists (which should also be prefix-stripped for comparison) against this vector.
    *   **`bidser::tasks(bids_project_object)`:** Returns a character vector of unique task names (e.g., `"rest"`, `"nback"`) found across the dataset. These are *without* the `"task-"` prefix.
        *   *Engineer's Task:* Compare `fmri_config$spec$dataset$tasks` against this.
    *   **`bidser::runs(bids_project_object, subid, task, session)`:** (Conceptual - `bidser` might not have a direct `runs()` function that lists all available run numbers for a given subject/task/session. `search_files` might be used, or a more specific helper within `bidser`.)
        *   *Engineer's Task:* To validate specified runs, the engineer might need to use `bidser::search_files()` with appropriate filters (e.g., for `_bold.nii.gz` files for a given subject/task) and then parse the `run-` entity from the returned filenames or their encoded attributes.
        *   Alternatively, if `fmri_config$spec$dataset$runs` is provided, it can be a primary filter. The check is more about whether *any* data can be found for these runs, rather than an exhaustive list of all possible runs.

3.  **Validating Event Column Existence (DSL-203):**
    *   **`bidser::read_events(bids_project_object, subid, task, run, session)`:** This function reads the `events.tsv` file(s) matching the criteria and returns a nested tibble. The actual event data is in the `data` list-column.
        *   *Engineer's Task:* Call this for a *representative* subject/task/run (e.g., the first valid combination found). Access the first data frame (e.g., `read_events_output$data[[1]]`) and check if `fmri_config$spec$events$onset_column`, `duration_column`, and `block_column` exist as column names in that data frame.

4.  **Validating `variables` to BIDS Column Mapping (DSL-204):**
    *   This involves checking two sources: `events.tsv` and `confounds.tsv`.
    *   **For event-derived variables:** Use the representative `events.tsv` data frame loaded in DSL-203. For each entry in `fmri_config$variables` where the `bids_column` is expected to come from an event file, check if that `bids_column` exists in the representative `events_df`.
    *   **For confound-derived variables (NuisanceSource):**
        *   **`bidser::confound_files(bids_project_object, subid, task, run, session)`:** Finds paths to potential confound files.
        *   **`bidser::read_confounds(bids_project_object, ...)`:** Reads these files. This function also handles selection of specific confound variables (`cvars` argument).
        *   *Engineer's Task:* For a representative subject/task/run, load the confounds. For each entry in `fmri_config$variables` with `role: NuisanceSource`, check if its `bids_column` exists in the columns of the loaded confound data.

5.  **Validating Confound Group Resolution & Baseline Confounds (DSL-205):**
    *   Use `bidser::read_confounds()` for a representative subject/task/run to get the *full list of available confound column names* from the BIDS dataset.
    *   *Engineer's Task:*
        *   For each group in `fmri_config$confound_groups`:
            *   If `select_by_pattern` is used, apply `grep()` for each pattern against the list of available confound names from BIDS. Collect all matches.
            *   If `select_by_bids_column` is used, check if these explicit column names exist in the available BIDS confounds.
        *   For each model in `fmri_config$models`:
            *   For each group name listed in `model$baseline$include_confound_groups`:
                *   Check if the group name is a defined key in `fmri_config$confound_groups`.
                *   Verify that the resolution of this group (from the step above) yields at least one actual confound column from BIDS.

**General `bidser` Usage Notes for the Engineer:**

*   **Entity Naming:** `bidser` functions often return BIDS entity values *without* their prefixes (e.g., `participants()` returns `"01"`, not `"sub-01"`). The DSL stores user input which might or might not have prefixes. The engineer needs to be consistent when comparing, likely by stripping prefixes from DSL values or adding them to `bidser` outputs before comparison. The `generate_bids_filename` and `generate_bids_path` helpers in `mock_bids.R` (which might be internalized or moved to `bidser`) handle prefixing correctly.
*   **File Searching with `search_files()`:** This is a versatile function.
    *   `search_files(proj, regex = "pattern", subid = "01", task = "rest", full_path = TRUE, fmriprep = TRUE/FALSE, ...)`
    *   It matches against BIDS entities stored as attributes on file nodes in `bidser`'s internal tree. The keys for these attributes are usually the short BIDS keys (e.g., `sub`, `ses`, `task`, `run`, `desc`, `space`, `suffix`, `kind`).
    *   The `fmriprep` argument in `search_files` can be used to target raw BIDS (`fmriprep = FALSE`) or derivatives (`fmriprep = TRUE`).
*   **Encoding Filenames with `bidser::encode(filename_string)`:** This utility can parse a BIDS filename string and return a list of its constituent entities. This can be useful for robustly extracting metadata if needed, though direct attribute access from `search_files` results is often simpler.
*   **Error Handling:** Wrap calls to `bidser` functions in `tryCatch` or check return values (e.g., for `NULL` or zero-length vectors/data frames) to handle missing data or BIDS inconsistencies gracefully and report them via the `ValidationErrors` object.

By using these `bidser` functions, the engineer can efficiently query the BIDS dataset structure and content to perform the necessary semantic validations required by the DSL implementation.