Okay, this is a substantial and well-structured R package for fMRI regression modeling. Let's break down its architecture and philosophy, focusing on the core modeling aspects.

## fmrireg: Architectural & Philosophical Overview

**Core Philosophy:**

1.  **Modularity & Composability:** The library is built on distinct, interacting components: HRF definitions, event specifications, baseline modeling, and finally, the full GLM. This allows for flexibility and reuse.
2.  **Formula-Driven Specification (for Events):** The primary way to define event-related regressors is through R's formula interface (e.g., `onset ~ hrf(condition) + hrf(modulator)`), which is intuitive for R users familiar with `lm()` or `glm()`.
3.  **Explicit HRF Handling:** Hemodynamic Response Functions (HRFs) are first-class citizens. Users can choose from predefined HRFs, create custom ones, or use basis sets. HRFs can be decorated (lagged, blocked, normalized).
4.  **Structured Event Representation:** Events are internally represented in a structured way (`event_term`, `event` objects) that separates the raw event attributes (onsets, durations, factor levels, modulator values) from their convolution with HRFs.
5.  **Clear Naming Conventions:** There's a significant emphasis on generating unique, informative, and predictable column names for the final design matrix. This is crucial for interpreting model outputs and defining contrasts.
6.  **Efficiency for Core Operations:** Computationally intensive parts, like regressor evaluation (convolution), are offloaded to C++ (`evaluate_regressor_cpp`).
7.  **Separation of Specification and Realization:**
    *   `hrfspec` objects define *how* a term should be constructed.
    *   `event_term` objects are the *realization* of that specification with actual data.
    *   The `event_model` object orchestrates this process and assembles the final design matrix.
8.  **Configuration via DSL (fmriDSL):** For larger projects or standardized pipelines, a YAML-based Domain Specific Language (DSL) is provided (`fmridsl.R`, `fmridsl_builder.R`) to define entire analysis configurations, which then build `fmri_model` objects.

---

**Key Architectural Components & Workflow:**

Here's a conceptual flow from user input to a design matrix:

```mermaid
graph TD
    A[User Input: Formula or List of hrfspecs + Data] --> B(parse_event_model);
    B --> C{hrfspec Objects};
    C --> D(realise_event_terms);
    D -- calls --> E{construct S3 generic};
    E -- dispatches to --> F(construct.hrfspec);
    F -- calls --> G(construct_event_term);
    G -- creates --> H[event_term Objects];
    H -- contains --> I[event Objects];
    H --> J(build_event_model_design_matrix);
    J -- calls --> K(convolve.event_term);
    K -- uses --> L[HRF Object from hrfspec];
    K -- uses --> M[sampling_frame];
    K -- may use --> N(evaluate.Reg internally);
    J --> O[Final Event Design Matrix (DM)];

    P[User Input: Baseline Spec + Sampling Frame] --> Q(baseline_model);
    Q --> R[Baseline DM];

    O & R --> S(fmri_model);
    S --> T[Full Design Matrix];
```

1.  **HRF System (`hrf.R`, `hrf-*.R`, `hrf_decorators.R`):**
    *   **`HRF` S3 Class:** The fundamental unit. It's a function `f(t)` with attributes like `name`, `nbasis`, `span`, `params`.
    *   **Constructors:**
        *   `as_hrf()`: Core internal constructor.
        *   Predefined objects: `HRF_SPMG1`, `HRF_GAMMA`, etc.
        *   Generator functions: `hrfspline_generator()`, `hrf_tent_generator()` create `HRF` objects based on parameters (like `nbasis`).
    *   **Access:** `getHRF(name, ...)` is the central way to get an HRF instance, applying parameters like `nbasis`, `span`.
    *   **Manipulation:**
        *   Decorators: `lag_hrf()`, `block_hrf()`, `normalise_hrf()` take an `HRF` and return a new, modified `HRF`.
        *   Combining: `bind_basis()` (or `hrf_set()`) takes multiple `HRF` objects and creates a single multi-basis `HRF` object. Its function will `cbind` the outputs of the component HRFs.
    *   **Evaluation:** `evaluate.HRF(hrf_obj, grid, duration, ...)` calculates the HRF's value(s) over a time `grid`, handling event `duration` and `amplitude`.

2.  **Low-Level Regressor Object (`regressor.R`, `reg-constructor.R`, `reg-methods.R`):**
    *   **`Reg` S3 Class:** Represents a sequence of events convolved with a *single* HRF object.
    *   Stores `onsets`, `duration` (vector), `amplitude` (vector), the `hrf` object itself, and `span`.
    *   Created by `Reg()` (internal) or its public alias `regressor()`.
    *   **Evaluation:** `evaluate.Reg(reg_obj, grid, precision, method)` is the workhorse for generating the time course.
        *   `prep_reg_inputs()`: Pre-processes onsets, durations, amplitudes, and memoizes the HRF evaluation on a fine grid.
        *   Dispatches to C++ engines (`eval_fft`, `eval_conv`) or R implementations (`eval_Rconv`, `eval_loop`).
    *   This object seems to be used internally by `convolve.event_term` for each specific condition derived from an `event_term`.

3.  **Event Specification & Model Generation (The Main Pipeline):**
    This is the core of generating the event-related part of the design matrix.

    *   **User Input:**
        *   A `formula_or_list`:
            *   Formula: e.g., `onsets ~ hrf(condition) + hrf(Modulator, basis="spmg3")`
            *   List: A list of `hrfspec` objects.
        *   `data`: A `data.frame` with event variables.
        *   `block`: Formula or vector for block/run identifiers.
        *   `sampling_frame`: Defines scan timing (TR, number of scans per block). Crucial for convolution.

    *   **Stage 1: Parsing (`event_model_helpers.R::parse_event_model`)**
        *   Normalizes `formula_or_list` input.
        *   If formula, `parse_event_formula` is called:
            *   LHS (e.g., `onsets`) is evaluated in `data` to get actual onset times.
            *   RHS is traversed by `find_and_eval_hrf_calls` to find calls to `hrf()`, `trialwise()`, `afni_hrf()`, `covariate()`.
            *   Each such call is evaluated in an environment derived from the formula's environment and `data`. This evaluation produces an `hrfspec` (or `afni_hrfspec`, `covariatespec`) object.
        *   If list, it's assumed to be a list of `hrfspec` objects already. `data$onset` is used for onsets.
        *   Validates and canonicalizes `blockids` and `durations`.
        *   **Output:** A list containing `spec_tbl` (a tibble of `hrfspec` objects), `onsets`, `durations`, `blockids`, `data`, `formula_env`, and `interface` type.

    *   **Stage 2: Realizing Event Terms (`event_model_helpers.R::realise_event_terms`)**
        *   Iterates through each `hrfspec` in `parsed_spec$spec_tbl`.
        *   For each `hrfspec`:
            *   A unique `term_uid` (e.g., "t01") is generated.
            *   A `term_tag` is generated by `naming-utils.R::make_term_tag()`. This tag is crucial for column naming. It's derived from `hrfspec$id` (if provided in `hrf(..., id="my_id")`) or auto-generated from variable names (e.g., `hrf(cond)` -> `condition`). Handles clashes by appending `#1`, `#2`.
            *   Calls the S3 generic `construct(hrfspec, model_spec)`. This dispatches to methods like `construct.hrfspec`.
        *   **`construct.hrfspec` (from `hrf-formula.R`):**
            *   Calls `construct_event_term(hrfspec, model_spec)`.
        *   **`construct_event_term` (from `event_model_helpers.R`):**
            *   This is a key function. It takes an `hrfspec` and the `model_spec` (which includes `data`, `onsets`, etc.).
            *   It evaluates the variables/expressions listed in `hrfspec$vars` (which are quosures) within the context of `data` and the formula environment.
            *   For each variable/expression, it creates an `event` object using the public wrappers (`event_factor`, `event_variable`, `event_basis`, etc., from `event-classes.R`).
            *   `event` objects (class `c("event", "event_seq")`) store the actual onsets, durations, block IDs, and the *payload* (`$value` as a numeric matrix, e.g., indicator columns for factors, numeric values for modulators, basis matrix for `ParametricBasis`). Metadata (like factor levels or the `ParametricBasis` object) is stored in `$meta`.
            *   These `event` objects are collected into a list.
            *   It then calls the public `event_term()` constructor with this list of evaluated event variables and shared timing info.
        *   **`event_term` object (from `event_vector.R`):**
            *   Represents a single term in the event model (e.g., `hrf(condition, modulator)`).
            *   Contains a list of the constituent `event` objects (one for `condition`, one for `modulator`).
            *   Also stores shared `onsets`, `blockids`, `durations`, and an `event_table` (a descriptive table of event combinations).
            *   The `hrfspec` is attached as an attribute to the `event_term`. The `term_tag` and `uid` are also attached.
        *   **Output:** A named list of `event_term` objects, where names are the `term_tag`s.

    *   **Stage 3: Building the Design Matrix (`event_model_helpers.R::build_event_model_design_matrix`)**
        *   Takes the list of `event_term` objects.
        *   Iterates through each `event_term`:
            *   Calls `convolve.event_term(event_term_obj, hrf = hrfspec$hrf, sampling_frame, ...)`.
        *   **`convolve.event_term` (from `event_vector.R`):**
            *   This is where the actual convolution happens and where column names are finalized.
            *   It retrieves the `term_tag` attribute from the `event_term`.
            *   It gets the unconvolved design matrix for the term using `design_matrix.event_term()`. The columns of this intermediate matrix represent the different conditions/levels within the term (e.g., "Cond.A_Modulator", "Cond.B_Modulator").
            *   For each column of this intermediate design matrix (representing a specific condition), it creates a temporary `Reg` object using the onsets, the column values as amplitudes, and the HRF from the `hrfspec`.
            *   It evaluates these `Reg` objects over the `sampling_frame` (likely using `evaluate.Reg` which can dispatch to C++).
            *   The resulting convolved time courses are `cbind`ed.
            *   **Crucially, it calls `naming-utils.R::make_column_names(term_tag, base_condition_names, nbasis)` to generate the final, globally unique, and structured column names for this term's part of the design matrix.**
        *   The convolved matrices from all terms are `cbind`ed to form the full event design matrix.
        *   Attributes `term_spans` (indicating which columns belong to which term) and `col_indices` are attached to the final design matrix.

    *   **Final `event_model` Object (from `event_model.R::event_model()`):**
        *   A list containing:
            *   `terms`: The list of realized `event_term` objects.
            *   `design_matrix`: The final convolved design matrix (a tibble).
            *   `blockids`, `sampling_frame`.
            *   `contrasts`: (Note: contrasts are now primarily defined within `hrfspec`s and applied later).
            *   `model_spec`: Information about the original call.

4.  **Naming Scheme (`naming-utils.R` and used in `convolve.event_term`):**
    *   **Grammar:** `term_tag` + `_` + `condition_tag` + `_b##` (optional basis suffix).
    *   `term_tag`:
        *   From `hrf(..., id = "my_id")` or `hrf(..., name = "my_name")`.
        *   Auto-generated from variable names in `hrf(...)` if no `id`/`name` (e.g., `hrf(condition)` -> `condition`; `hrf(Poly(RT,2))` -> `Poly_RT`).
        *   `make_term_tag()` ensures uniqueness (e.g., `condition`, `condition#1`).
        *   Sanitized (dots to underscores) by `sanitize()`.
    *   `condition_tag`:
        *   Derived from the `event_term`'s internal structure. This happens via `conditions.event_term()`.
        *   For a simple factor `hrf(condition)` where `condition` has levels "A", "B": `condition.A`, `condition.B`.
        *   For an interaction `hrf(factorA, factorB)`: `factorA.L1_factorB.X1`, etc. (uses `_` as separator).
        *   For a continuous modulator `hrf(Poly(RT,2))`: `01`, `02` (from `levels.ParametricBasis`).
        *   For `hrf(Ident(X,Y))`: `X`, `Y`.
    *   `_b##`:
        *   Added by `add_basis()` if `nbasis(hrf) > 1`.
        *   Uses `basis_suffix()` which calls `zeropad()` (e.g., `_b01`, `_b02`).
    *   `make_column_names(term_tag, condition_tags, nbasis)`: Assembles the final names. If `term_tag` is `NULL` (e.g., for `hrf(Ident(X,Y))` without an `id`), it just uses `condition_tags` (potentially with basis suffix), resulting in direct column names like `X_b01`, `Y_b01`.
    *   Global uniqueness of final column names is ensured in `build_event_model_design_matrix` by a final `make.names(..., unique=TRUE)` pass if any initial clashes remain after the structured naming.

5.  **Baseline Model (`baseline_model.R`):**
    *   Models drift, block intercepts, and nuisance regressors.
    *   `baseline_model()` constructor takes basis type (`constant`, `poly`, `bs`, `ns`), degree, `sampling_frame`, intercept type (`runwise`, `global`, `none`), and an optional `nuisance_list`.
    *   Internally creates `baseline_term` objects for drift, block, and nuisance components.
    *   `design_matrix.baseline_model()` combines these.

6.  **Full fMRI Model (`fmri_model.R`):**
    *   `fmri_model(event_model, baseline_model)`: Simply combines the two.
    *   `design_matrix.fmri_model()`: `cbind`s the design matrices from the event and baseline models.
    *   `terms()`: Returns combined terms.
    *   Provides methods for `contrast_weights`, `conditions`, plotting, etc.

7.  **Contrasts (`contrast.R`):**
    *   Defined using functions like `contrast()`, `pair_contrast()`, `column_contrast()`, `poly_contrast()`.
    *   These create `contrast_spec` objects.
    *   `contrast_weights.event_model()` and `contrast_weights.fmri_model()`:
        *   Iterate through terms.
        *   For each term, retrieve its defined contrasts (from `hrfspec`).
        *   Call `contrast_weights()` on the specific `contrast_spec` and the `event_term`.
        *   `contrast_weights.pair_contrast_spec()` (etc.): Evaluates the contrast logic against the (expanded) conditions of the term, producing a named vector of weights.
        *   These term-local weights are then mapped to the full design matrix using `col_indices` from the `event_model`.

8.  **Model Fitting (`fmrilm.R`, `fmrirlm.R`):**
    *   `fmri_lm()` / `fmri_rlm()`: High-level functions that take a formula, dataset, etc.
    *   They internally call `create_fmri_model()` to build the `fmri_model` object.
    *   Then, `fmri_lm_fit()` is called, which handles different strategies (`runwise`, `chunkwise`) and can use `lm.fit` or robust methods.
    *   The "fast path" (`use_fast_path=TRUE`) uses pre-projectors (`.fast_preproject`, `.fast_lm_matrix`) for potentially faster OLS estimation.

---

**Diagram of Key Objects in Event Modeling:**

```mermaid
graph LR
    subgraph UserInput
        direction LR
        FMLA[Formula: onset ~ hrf(V1, V2, id="T1")];
        DATA[data.frame: V1, V2, onset, block];
        SFRAME[sampling_frame];
    end

    subgraph ParsingAndRealization
        direction TB
        HRFSPEC[hrfspec for "T1" \n vars: (V1, V2)\n id: "T1"\n hrf: HRF_obj];
        ET_CONSTRUCT[construct_event_term];
        EV1[event for V1 \n (e.g., factor levels)];
        EV2[event for V2 \n (e.g., numeric values)];
        ET[event_term "T1"\n Events: (EV1, EV2)\n Onsets, Durations...\n Attr: hrfspec, term_tag="T1"];
    end

    subgraph ConvolutionAndNaming
        direction TB
        CONV[convolve.event_term];
        DM_INTERMEDIATE[Intermediate DM for "T1" \n Cols: V1.L1_V2, V1.L2_V2...];
        REG[Reg Objects \n (for each col of DM_INTERMEDIATE)];
        EVAL_REG[evaluate.Reg \n (uses C++)];
        MAKE_CNAMES[make_column_names \n term_tag="T1"\n cond_tags from DM_INTERMEDIATE cols \n nbasis from HRF_obj];
        CONVOLVED_TERM_DM[Convolved DM for "T1" \n Cols: T1_V1.L1_V2_b01...];
    end
    
    subgraph FinalAssembly
        direction LR
        EM[event_model Object];
        BM[baseline_model Object];
        FM[fmri_model Object];
    end

    FMLA & DATA & SFRAME -->|"event_model()" calls `parse_event_model`| HRFSPEC;
    HRFSPEC -->|"realise_event_terms" calls `construct` -> `construct_event_term`| ET_CONSTRUCT;
    ET_CONSTRUCT --> EV1;
    ET_CONSTRUCT --> EV2;
    ET_CONSTRUCT --> ET;
    
    ET -->|"build_event_model_design_matrix" calls| CONV;
    CONV --> DM_INTERMEDIATE;
    DM_INTERMEDIATE --> REG;
    REG --> EVAL_REG;
    CONV --> MAKE_CNAMES;
    MAKE_CNAMES --> CONVOLVED_TERM_DM;
    EVAL_REG --> CONVOLVED_TERM_DM;
    
    CONVOLVED_TERM_DM -->|cbind with other terms' DMs| EM;
    UserInput --> BM; %% Simplified path for baseline
    EM & BM --> FM;

```

This overview highlights the journey from specifying an event model using formulas to the creation of a convolved design matrix with structured column names. The modularity allows different parts (HRFs, event definitions, baseline) to be developed and tested somewhat independently, while the pipeline ensures they are combined correctly. The DSL provides a higher-level abstraction for users who prefer configuration files over programmatic R scripting.Okay, this is a sophisticated R package for fMRI regression modeling. Here's a "lay of the land" breakdown focusing on its philosophy, architecture, and key components, particularly around event modeling, HRFs, and naming.

## fmrireg: Philosophy and Architectural Overview

**Core Philosophy:**

1.  **Modularity and Composability:** The library is designed with distinct, yet interacting, components. You define Hemodynamic Response Functions (HRFs), specify experimental events, model baseline signal components (drift, nuisance), and then combine these into a full regression model. This promotes flexibility and reusability.
2.  **Formula-Driven Event Specification:** A key design choice is the use of R's formula interface (e.g., `onset ~ hrf(condition) + hrf(modulator)`) for defining the event-related part of the model. This is intuitive for R users familiar with `lm()` or `glm()`.
3.  **Explicit and Flexible HRF Handling:** HRFs are first-class citizens. The library provides:
    *   Predefined canonical HRFs (SPM-like, Gaussian, Gamma).
    *   Support for multi-basis HRFs (e.g., SPMG1 + temporal derivative = SPMG2).
    *   Generators for flexible basis sets (B-splines, Fourier, Tent).
    *   Decorators to modify HRFs (e.g., `lag_hrf`, `block_hrf`, `normalise_hrf`).
    *   A clear mechanism (`as_hrf`, `bind_basis`) for creating and combining HRF objects.
4.  **Structured Event Representation:** Raw event information (onsets, durations, conditions, modulators) is processed into structured `event` and `event_term` objects. This separates the *specification* of events from their *realization* as regressors after convolution.
5.  **Systematic Design Matrix Naming:** A significant effort is dedicated to generating unique, informative, and predictable column names for the final design matrix. This is crucial for interpreting model outputs and defining contrasts correctly and robustly.
6.  **Efficiency for Core Computations:** Computationally intensive tasks, particularly the convolution of event sequences with HRFs (`evaluate.Reg`), are implemented in C++ for speed.
7.  **Separation of Concerns:**
    *   The `event_model` primarily deals with generating regressors related to experimental events.
    *   The `baseline_model` handles regressors for signal drift, run-specific intercepts, and other nuisance covariates.
    *   The `fmri_model` combines these two into a complete model specification.
8.  **Configuration via DSL (Domain Specific Language):** For more complex or standardized pipelines, `fmrireg` offers a YAML-based DSL (`fmridsl.R`, `fmridsl_builder.R`). This allows users to define entire analysis configurations, which are then parsed and used to build `fmri_model` objects programmatically.

---

**Key Architectural Components & Workflow for fMRI Modeling:**

Let's trace the typical workflow for generating an fMRI design matrix:

**I. Hemodynamic Response Function (HRF) System:**
   (Primarily in `hrf.R`, `hrf-functions.R`, `hrf_decorators.R`)

*   **Core Object:** `HRF` (S3 class). An HRF object is essentially a function `f(t)` that returns the HRF's amplitude at time `t`, augmented with attributes like `name`, `nbasis` (number of basis functions, e.g., 1 for SPMG1, 2 for SPMG2), `span` (duration), and `params`.
*   **Creation:**
    *   `as_hrf(function, ...)`: The fundamental constructor.
    *   Predefined objects: `HRF_SPMG1`, `HRF_GAMMA`, etc.
    *   Generator functions: `hrfspline_generator(nbasis=5)`, `hrf_tent_generator(nbasis=7)` return `HRF` objects.
    *   `getHRF(name, nbasis, ...)`: A central dispatcher to retrieve or generate an HRF instance by name, applying parameters.
*   **Manipulation:**
    *   `lag_hrf(hrf, lag)`: Shifts an HRF in time.
    *   `block_hrf(hrf, width)`: Convolves an HRF with a boxcar of given `width` to model sustained events.
    *   `normalise_hrf(hrf)`: Scales an HRF to have a peak amplitude of 1.
    *   `bind_basis(hrf1, hrf2, ...)` or `hrf_set(...)`: Combines multiple `HRF` objects into a single multi-basis `HRF`. The resulting function will `cbind` the outputs of its constituent HRFs.
*   **Evaluation:** `evaluate.HRF(hrf_obj, grid, duration, amplitude, ...)`: Calculates the HRF's value(s) over a time `grid`, considering event `duration` (for block HRFs) and `amplitude`.

**II. Event Specification and Model Generation (The Core Pipeline):**
    (Primarily in `event_model.R`, `event_model_helpers.R`, `hrf-formula.R`, `event-classes.R`, `event_vector.R`)

This is the most intricate part, responsible for translating user-specified experimental designs into a convolved design matrix.

```mermaid
graph TD
    A[User: Formula & Data \n e.g., onset ~ hrf(cond) + hrf(RT, id="RespTime")] --> B{1. Parse Inputs \n (parse_event_model)};
    B -- creates list of --> C[hrfspec Objects \n (Specifications for each term)];
    C --> D{2. Realise Event Terms \n (realise_event_terms)};
    D -- for each hrfspec, calls --> E[construct (S3 generic)];
    E -- dispatches to --> F[construct.hrfspec];
    F -- calls --> G[construct_event_term \n (Evaluates vars, creates 'event' objects)];
    G -- produces --> H[event_term Object \n (Contains 'event' objects, onsets, etc. \n Attributes: hrfspec, term_tag)];
    H --> I{3. Build Design Matrix \n (build_event_model_design_matrix)};
    I -- for each event_term, calls --> J[convolve.event_term];
    J -- uses HRF from event_term's hrfspec --> K[HRF Object];
    J -- uses --> L[sampling_frame];
    J -- internally may use --> M[Reg Object & evaluate.Reg];
    J -- finalizes col names via --> N[make_column_names];
    J -- produces --> O[Convolved Matrix for Term];
    O --> P[Final Event Design Matrix \n (Attributes: term_spans, col_indices)];
    P --> Q[event_model Object];
```

1.  **User Input:**
    *   `formula_or_list`: Typically a formula like `onsets ~ hrf(condition) + hrf(modulator, basis="spmg3")`. The LHS specifies the onset variable. The RHS uses special functions:
        *   `hrf(...)`: Defines an event-related term. Arguments specify variables from `data`, the `basis` HRF, an optional `id` or `name` for the term, `contrasts`, etc.
        *   `trialwise(...)`: A wrapper around `hrf()` for trial-specific regressors.
        *   `covariate(...)`: For regressors that are *not* convolved with an HRF (e.g., motion parameters).
        *   `afni_hrf(...)` / `afni_trialwise(...)`: For AFNI-specific HRF models.
    *   `data`: A `data.frame` containing columns referenced in the formula.
    *   `block`: A formula (e.g., `~ run`) or vector defining the run/block structure.
    *   `sampling_frame`: An object defining scan timing (TR, scans per block). Essential for convolution.
    *   `durations`: Event durations.

2.  **Stage 1: Parsing Inputs (`parse_event_model` in `event_model_helpers.R`)**
    *   Normalizes the `formula_or_list` input.
    *   If a formula is provided, `parse_event_formula` is called:
        *   It evaluates the LHS (e.g., `onsets`) in the `data` environment to get actual onset times.
        *   It traverses the RHS using `find_and_eval_hrf_calls` to identify calls to `hrf()`, `trialwise()`, etc.
        *   Each such call is *evaluated* in an environment that combines the formula's environment and the `data`. This evaluation creates an **`hrfspec` object** (defined in `hrf-formula.R`). An `hrfspec` is a list that captures the *specification* of a term: the variables involved (as quosures), the chosen HRF object, contrast definitions, subsetting expressions, and any explicit `id`.
    *   If a list of `hrfspec` objects is provided directly, this stage is simpler.
    *   It also validates and canonicalizes `blockids` and `durations`.
    *   **Output:** A list containing `spec_tbl` (a tibble of `hrfspec` objects), processed `onsets`, `durations`, `blockids`, `data`, `formula_env`, and the `interface` type.

3.  **Stage 2: Realizing Event Terms (`realise_event_terms` in `event_model_helpers.R`)**
    *   This stage iterates through each `hrfspec` from the `spec_tbl`.
    *   For each `hrfspec`:
        *   A unique `term_uid` (e.g., "t01", "t02") is generated.
        *   A **`term_tag`** is generated using `make_term_tag()` (from `naming-utils.R`). This tag is crucial for the final column naming. It's derived from `hrfspec$id` (if the user provided `hrf(..., id="my_id")`) or auto-generated from the variable names in the `hrfspec` (e.g., `hrf(condition)` might yield `condition`; `hrf(Poly(RT,2))` might yield `Poly_RT`). `make_term_tag` ensures uniqueness by appending `#1`, `#2`, etc., if clashes occur.
        *   The S3 generic `construct(hrfspec_object, model_spec)` is called. This dispatches to a method like `construct.hrfspec` (defined in `hrf-formula.R`).
    *   **`construct.hrfspec` method:**
        *   Its main job is to call `construct_event_term(hrfspec_object, model_spec)`.
    *   **`construct_event_term` (in `event_model_helpers.R`):**
        *   This is a key function that takes the *specification* (`hrfspec`) and the *data context* (`model_spec` which includes `data`, `onsets`, `blockids`, `durations`, `formula_env`).
        *   It evaluates the variables/expressions listed in `hrfspec$vars` (which are stored as quosures). For example, if `hrfspec$vars` contains `quo(condition)` and `quo(Poly(RT,2))`, it evaluates `condition` and `Poly(RT,2)` in the provided data environment.
        *   For each evaluated variable/expression, it creates an **`event` object**. These `event` objects (defined in `event-classes.R`, with class `c("event", "event_seq")`) are the fundamental building blocks. They store the actual onsets, durations, block IDs relevant to *that specific variable sequence*, and its payload (e.g., for a factor, `$value` is a matrix of integer codes; for a numeric modulator, `$value` is a matrix of its values; for a `ParametricBasis` like `Poly(RT,2)`, `$value` is the basis matrix). Metadata (like factor levels or the `ParametricBasis` object) is stored in `event$meta`. `event` objects are created by public wrappers like `event_factor()`, `event_variable()`, `event_basis()`.
        *   The list of these processed `event` objects, along with the shared `onsets`, `blockids`, and `durations` (subsetted if `hrfspec$subset` was used), are passed to the **`event_term()` constructor** (defined in `event_vector.R`).
    *   **`event_term` object:**
        *   Represents a single realized term from the model formula (e.g., `hrf(condition, modulator)` becomes one `event_term`).
        *   It contains a list of the constituent `event` objects (e.g., one `event` for `condition`, one `event` for `modulator`).
        *   It also stores the (potentially subsetted) `onsets`, `blockids`, `durations`, and an `event_table` (a descriptive table of all unique combinations of event levels/values within that term).
        *   The original `hrfspec` and the generated `term_tag` and `uid` are attached as attributes to this `event_term` object.
    *   **Output:** A named list of `event_term` objects. The names of this list are the unique `term_tag`s.

4.  **Stage 3: Building the Final Design Matrix (`build_event_model_design_matrix` in `event_model_helpers.R`)**
    *   Takes the named list of `event_term` objects and the `sampling_frame`.
    *   Iterates through each `event_term`:
        *   Retrieves the HRF object from the `hrfspec` attached to the `event_term`.
        *   Calls **`convolve.event_term(event_term_object, hrf, sampling_frame, ...)`** (defined in `event_vector.R`).
    *   **`convolve.event_term`:**
        *   This is where the actual convolution happens and where the **final column names are generated**.
        *   It retrieves the `term_tag` attribute from the `event_term`.
        *   It gets the *unconvolved* design matrix for the term using `design_matrix.event_term()`. The columns of this intermediate matrix represent the different specific conditions or basis functions within the term (e.g., for `hrf(condition, Poly(RT,2))`, columns might represent `condition.A_Poly.RT.1`, `condition.A_Poly.RT.2`, `condition.B_Poly.RT.1`, etc.). The names of these columns are the `condition_tag`s (or parts of them).
        *   For each column of this intermediate "unconvolved design matrix" (which effectively represents a single "stick function" to be convolved):
            *   It creates a temporary `Reg` object (from `reg-constructor.R`). This `Reg` object gets the overall event `onsets`, the HRF (from the `hrfspec`), and uses the values from the current column of the intermediate design matrix as its `amplitude`.
            *   It evaluates this `Reg` object over the `sampling_frame` using `evaluate.Reg()`. This function can dispatch to efficient C++ implementations (`eval_fft`, `eval_conv`).
        *   The resulting convolved time courses (one for each column of the intermediate design matrix, potentially expanded further if the HRF itself is multi-basis) are `cbind`ed.
        *   **Crucially, it calls `make_column_names(term_tag, base_condition_names, nbasis(hrf))` (from `naming-utils.R`) to generate the final, globally unique, and structured column names for this term's part of the overall design matrix.**
    *   The convolved matrices from all terms are `cbind`ed to form the full event-related design matrix.
    *   Attributes like `term_spans` (indicating which columns belong to which original term) and `col_indices` are attached to this final design matrix.

5.  **Final `event_model` Object Construction (by `event_model()` in `event_model.R`)**
    *   This is a list wrapper around the results of the pipeline.
    *   Contains:
        *   `terms`: The list of realized `event_term` objects.
        *   `design_matrix`: The final convolved design matrix (a tibble) with all attributes.
        *   `blockids`, `sampling_frame`.
        *   `contrasts`: (Note: In the newer design, contrasts are primarily defined within `hrfspec`s and applied later during statistical analysis rather than being a direct part of the `event_model` object's main structure, though `contrast_weights.event_model` can extract them).
        *   `model_spec`: Stores some information about the original call for reproducibility/inspection.

**III. Naming Scheme (`naming-utils.R`):**

This is a critical part for interpretability and programmatic access to design matrix columns.

*   **Goal:** Create unique, informative, and predictable column names.
*   **Structure:** `term_tag` + `_` + `condition_tag` + (optional)`_b##`
    *   **`term_tag`**:
        *   User-defined via `hrf(..., id = "my_term_id")`.
        *   Auto-generated by `make_term_tag()` based on variable names in `hrf(...)` (e.g., `hrf(condition)` -> `condition`; `hrf(RT, accuracy)` -> `RT_accuracy`; `hrf(Poly(RT,2))` -> `Poly_RT`).
        *   `make_term_tag()` ensures uniqueness by appending `#1`, `#2` if the base tag clashes with an existing one for another term.
        *   Sanitized by `sanitize()` (e.g., dots become underscores).
        *   If `term_tag` ends up being `NULL` (e.g., for `hrf(Ident(Var1, Var2))` without an `id`), then `make_column_names` omits the prefix, leading to direct names like `Var1_b01`.
    *   **`condition_tag`**:
        *   Represents the specific combination of factor levels or continuous regressor components *within* a term.
        *   Generated by `conditions.event_term()`, which combines:
            *   For factors: `VarName.LevelName` (e.g., `condition.A`).
            *   For continuous/basis (`ParametricBasis` like `Poly`, `BSpline`, `Ident`): The "levels" of the basis (e.g., `01`, `02` for `Poly(X,2)`; `Var1`, `Var2` for `Ident(Var1,Var2)`). These come from `levels.ParametricBasis` via `columns.ParametricBasis`.
        *   Interactions within a term are formed by `_` (underscore) joining these tokens (e.g., `condition.A_Task.Go`, `factorA.L1_PolyRT.01`). This is handled by `make_cond_tag()`.
    *   **`_b##` (Basis Suffix)**:
        *   Added by `add_basis()` if `nbasis(hrf_object) > 1`.
        *   Uses `basis_suffix()` which calls `zeropad()` for consistent numbering (e.g., `_b01`, `_b02`). For `nbasis=1`, no suffix is added.
*   **Assembly:** `make_column_names(term_tag, condition_tags, nbasis)` is the central function in `naming-utils.R` that puts these pieces together.
*   **Global Uniqueness:** While `make_term_tag` ensures `term_tag`s are unique, `build_event_model_design_matrix` performs a final `make.names(..., unique=TRUE)` on all generated column names as a safeguard against unforeseen clashes.

**IV. Baseline Model (`baseline_model.R`):**

*   Models non-event-related variance.
*   `baseline_model(basis, degree, sframe, intercept, nuisance_list)` is the constructor.
    *   `basis`: "constant", "poly", "bs" (B-spline), "ns" (natural spline) for drift modeling.
    *   `intercept`: "runwise" (separate intercept per run), "global" (one intercept), or "none".
    *   `nuisance_list`: A list of matrices (one per run) containing nuisance regressors (e.g., motion parameters).
*   Creates `baseline_term` objects for drift, block intercepts, and nuisance regressors.
*   `design_matrix.baseline_model()` combines these into the baseline design matrix.

**V. Full fMRI Model (`fmri_model.R`):**

*   `fmri_model(event_model, baseline_model)`: A simple container that combines the event-related and baseline model components.
*   `design_matrix.fmri_model()`: `cbind`s the design matrices from the `event_model` and `baseline_model`.
*   Provides methods for `terms()`, `conditions()`, `contrast_weights()`, plotting, etc., which typically delegate to or combine results from its constituent models.

**VI. Configuration DSL (`fmridsl.R`, `fmridsl_builder.R`):**

*   Provides a higher-level way to specify analyses using YAML configuration files.
*   `load_fmri_config(yaml_file)`:
    *   Parses the YAML.
    *   Applies default values.
    *   Performs schema validation (`parse_and_validate_config`).
    *   Performs context-dependent validation against a BIDS dataset (e.g., checking if specified subjects/tasks exist, if event columns are present) and infers variable roles (`build_config_from_ior`).
    *   Returns a validated `fmri_config` object.
*   `build_fmri_model_from_config(config, subject_id)`:
    *   Takes the `fmri_config` object.
    *   Loads subject-specific event and confound data.
    *   Constructs the `sampling_frame`, `baseline_model`, and `event_model` (by creating `hrfspec`s from the config and then following a similar pipeline to the direct `event_model` call).
    *   Returns a complete `fmri_model` object for that subject.

---

**Summary Diagram of Overall Model Building:**

```mermaid
graph TD
    subgraph User Input
        U_FORM[Formula & Data OR DSL Config File]
    end

    subgraph Core Processing
        HRF_SYS[HRF System: Define/Select/Decorate HRFs]
        EV_MODEL_PIPE[Event Model Pipeline: hrfspec -> event_term -> convolve -> DM_event]
        BASE_MODEL[Baseline Model: Drift, Intercept, Nuisance -> DM_baseline]
        NAMING[Naming Utilities: Systematic Column Names]
    end

    subgraph Output
        FMRI_MODEL[fmri_model Object]
        FULL_DM[Full Design Matrix]
        FIT[Model Fitting (fmri_lm, etc.)]
        CONTRASTS[Contrast Definition & Stats]
    end

    U_FORM --> HRF_SYS;
    U_FORM --> EV_MODEL_PIPE;
    U_FORM --> BASE_MODEL;

    HRF_SYS --> EV_MODEL_PIPE;
    NAMING --> EV_MODEL_PIPE;
    
    EV_MODEL_PIPE --> FMRI_MODEL;
    BASE_MODEL --> FMRI_MODEL;
    FMRI_MODEL --> FULL_DM;
    FULL_DM --> FIT;
    U_FORM --> CONTRASTS; %% User defines contrasts
    FIT & CONTRASTS --> CONTRASTS; %% Apply contrasts to fit
```

This architecture allows for considerable flexibility in defining complex fMRI models while maintaining a structured approach to regressor generation and naming. The separation of specification (e.g., `hrfspec`) from realization (e.g., `event_term`, convolved regressors) is a powerful design pattern.