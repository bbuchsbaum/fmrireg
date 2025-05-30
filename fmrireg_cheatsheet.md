## Micro-DSL (v2.6) & Output Format
Your output must be **Pure Markdown** that implies the structure defined by the Micro-DSL. **Do NOT output EBNF token names** (e.g., `H1_TOKEN`). Generate the actual Markdown (e.g., `# My Title`, `@f my_func(...)`).
**1. Markup Tokens (Implicitly Handled by Markdown):**
`H1` = `# text`
`H2` = `## text`
`NL` = Newline (use sparingly, primarily to end logical entries or separate blocks)
`HR` = `---` (use to separate major logical groups within a section, or between sections)
`Bul` = `- ` (dash + space for general bullets in Header/Legend/Deps)
`IndBul`= `  - ` (two spaces + dash + space for indented `Desc` lines under an `@sigil` entry)
**2. DSL Sigils:**
`@f` = Function `@d` = Data object (e.g., from `data()`)
`@g` = S4/S3 Generic `@m` = S4/S3 Method (class dispatch)
`@c` = R6/R7 constructor (if present)
**3. Other DSL Symbols:**
`|` = Separates multiple function names (e.g., `name1|name2`)
`[...]` = Used for:
* Constructor variants: `ConstructorName[VariantA|VariantB]`
* Method dispatch classes: `methodName[ClassA|ClassB]`
`(...)` = Parameter list in an entry signature.
`param?` = Optional parameter.
`param?=val`= Optional parameter with a default value.
`|` = Separates signature from short description (after params or name if no params).
`->` = Separates short description from return type. Omit if no return value (side-effect).
`!` = Prefix for inline notes (e.g., `!dep`, `!note: text`). There is no space between `!` and the note keyword.
**4. Document Skeleton:**
`Legend?` (H2 "Legend:" + type abbreviation table)
`Header` (H1 PackageName; optional H2 "Core Purpose:", H2 "Key Objects & Concepts:")
`Sections+` (H2 `1.` Title; H2 Title (unnumbered ok); optional H3 "Usage:")
`Entries*` (@sigil lines + optional indented bullets)
`Deps?` (H2 "Key R Dependencies:")
**5. Entry Line Structure:**
`@sigil name(|name)*[variant|ClassA|ClassB]? (alias alt1,alt2)? (param1?=val, param2?, ...)? | Short, pithy description -> ReturnTypeAbbr !note_type: Optional note text`
* **Rules for Entry Line:**
* Omit `()` if no parameters.
* Omit `-> ReturnTypeAbbr` if function has no return value (side-effect only).
* Bundle identical signatures using `name1|name2`.
* Use `ConstructorName[VariantA|VariantB]` for constructor subtypes.
* Use `methodName[DispatchClassA|DispatchClassB]` for S4/S3 methods.
* Notes (`!notetype: text` or `!notetype`) are optional postfixes. Ensure no leading space (e.g., `!ok` not `! ok`).
* Truncate parameter list with `...` if it exceeds 8 parameters (e.g., `(param1, param2, ..., param8, ...)`).
* Example of grouping aliases and optional params: `@f read_csv|read_tsv (file, col_types?="auto", ...) | Parse delimited file -> tib`
**6. Indented Description Bullets (`Desc` lines under an Entry):**
* Format: ` - param_name : type_abbr (constants: val1, "val2" | key_funcs: fnA, fnB)? Brief, essential clarification.`
* Include `(constants: ...)` for params that take a small, fixed set of string literals.
* Include `(key_funcs: ...)` for params that expect specific functions from the package as input.
* **Only include if adding significant clarity** beyond the signature. Omit for common/obvious parameters (e.g., `x`, `...`) or standard defaults (e.g., `drop=TRUE`).
* If an `Entry` already fits in ≤ 110 characters, do not add `Desc` lines unless they would prevent ambiguity.
* Can also be plain text for general notes: ` - General descriptive point.`
**7. Type Abbreviations:**
* Use very short (1-4 letter) type abbreviations (e.g., `NS` for NeuroSpace, `iv` for integer vector, `chr` for character, `log` for logical, `mat` for matrix, `obj` for generic S4/R6 object).
* Reuse type abbreviations whenever identical across entries; do not invent synonyms (e.g., use `int` consistently, not `int`, `intv`, `iv` interchangeably).
**7. Type Abbreviations (deprecated – see next subsection):**
* Use very short (1-4 char) codes, reusing consistently (e.g., `int` not `intv`).
* Built-in codes: int, dbl, num, chr, lgl, lst, vec, df, tib, tbl, mat, arr, fn, env, obj, NS
* Provide a `## Legend:` block only when introducing abbreviations beyond this list.
## Compression Heuristics & Content Selection:
1. **Focus:** Public, user-facing API. Omit internal helpers, unexported symbols, and direct S4 slot accessors (like `slotNames` or methods that just return a slot value if there's already a clear getter). Include only symbols present in NAMESPACE export list.
2. **Grouping:**
* Group trivial getters or functions with identical signatures and purpose using `name1|name2`.
* Group constructors with identical fields but different return subtypes using `ConstructorName[VariantA|VariantB]`.
* Group S3/S4 methods with identical implementations/docs using `methodName[ClassA|ClassB]`.
3. **Methods:** Define generics with `@g`. Emit methods (`@m`) **only** when their behavior, parameters, or return type significantly differ from the generic, or to explicitly list key supported classes.
4. **Omissions:** Skip indented parameter descriptions (`Desc` lines) for obvious defaults (e.g., `drop=TRUE`, `smooth=FALSE`) or very common arguments like `x` or `...` unless they have package-specific meaning. The cheatsheet is not full documentation.
5. **Notes:** Use `!` notes sparingly for critical info (e.g., `!dep`, `!imp`, `!retlist`, `!side`).
6. **Re-exports:** Skip functions and generics re-exported from other packages (e.g., if `dplyr::filter` is re-exported, do not list it).
## Output Contract:
* **Pure Markdown only.** Adhere strictly to the Micro-DSL v2.6.
* No commentary, no intro/outro paragraphs, no code fences unless part of a `CODE_BLOCK_TOKEN` within a `BlockContent` (rarely needed for cheatsheets).
* If any line violates the DSL, regenerate until fully compliant—no prose explanations.
* Generation stops at first line that begins with a second `#` at H1 depth (e.g., if `\\n#` is used as a stop sequence).
* Use a single blank line to separate `Entry` blocks if it aids readability, but avoid excessive blank lines. Use `---` (`HR_TOKEN`) to separate major thematic groups within a section or at the end of sections.
* All parens/brackets/pipes must be balanced.
## Self-Check (Mental Step - Crucial):
Before finalizing, review your output against these critical checks:
1. Does every content line belong to a defined DSL structure (Header, Legend, Section, Entry, Desc, Deps, Block bullet)?
2. Is the DSL syntax for `Entry` lines (sigils, names, params, `|`, `->`, `!`) correctly used? No bare `->` tokens.
3. (reserved)
4. For any abbreviation NOT in the built-in list, is it defined in `## Legend:`?
5. Have you omitted non-essential details and internal functions?

## Few-shot Exemplar
```markdown
# dummyPkg
## Legend:
- int : integer
- chr : character
## 1. Core Functions
@f add (x, y) | Sum two ints -> int
```

## Formal Grammar "v2.6 Micro-EBNF"
Cheatsheet ::= Header Legend? Section+ Deps?
Legend ::= H2_TOKEN TEXT_CONTENT NEWLINE_TOKEN Block
Header ::= H1_TOKEN TEXT_CONTENT NEWLINE_TOKEN+ (H2Section)*
H2Section ::= H2_TOKEN TEXT_CONTENT NEWLINE_TOKEN Block
Section ::= H2_TOKEN (NUMBER_TOKEN PERIOD_TOKEN)? TEXT_CONTENT NEWLINE_TOKEN UsageBlock? Entry+ HR_TOKEN?
UsageBlock ::= H3_TOKEN "Usage:" NEWLINE_TOKEN Block
Entry ::= Sigil_TOKEN EntryIdent ParamList? Bar_TOKEN TEXT_CONTENT ArrowReturn Note? NEWLINE_TOKEN Desc*
Sigil_TOKEN ::= AT_F_TOKEN | AT_D_TOKEN | AT_G_TOKEN | AT_M_TOKEN // Lexer provides @f, @d, @g, @m
EntryIdent ::= IdentGroup MethodOrVariantClass? AliasSpec?
IdentGroup ::= IDENT_TOKEN ("|" IDENT_TOKEN)* // For "foo|bar"
MethodOrVariantClass ::= LBRACKET_TOKEN IDENT_TOKEN (PIPE_TOKEN IDENT_TOKEN)* RBRACKET_TOKEN // For "[ClassA|ClassB]" or "[variantA|variantB]"
AliasSpec ::= LPAREN_TOKEN ALIAS_KEYWORD_TOKEN IDENT_TOKEN (COMMA_TOKEN IDENT_TOKEN)* RPAREN_TOKEN // For "(alias alt1, alt2)"
ParamList ::= LPAREN_TOKEN Param (COMMA_TOKEN Param)* RPAREN_TOKEN
Param ::= IDENT_TOKEN (EQUALS_TOKEN DefaultValue)? OPTIONAL_MARKER_TOKEN?
DefaultValue ::= LITERAL_TOKEN | IDENT_TOKEN
ArrowReturn ::= (ARROW_TOKEN IDENT_TOKEN)?
Note ::= EXCLAMATION_TOKEN NOTETYPE_TOKEN (COLON_TOKEN TEXT_CONTENT)? NEWLINE_TOKEN?
Desc ::= INDENT_TOKEN BULLET_MARKER_TOKEN (ParamDesc | TEXT_CONTENT) NEWLINE_TOKEN
ParamDesc ::= IDENT_TOKEN COLON_TOKEN TYPE_ABBR_TOKEN ParamExtra? TEXT_CONTENT?
ParamExtra ::= LPAREN_TOKEN (ConstantsSpec | KeyFuncsSpec) RPAREN_TOKEN
ConstantsSpec ::= "constants:" (IDENT_TOKEN|LITERAL_TOKEN) (COMMA_TOKEN (IDENT_TOKEN|LITERAL_TOKEN))*
KeyFuncsSpec ::= "key_funcs:" IDENT_TOKEN (COMMA_TOKEN IDENT_TOKEN)*
Deps ::= H2_TOKEN "Key R Dependencies:" NEWLINE_TOKEN Block
Block ::= (Bullet | TEXT_CONTENT | CODE_BLOCK_TOKEN)* (NEWLINE_TOKEN | EOF_TOKEN)
Bullet ::= BULLET_MARKER_TOKEN TEXT_CONTENT NEWLINE_TOKEN
/* --- LEXER-IMPLIED TOKENS (Illustrative) ---
All previous tokens from v2.4, plus:
AT_G_TOKEN, AT_M_TOKEN // @g, @m
// The lexer provides LBRACKET_TOKEN, IDENT_TOKEN, PIPE_TOKEN, RBRACKET_TOKEN.
// The parser, guided by the Sigil_TOKEN (@f for constructor variants, @m for method classes),
// will interpret the content of MethodOrVariantClass appropriately.
*/

---

# fmrireg

## 1. Exported Data Objects: HRF Types & Benchmark Datasets
### Usage:
- Use HRF constants (e.g., `HRF_SPMG1`, `HRF_GAMMA`, `HRF_BSPLINE`) as the `basis` argument in `hrf()` or directly in model formulas.
- Benchmark datasets can be loaded with `load_benchmark_dataset()`; see available names with `list_benchmark_datasets()`.
@d HRF_SPMG1 | SPM canonical HRF (single basis) -> fn
@d HRF_SPMG2 | SPM canonical HRF + temporal derivative (2 basis) -> fn
@d HRF_SPMG3 | SPM canonical HRF + temporal & dispersion derivatives (3 basis) -> fn
@d HRF_GAMMA | Gamma function HRF -> fn
@d HRF_GAUSSIAN | Gaussian function HRF -> fn
@d HRF_BSPLINE | B-spline basis HRF (5 basis) -> fn
@d HRF_FIR | FIR basis HRF (12 basis) -> fn
@d fmri_benchmark_datasets | List of simulated fMRI datasets for benchmarking -> lst

---

## 2. HRF & Basis Constructors
### Usage:
- Use `hrf()` in model formulas to specify event regressors with a chosen basis (e.g., `hrf(condition, basis=HRF_SPMG1)`).
- For custom HRFs, use `as_hrf()` to wrap a function.
- Use `bind_basis()` to combine multiple HRFs into a basis set.
@f hrf (..., basis="spmg1", onsets?=NULL, durations?=NULL, prefix?=NULL, subset?=NULL, precision?=.3, nbasis?=1, contrasts?=NULL, id?=NULL, name?=NULL, lag?=0, summate?=TRUE) | Specify HRF regressor for model formula -> obj
@f as_hrf (f, name?=NULL, nbasis?=1, span?=24, params?=list()) | Wrap function as HRF object -> fn
@f bind_basis (...) | Combine multiple HRFs into a basis set -> fn
@f gen_hrf (hrf, lag?=0, width?=0, precision?=.1, summate?=TRUE, normalize?=FALSE, name?=NULL, span?=NULL, ...) | Compose HRF with lag/block/normalization -> fn
@f lag_hrf (hrf, lag) | Lag an HRF object -> fn
@f block_hrf (hrf, width, precision?=0.1, half_life?=Inf, summate?=TRUE, normalize?=FALSE) | Block (boxcar) HRF decorator -> fn
@f normalise_hrf (hrf) | Normalize HRF to peak 1 -> fn
@f hrf_set (...) | Combine HRFs into a set -> fn
@f empirical_hrf (t, y, name?="empirical_hrf") | Create HRF from empirical data -> fn
@f hrf_library (fun, pgrid, ...) | Create HRF library from parameter grid -> fn

---

## 3. Model & Design Construction
### Usage:
- Build event models with `event_model()` using formulas and HRF specs.
- Combine event and baseline models with `fmri_model()`.
- Use `baseline_model()` for drift/block/nuisance terms.
- Extract design matrices with `design_matrix()`.
@f event_model (formula_or_list, data, block, sampling_frame, durations?=0, drop_empty?=TRUE, precision?=0.3, parallel?=FALSE, progress?=FALSE, ...) | Create event-based fMRI model -> obj
@f fmri_model (event_model, baseline_model) | Combine event and baseline into fmri_model -> obj
@f baseline_model (basis?="constant", degree?=1, sframe, intercept?="runwise", nuisance_list?=NULL) | Construct baseline model (drift, block, nuisance) -> obj
@f baseline (degree?=1, basis?="constant", name?=NULL, intercept?="runwise") | Baseline drift term specification -> obj
@f block (x) | Block variable specification -> obj
@f nuisance (x) | Nuisance term specification -> obj
@f design_matrix (x, ...) | Extract design matrix from model/term -> mat
@f design_map (x, ...) | Plot design matrix heatmap -> obj
@f correlation_map (x, ...) | Plot design matrix correlation heatmap -> obj
@f design_plot (fmrimod, ...) | Interactive design matrix plot (shiny) -> obj
@f term_matrices (x, ...) | Extract per-term design matrices -> lst

---

## 4. Dataset Constructors & Utilities
### Usage:
- Use `matrix_dataset()`, `fmri_dataset()`, `fmri_mem_dataset()`, `latent_dataset()` to create datasets.
- Use `get_data_matrix()` to extract data as matrix.
@f matrix_dataset (datamat, TR, run_length, event_table?=data.frame()) | Create matrix-format fMRI dataset -> obj
@f fmri_dataset (scans, mask, TR, run_length, event_table?=data.frame(), base_path?=".", censor?=NULL, preload?=FALSE, mode?="normal") | Create file-based fMRI dataset -> obj
@f fmri_mem_dataset (scans, mask, TR, run_length, event_table?=data.frame(), base_path?=".", censor?=NULL) | Create in-memory fMRI dataset -> obj
@f latent_dataset (lvec, TR, run_length, event_table?=data.frame()) | Create latent (reduced) dataset -> obj
@f get_data_matrix (x, ...) | Extract data matrix from dataset -> mat
@f get_mask (x, ...) | Extract mask from dataset -> obj
@f as.matrix_dataset (x, ...) | Convert to matrix_dataset -> obj
@f data_chunks (x, nchunks?=1, runwise?=FALSE, ...) | Iterator over data chunks -> obj

---

## 5. Model Fitting & Estimation
### Usage:
- Fit models with `fmri_lm()` (standard) or `fmri_rlm()` (robust).
- Use `estimate_betas()` for beta estimation (various methods).
- Use `glm_ols()` for condition-level, `glm_lss()` for single-trial estimation.
@f fmri_lm (formula, block, baseline_model?=NULL, dataset, durations?=0, drop_empty?=TRUE, robust?=FALSE, robust_options?=NULL, ar_options?=NULL, strategy?="runwise", nchunks?=10, use_fast_path?=FALSE, progress?=FALSE, ...) | Fit fMRI linear model (GLM) -> obj
@f fmri_rlm (formula, block, baseline_model?=NULL, dataset, durations?=0, drop_empty?=TRUE, strategy?="runwise", nchunks?=10, cor_struct?="iid", cor_iter?=1L, cor_global?=FALSE, ar_p?=NULL, ar1_exact_first?=FALSE, robust_psi?="huber", robust_k_huber?=1.345, robust_c_tukey?=4.685, robust_max_iter?=2L, robust_scale_scope?="run", ...) | Fit robust fMRI linear model -> obj
@f estimate_betas (x, fixed?=NULL, ran, block, method?="mixed", basemod?=NULL, ...) | Estimate betas (various methods) -> obj
  - method : chr (constants: "mixed", "mixed_cpp", "lss", "lss_naive", "lss_cpp", "pls", "pls_global", "ols") Estimation method.
@f glm_ols (dataset, model_obj, basis_obj, basemod?=NULL, block?=~1, progress?=TRUE, ...) | OLS estimation (condition-level) -> obj
@f glm_lss (dataset, model_obj, basis_obj, basemod?=NULL, block?=~1, use_cpp?=TRUE, progress?=TRUE, ...) | LSS estimation (single-trial) -> obj
@f fmri_lm_control (robust_options?=list(), ar_options?=list()) | Create config for model fitting -> obj
@f fmri_lm_config (robust?=FALSE, ar_options?=NULL, weights?=NULL, method?="auto") | Wrapper for fmri_lm_control -> obj
@f read_fmri_config (file_name, base_path?=NULL) | Read fMRI config file -> obj

---

## 6. HRF Evaluation & Manipulation
### Usage:
- Evaluate HRFs or regressors at time points with `evaluate()`.
- Use `reconstruction_matrix()` to get basis-to-time mapping.
- Use `hrf_from_coefficients()` to create HRF from basis weights.
@f evaluate (x, grid, ...) | Evaluate HRF or regressor at time points -> vec|mat
@f reconstruction_matrix (hrf, sframe, ...) | Get HRF basis-to-time matrix -> mat
@f hrf_from_coefficients (hrf, h, ...) | Combine HRF basis with coefficients -> fn
@f penalty_matrix (x, ...) | Generate penalty matrix for HRF basis -> mat

---

## 7. Event & Regressor Constructors
### Usage:
- Use `regressor()` to create event regressors; `single_trial_regressor()` for single events.
- Use `event_term()` to build event terms from variables, onsets, etc.
@f regressor (onsets, hrf?=HRF_SPMG1, duration?=0, amplitude?=1, span?=40, summate?=TRUE) | Create event regressor -> obj
@f single_trial_regressor (onsets, hrf?=HRF_SPMG1, duration?=0, amplitude?=1, span?=24) | Single-trial regressor -> obj
@f event_term (evlist, onsets, blockids, durations?=0, subset?=NULL) | Create event term from variables -> obj
@f event_factor (fac, name, onsets, blockids?=1, durations?=0, subset?=NULL) | Event term from factor -> obj
@f event_variable (vec, name, onsets, blockids?=1, durations?=0, subset?=NULL) | Event term from numeric vector -> obj
@f event_matrix (mat, name, onsets, blockids?=1, durations?=0, subset?=NULL) | Event term from matrix -> obj
@f event_basis (basis, name?=NULL, onsets, blockids?=1, durations?=0, subset?=NULL) | Event term from basis -> obj
@f trialwise (basis?="spmg1", lag?=0, nbasis?=1, add_sum?=FALSE, label?="trial") | Trialwise HRF spec for single-trial GLM -> obj

---

## 8. Model Inspection, Extraction, and Visualization
### Usage:
- Use `coef()`, `stats()`, `standard_error()`, `p_values()` to extract results from fitted models.
- Use `fitted_hrf()` to get fitted HRF shapes.
- Use `plot_contrasts()` to visualize contrast weights.
@f coef (object, type?="betas", include_baseline?=FALSE, recon?=FALSE, ...) | Extract coefficients from fit -> mat|tib
@f stats (x, type?="estimates", ...) | Extract test statistics from fit -> mat|tib
@f standard_error (x, type?="estimates", ...) | Extract standard errors from fit -> mat|tib
@f p_values (x, ...) | Extract p-values from fit -> mat|tib
@f fitted_hrf (x, sample_at, ...) | Get fitted HRF from fit -> mat|tib
@f plot_contrasts (x, ...) | Visualize contrast weights -> obj
@f autoplot (object, ...) | Plot regressor or model object -> obj

---

## 9. Contrast Specification & Utilities
### Usage:
- Use `contrast()`, `pair_contrast()`, `unit_contrast()`, `oneway_contrast()`, `interaction_contrast()`, `poly_contrast()`, `column_contrast()` to define contrasts.
- Use `contrast_weights()` to compute weights for a model term.
- Use `contrast_set()` to group multiple contrasts.
@f contrast (form, name, where?=NULL) | Define linear contrast -> obj
@f pair_contrast (A, B, name, where?=NULL) | Pairwise contrast (A vs B) -> obj
@f unit_contrast (A, name, where?=NULL) | Unit contrast (sum-to-1) -> obj
@f one_against_all_contrast (levels, facname, where?=NULL) | Each level vs all others -> obj
@f one_way_contrast (A, name, where?=NULL) | One-way contrast (main effect) -> obj
@f interaction_contrast (A, name, where?=NULL) | Interaction contrast -> obj
@f poly_contrast (A, name, where?=NULL, degree?=1, value_map?=NULL) | Polynomial trend contrast -> obj
@f column_contrast (pattern_A, pattern_B?=NULL, name, where?=NULL) | Contrast by column regex -> obj
@f contrast_set (...) | Group multiple contrast specs -> obj
@f contrast_weights (x, ...) | Compute contrast weights for term -> lst
@f Fcontrasts (x, ...) | Generate F-contrasts for term -> lst
@f to_glt (x, ...) | Convert contrast to AFNI GLT -> obj
@f write_glt (x, ...) | Write GLT to file

---

## 10. Event/Design Inspection & Extraction
### Usage:
- Use `cells()`, `conditions()`, `event_table()`, `elements()`, `onsets()`, `durations()`, `amplitudes()`, `blockids()`, `blocklens()`, `samples()`, `split_by_block()`, `split_onsets()` for design/model inspection.
@f cells (x, ...) | Get experimental cells (factor level combos) -> tib
@f conditions (x, ...) | Get condition labels for term/model -> chr
@f event_table (x, ...) | Get event table for term/model -> tib
@f elements (x, ...) | Get ordered elements of term/variable -> vec
@f onsets (x, ...) | Get event onsets -> vec
@f durations (x, ...) | Get event durations -> vec
@f amplitudes (x, ...) | Get event amplitudes -> vec
@f blockids (x, ...) | Get block/run indices -> vec
@f blocklens (x, ...) | Get block/run lengths -> vec
@f samples (x, ...) | Get sampling times -> vec
@f split_by_block (x, ...) | Split data by block/run -> lst
@f split_onsets (x, sframe, global?=FALSE, blocksplit?=FALSE, ...) | Split onsets by factor/block -> lst
@f longnames (x, ...) | Get full condition names (term+level+basis) -> chr
@f shortnames (x, ...) | Get short condition names (levels only) -> chr

---

## 11. Simulation & Benchmarking
### Usage:
- Use `simulate_bold_signal()`, `simulate_noise_vector()`, `simulate_simple_dataset()`, `simulate_fmri_matrix()` for simulation.
- Use `load_benchmark_dataset()`, `list_benchmark_datasets()`, `get_benchmark_summary()`, `create_design_matrix_from_benchmark()` for benchmarking.
@f simulate_bold_signal (ncond, hrf?=HRF_SPMG1, nreps?=12, amps?=rep(1,ncond), ampsd?=0, isi?=c(3,6), TR?=1.5) | Simulate fMRI BOLD signal -> lst
@f simulate_noise_vector (n, TR?=1.5, ar?=c(0.3), ma?=c(0.5), sd?=1, drift_freq?=1/128, drift_amplitude?=2, physio?=TRUE, seed?=NULL) | Simulate fMRI noise vector -> vec
@f simulate_simple_dataset (ncond, nreps?=12, TR?=1.5, snr?=0.5, hrf?=HRF_SPMG1, seed?=NULL) | Simulate full fMRI dataset (signal+noise) -> lst
@f simulate_fmri_matrix (...) | Simulate fMRI time series matrix with per-column variability -> lst
@f load_benchmark_dataset (dataset_name?="BM_Canonical_HighSNR") | Load benchmark dataset -> lst
@f list_benchmark_datasets () | List available benchmark datasets -> df
@f get_benchmark_summary (dataset_name) | Get summary of benchmark dataset -> lst
@f create_design_matrix_from_benchmark (dataset_name, hrf, include_intercept?=TRUE) | Create design matrix from benchmark dataset -> mat
@f evaluate_method_performance (dataset_name, estimated_betas, method_name?="Unknown") | Evaluate method on benchmark dataset -> lst

---

## 12. AFNI Integration
### Usage:
- Use `afni_hrf()` and `afni_trialwise()` to specify AFNI-compatible HRF terms.
- Use `afni_lm()` to set up AFNI 3dDeconvolve models.
- Use `build_afni_stims()` to generate AFNI stimulus files.
@f afni_hrf (..., basis, onsets?=NULL, durations?=NULL, prefix?=NULL, subset?=NULL, nbasis?=1, contrasts?=NULL, id?=NULL, lag?=0, precision?=0.3, summate?=TRUE, start?=NULL, stop?=NULL) | AFNI HRF specification for 3dDeconvolve -> obj
@f afni_trialwise (label, basis, onsets?=NULL, durations?=0, subset?=NULL, id?=NULL, start?=0, stop?=22, precision?=0.3, summate?=TRUE) | AFNI trialwise HRF spec for -stim_times_IM -> obj
@f afni_lm (fmri_mod, dataset, working_dir?=".", polort?=-1, jobs?=1, censor?=NULL, options?=list()) | Prepare AF