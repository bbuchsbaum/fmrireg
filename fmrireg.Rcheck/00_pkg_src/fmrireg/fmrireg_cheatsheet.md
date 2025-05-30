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

## 1. Core Data Objects & Constants
### Usage:
- Use exported HRF objects (e.g., `HRF_SPMG1`, `HRF_SPMG2`, `HRF_SPMG3`, `HRF_GAMMA`, `HRF_GAUSSIAN`, `HRF_BSPLINE`, `HRF_FIR`) as the `basis` argument in `hrf()` or directly in model specification.
- Use AFNI HRF objects (e.g., `AFNI_HRF`) for AFNI-specific workflows.
- Benchmark datasets: use `load_benchmark_dataset()`, `list_benchmark_datasets()`, `get_benchmark_summary()` for simulation/testing.

@d HRF_SPMG1 | SPM canonical HRF (single basis) -> fn
@d HRF_SPMG2 | SPM canonical HRF + temporal derivative (2 basis) -> fn
@d HRF_SPMG3 | SPM canonical HRF + temporal & dispersion derivatives (3 basis) -> fn
@d HRF_GAMMA | Gamma function HRF -> fn
@d HRF_GAUSSIAN | Gaussian function HRF -> fn
@d HRF_BSPLINE | B-spline basis HRF (default 5 basis) -> fn
@d HRF_FIR | FIR basis HRF (default 12 basis) -> fn
@d AFNI_HRF | AFNI-specific HRF object -> obj
@d fmri_benchmark_datasets | List of simulated benchmark fMRI datasets -> lst

---

## 2. Model Construction & Design
### Usage:
- Build event models with `event_model(formula, data, block, sampling_frame, ...)`.
- Combine with a baseline model via `fmri_model(event_model, baseline_model)`.
- Use `baseline_model()` to specify drift, block, and nuisance terms.
- For AFNI workflows, use `afni_hrf()` or `afni_trialwise()` as basis in `event_model`.
- Use exported HRF objects or pass a custom HRF via `basis` in `hrf()`.

@f event_model (formula_or_list, data, block, sampling_frame, durations?=0, drop_empty?=TRUE, ...) | Construct event-based fMRI model -> obj
  - formula_or_list : formula or lst Model formula or list of hrfspecs.
  - block : formula or vec Block/run structure.
  - sampling_frame : obj (key_funcs: sampling_frame) Temporal/block structure.
  - durations : num (default: 0) Event durations.

@f fmri_model (event_model, baseline_model) | Combine event and baseline models -> obj
  - event_model : obj (key_funcs: event_model) Event model.
  - baseline_model : obj (key_funcs: baseline_model) Baseline model.

@f baseline_model (basis?="constant", degree?=1, sframe, intercept?="runwise", nuisance_list?=NULL) | Construct baseline model -> obj
  - basis : chr (constants: "constant", "poly", "bs", "ns") Drift basis.
  - degree : int Degree of basis.
  - sframe : obj (key_funcs: sampling_frame) Sampling frame.
  - intercept : chr (constants: "runwise", "global", "none") Intercept type.

@f hrf (..., basis?="spmg1", onsets?=NULL, durations?=NULL, prefix?=NULL, subset?=NULL, precision?=0.3, nbasis?=1, contrasts?=NULL, id?=NULL, name?=NULL, lag?=0, summate?=TRUE) | Specify HRF term for model formula -> obj
  - basis : chr or fn (constants: "spmg1", "spmg2", "spmg3", "gamma", "gaussian", "bspline", "tent", "fourier", "fir", HRF_SPMG1, ...) HRF basis.

@f afni_hrf (..., basis, onsets?=NULL, durations?=NULL, prefix?=NULL, subset?=NULL, nbasis?=1, contrasts?=NULL, id?=NULL, lag?=0, precision?=0.3, summate?=TRUE, start?=NULL, stop?=NULL) | AFNI HRF spec for AFNI workflows -> obj

@f afni_trialwise (label, basis, onsets?=NULL, durations?=0, subset?=NULL, id?=NULL, start?=0, stop?=22, precision?=0.3, summate?=TRUE) | AFNI trialwise HRF spec -> obj

@f covariate (..., data, id?=NULL, prefix?=NULL, subset?=NULL) | Add non-convolved covariate term -> obj

@f sampling_frame (blocklens, TR, start_time?=TR/2, precision?=0.1) | Define scan/block structure -> obj

---

## 3. Data Input, Simulation, and Benchmarking
### Usage:
- Use `fmri_dataset()`, `matrix_dataset()`, `fmri_mem_dataset()`, `latent_dataset()` to create datasets.
- For simulation/benchmarking, use `simulate_bold_signal()`, `simulate_noise_vector()`, `simulate_simple_dataset()`, `simulate_fmri_matrix()`.
- Access built-in datasets with `load_benchmark_dataset()`, `list_benchmark_datasets()`, `get_benchmark_summary()`.

@f fmri_dataset (scans, mask, TR, run_length, event_table?=df, base_path?=".", censor?=NULL, preload?=FALSE, mode?="normal") | Create fMRI dataset from files -> obj
@f matrix_dataset (datamat, TR, run_length, event_table?=df) | Create dataset from matrix -> obj
@f fmri_mem_dataset (scans, mask, TR, run_length, event_table?=df, base_path?=".", censor?=NULL) | In-memory fMRI dataset -> obj
@f latent_dataset (lvec, TR, run_length, event_table?=df) | Latent component dataset -> obj

@f simulate_bold_signal (ncond, hrf?=HRF_SPMG1, nreps?=12, amps?=rep(1,ncond), isi?=c(3,6), TR?=1.5, ...) | Simulate BOLD time series -> lst
@f simulate_noise_vector (n, TR?=1.5, ar?=c(0.3), ma?=c(0.5), sd?=1, ...) | Simulate fMRI noise vector -> vec
@f simulate_simple_dataset (ncond, nreps?=12, TR?=1.5, snr?=0.5, hrf?=HRF_SPMG1, seed?=NULL) | Simulate complete fMRI dataset -> lst
@f simulate_fmri_matrix (n, total_time, TR, hrf, n_events, ...) | Simulate multi-column fMRI matrix -> lst

@f load_benchmark_dataset (dataset_name?="BM_Canonical_HighSNR") | Load benchmark dataset -> lst
@f list_benchmark_datasets () | List available benchmark datasets -> df
@f get_benchmark_summary (dataset_name) | Get summary of benchmark dataset -> lst
@f create_design_matrix_from_benchmark (dataset_name, hrf, include_intercept?=TRUE) | Create design matrix from benchmark -> mat
@f evaluate_method_performance (dataset_name, estimated_betas, method_name?="Unknown") | Evaluate method on benchmark -> lst

---

## 4. Model Fitting & Estimation
### Usage:
- Fit models with `fmri_lm()` (OLS), `fmri_rlm()` (robust), or `fmri_latent_lm()` (latent).
- For beta estimation, use `estimate_betas()` (supports methods: "mixed", "lss", "pls", "ols", etc.).
- Use `glm_ols()` for condition-level, `glm_lss()` for single-trial estimation.
- Use `gen_afni_lm()` and `afni_lm()` for AFNI 3dDeconvolve workflows.

@f fmri_lm (formula, block, baseline_model?=NULL, dataset, durations?=0, drop_empty?=TRUE, robust?=FALSE, strategy?="runwise", nchunks?=10, use_fast_path?=FALSE, progress?=FALSE, ...) | Fit fMRI linear model (OLS) -> obj
@f fmri_rlm (formula, block, baseline_model?=NULL, dataset, durations?=0, drop_empty?=TRUE, strategy?="runwise", nchunks?=10, ...) | Fit robust fMRI linear model -> obj
@f fmri_latent_lm (formula, block, baseline_model?=NULL, dataset, durations, drop_empty?=TRUE, robust?=FALSE, autocor?="none", bootstrap?=FALSE, nboot?=1000, ...) | Fit model to latent dataset -> obj

@f estimate_betas (x, fixed?=NULL, ran, block, method?="mixed", basemod?=NULL, ...) | Estimate beta coefficients -> lst
  - method : chr (constants: "mixed", "mixed_cpp", "lss", "lss_naive", "lss_cpp", "pls", "pls_global", "ols") Estimation method.

@f glm_ols (dataset, model_obj, basis_obj, basemod?=NULL, block?=~1, progress?=TRUE, ...) | OLS condition-level estimation -> lst
@f glm_lss (dataset, model_obj, basis_obj, basemod?=NULL, block?=~1, use_cpp?=TRUE, progress?=TRUE, ...) | LSS single-trial estimation -> lst

@f gen_afni_lm (x, ...) | Generate AFNI linear model from config -> obj
@f afni_lm (fmri_mod, dataset, working_dir?=".", polort?=-1, jobs?=1, censor?=NULL, options?=list()) | Prepare AFNI 3dDeconvolve model -> obj

---

## 5. Design Matrix, Regressors, and Visualization
### Usage:
- Extract design matrices with `design_matrix()`, visualize with `design_map()` or `correlation_map()`.
- Use `regressor()` or `single_trial_regressor()` to create regressor objects.
- Evaluate regressors with `evaluate()`.
- Plot regressors with `autoplot()`.

@f design_matrix (x, ...) | Extract design matrix -> mat
@f design_map (x, ...) | Plot design matrix heatmap -> obj
@f correlation_map (x, ...) | Plot design matrix correlation heatmap -> obj
@f regressor (onsets, hrf?=HRF_SPMG1, duration?=0, amplitude?=1, span?=40, summate?=TRUE) | Create regressor object -> obj
@f single_trial_regressor (onsets, hrf?=HRF_SPMG1, duration?=0, amplitude?=1, span?=24) | Single-trial regressor -> obj
@f evaluate (x, grid, ...) | Evaluate regressor or HRF -> vec|mat
@f autoplot (object, grid?=NULL, precision?=0.1, method?="conv", ...) | Plot regressor object -> obj

---

## 6. Contrasts & Statistical Inference
### Usage:
- Define contrasts with `contrast()`, `unit_contrast()`, `pair_contrast()`, `oneway_contrast()`, `interaction_contrast()`, `column_contrast()`, `poly_contrast()`, `contrast_set()`, `pairwise_contrasts()`, `one_against_all_contrast()`.
- Compute contrast weights with `contrast_weights()`, F-contrasts with `Fcontrasts()`.
- Estimate contrasts with `estimate_contrast()`, `to_glt()`, `write_glt()`.
- Visualize with `plot_contrasts()`.

@f contrast (form, name, where?=NULL) | Define linear contrast -> obj
@f unit_contrast (A, name, where?=NULL) | Define unit contrast -> obj
@f pair_contrast (A, B, name, where?=NULL) | Define pairwise contrast -> obj
@f oneway_contrast (A, name, where?=NULL) | Define one-way contrast -> obj
@f interaction_contrast (A, name, where?=NULL) | Define interaction contrast -> obj
@f column_contrast (pattern_A, pattern_B?=NULL, name, where?=NULL) | Define column-based contrast -> obj
@f poly_contrast (A, name, where?=NULL, degree?=1, value_map?=NULL) | Define polynomial contrast -> obj
@f contrast_set (...) | Set of contrasts -> lst
@f pairwise_contrasts (levels, facname, where?=NULL, name_prefix?="con") | All pairwise contrasts -> obj
@f one_against_all_contrast (levels, facname, where?=NULL) | Each level vs all others -> obj

@f contrast_weights (x, ...) | Compute contrast weights -> lst
@f Fcontrasts (x, ...) | Compute F-contrasts -> lst
@f estimate_contrast (x, fit, colind, ...) | Estimate contrast from fit -> lst
@f to_glt (x, ...) | Convert contrast to AFNI GLT -> obj
@f write_glt (x, ...) | Write GLT to file

@f plot_contrasts (x, ...) | Visualize contrasts as heatmap -> obj

---

## 7. HRF Utilities & Advanced
### Usage:
- Use `gen_hrf()`, `gen_hrf_lagged()`, `gen_hrf_blocked()`, `hrf_set()`, `hrf_library()` for custom HRF construction.
- Decorate HRFs with `lag_hrf()`, `block_hrf()`, `normalise_hrf()`.
- List available HRFs with `list_available_hrfs()`.

@f gen_hrf (hrf, lag?=0, width?=0, precision?=0.1, summate?=TRUE, normalize?=FALSE, name?=NULL, span?=NULL, ...) | Compose HRF with lag/block/normalization -> fn
@f gen_hrf_lagged (hrf, lag?=2, normalize?=FALSE, ...) | Lagged HRF generator -> fn
@f gen_hrf_blocked (hrf?=hrf_gaussian, width?=5, precision?=0.1, half_life?=Inf, summate?=TRUE, normalize?=FALSE, ...) | Blocked HRF generator -> fn
@f hrf_set (..., name?="hrf_set") | Combine HRFs into basis set -> fn
@f hrf_library (fun, pgrid, ...) | HRF library from parameter grid -> fn
@f lag_hrf (hrf, lag) | Lag HRF object -> fn
@f block_hrf (hrf, width, precision?=0.1, half_life?=Inf, summate?=TRUE, normalize?=FALSE) | Block HRF object -> fn
@f normalise_hrf (hrf) | Normalize HRF object -> fn
@f list_available_hrfs (details?=FALSE) | List available HRFs -> df

---

## 8. Miscellaneous & Helper Functions
@f design_plot (fmrimod, term_name?=NULL, longnames?=FALSE, plot_title?=NULL, x_label?="Time (s)", y_label?="Amplitude", line_size?=1, color_palette?="viridis", facet_ncol?=2, theme_custom?=NULL, legend_threshold?=30, ...) | Interactive design plot (Shiny) -> obj
@f translate_legacy_pattern (pattern) | Convert legacy contrast regex to new naming -> chr
@f penalty_matrix (x, ...) | Generate penalty matrix for HRF/basis -> mat
@f hrf_smoothing_kernel (len, TR?=2, form?=onset~trialwise(), buffer_scans?=3L, normalise?=TRUE, method?="gram") | Compute HRF smoothing kernel -> mat

---

## 9. Data Accessors & Utilities
@f get_data (x, ...) | Extract data from dataset -> mat|obj
@f get_data_matrix (x, ...) | Extract data matrix from dataset -> mat
@f get_mask (x, ...) | Extract mask from dataset -> mat|vec
@f as.matrix_dataset (x, ...) | Convert to matrix_dataset -> obj
@f data_chunks (x, nchunks?=1, runwise?=FALSE, ...) | Iterator over data chunks -> obj
@f samples (x, ...) | Get sampling times -> vec
@f blockids (x, ...) | Get block/run indices -> vec
@f blocklens (x, ...) | Get block/run lengths -> vec
@f split_by_block (x, ...) | Split data by block/run -> lst

---

## 10. Statistical Output & Summaries
@f coef (object, ...) | Extract model coefficients -> mat|tib
@f stats (x, ...) | Extract test statistics -> mat|tib
@f standard_error (x, ...) | Extract standard errors -> mat|tib
@f p_values (x, ...) | Extract p-values -> mat|tib
@f fitted_hrf (x, sample_at, ...) | Compute fitted HRF at timepoints -> vec|mat

---

## 11. Event & Term Utilities
@f event_term (evlist, onsets, blockids, durations?=0, subset?=NULL) | Construct event term from variables -> obj
@f event_factor (fac, name, onsets, blockids?=1, durations?=0, subset?=NULL) | Categorical event sequence -> obj
@f event_variable (vec, name, onsets, blockids?=1, durations?=0, subset?=NULL) | Continuous event sequence -> obj
@f event_matrix (mat, name, onsets, blockids?=1, durations?=0, subset?=NULL) | Matrix event sequence -> obj
@f event_basis (basis, name?=NULL, onsets, blockids?=1, durations?=0, subset?=NULL) | Basis-modulated event sequence -> obj
@f elements (x, ...) | Extract elements/levels of event or term -> lst|vec
@f levels (x, ...) | Extract levels/columns of event or term -> vec
@f columns (x, ...) | Extract column labels of term -> vec
@f cells (x, ...) | Extract experimental cells from