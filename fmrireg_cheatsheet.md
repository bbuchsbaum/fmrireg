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
- Use exported HRF objects (e.g., `HRF_SPMG1`, `HRF_GAMMA`, `HRF_GAUSSIAN`, `HRF_SPMG2`, `HRF_SPMG3`, `HRF_BSPLINE`) as the `basis` argument in `hrf()` or directly in model formulas.
- Use exported basis constructors (e.g., `Poly`, `BSpline`, `Scale`, `RobustScale`, `Standardized`, `ScaleWithin`, `Ident`) to create basis objects for continuous regressors.
- Use contrast constructors (e.g., `contrast`, `pair_contrast`, `unit_contrast`, `oneway_contrast`, `interaction_contrast`, `column_contrast`, `poly_contrast`, `contrast_set`, `pairwise_contrasts`, `one_against_all_contrast`) to define contrasts for model terms.
- Use AFNI HRF objects (e.g., `AFNI_HRF`, `AFNI_SPMG1`, `AFNI_SPMG2`, `AFNI_SPMG3`, `AFNI_BLOCK`, `AFNI_dmBLOCK`, `AFNI_TENT`, `AFNI_CSPLIN`, `AFNI_POLY`, `AFNI_SIN`, `AFNI_GAM`, `AFNI_WAV`) for AFNI/3dDeconvolve workflows.

@d HRF_SPMG1 | SPM canonical HRF (single basis) -> fn
@d HRF_SPMG2 | SPM canonical HRF + temporal derivative (2 basis) -> fn
@d HRF_SPMG3 | SPM canonical HRF + temporal & dispersion derivatives (3 basis) -> fn
@d HRF_GAMMA | Gamma function HRF -> fn
@d HRF_GAUSSIAN | Gaussian HRF -> fn
@d HRF_BSPLINE | B-spline basis HRF -> fn
@d HRF_FIR | FIR basis HRF -> fn
@d AFNI_HRF | AFNI HRF object constructor -> fn
@d AFNI_SPMG1 | AFNI SPMG1 HRF -> fn
@d AFNI_SPMG2 | AFNI SPMG2 HRF -> fn
@d AFNI_SPMG3 | AFNI SPMG3 HRF -> fn
@d AFNI_BLOCK | AFNI BLOCK HRF -> fn
@d AFNI_dmBLOCK | AFNI dmBLOCK HRF -> fn
@d AFNI_TENT | AFNI TENT HRF -> fn
@d AFNI_CSPLIN | AFNI CSPLIN HRF -> fn
@d AFNI_POLY | AFNI POLY HRF -> fn
@d AFNI_SIN | AFNI SIN HRF -> fn
@d AFNI_GAM | AFNI GAM HRF -> fn
@d AFNI_WAV | AFNI WAV HRF -> fn

@d fmri_benchmark_datasets | List of simulated fMRI benchmark datasets -> lst

@f Poly (x, degree) | Polynomial basis constructor -> obj
@f BSpline (x, degree) | B-spline basis constructor -> obj
@f Scale (x) | Z-score basis constructor -> obj
@f RobustScale (x) | Robust (median/MAD) scaling basis -> obj
@f Standardized (x) | Standardized basis (mean/sd) -> obj
@f ScaleWithin (x, g) | Z-score within groups -> obj
@f Ident (...) | Identity basis for raw variables -> obj

@f contrast (form, name, where?=NULL) | General linear contrast spec -> obj
@f pair_contrast (A, B, name, where?=NULL) | Pairwise contrast spec -> obj
@f unit_contrast (A, name, where?=NULL) | Unit (sum-to-1) contrast spec -> obj
@f oneway_contrast (A, name, where?=NULL) | One-way F-contrast spec -> obj
@f interaction_contrast (A, name, where?=NULL) | Interaction F-contrast spec -> obj
@f column_contrast (pattern_A, pattern_B?=NULL, name, where?=NULL) | Column-based contrast spec -> obj
@f poly_contrast (A, name, where?=NULL, degree=1, value_map?=NULL) | Polynomial contrast spec -> obj
@f contrast_set (...) | Set of contrast specs -> lst
@f pairwise_contrasts (levels, facname, where?=NULL, name_prefix?="con") | All pairwise contrasts for levels -> obj
@f one_against_all_contrast (levels, facname, where?=NULL) | One-vs-all contrasts for levels -> obj

---

## 2. Model Construction & Data Import
### Usage:
- Use `sampling_frame()` to define scan timing and block structure.
- Use `fmri_dataset`, `matrix_dataset`, or `latent_dataset` to create datasets.
- Use `event_model()` to construct event models from formulas or lists.
- Use `baseline_model()` to define baseline/drift/nuisance models.
- Combine event and baseline models with `fmri_model()`.
- Use `read_fmri_config()` to import AFNI/fmrireg config files.
- Use `load_benchmark_dataset()`, `list_benchmark_datasets()`, `get_benchmark_summary()`, `create_design_matrix_from_benchmark()` for benchmark/simulated data.

@f sampling_frame (blocklens, TR, start_time?=TR/2, precision?=0.1) | Define scan timing and blocks -> obj
@f fmri_dataset (scans, mask, TR, run_length, event_table?=df, base_path?=".", censor?=NULL, preload?=FALSE, mode?="normal") | Create fMRI dataset from files -> obj
@f matrix_dataset (datamat, TR, run_length, event_table?=df) | Create dataset from matrix -> obj
@f latent_dataset (lvec, TR, run_length, event_table?=df) | Create dataset from latent variables -> obj

@f event_model (formula_or_list, data, block, sampling_frame, durations?=0, drop_empty?=TRUE, precision?=0.3, parallel?=FALSE, progress?=FALSE) | Build event model -> obj
@f baseline_model (basis?="constant", degree?=1, sframe, intercept?="runwise", nuisance_list?=NULL) | Build baseline/drift/nuisance model -> obj
@f fmri_model (event_model, baseline_model) | Combine event and baseline models -> obj

@f read_fmri_config (file_name, base_path?=NULL) | Read fmrireg/AFNI config file -> obj

@f load_benchmark_dataset (dataset_name?="BM_Canonical_HighSNR") | Load benchmark/simulated dataset -> lst
@f list_benchmark_datasets () | List available benchmark datasets -> df
@f get_benchmark_summary (dataset_name) | Get summary of benchmark dataset -> lst
@f create_design_matrix_from_benchmark (dataset_name, hrf, include_intercept?=TRUE) | Create design matrix from benchmark -> mat

---

## 3. HRF & Basis Utilities
### Usage:
- Use `hrf()` in model formulas to specify event convolution (basis can be string or HRF object).
- Use `gen_hrf()`, `gen_hrf_blocked()`, `gen_hrf_lagged()` for custom HRFs.
- Use `hrf_set()`, `hrf_library()` to build HRF basis sets.
- Use `as_hrf()` to wrap custom functions as HRF objects.
- Use `lag_hrf()`, `block_hrf()`, `normalise_hrf()` to decorate HRFs.

@f hrf (..., basis?="spmg1", onsets?=NULL, durations?=NULL, prefix?=NULL, subset?=NULL, precision?=0.3, nbasis?=1, contrasts?=NULL, id?=NULL, name?=NULL, lag?=0, summate?=TRUE) | HRF term spec for formulas -> obj
@f gen_hrf (hrf, lag?=0, width?=0, precision?=0.1, summate?=TRUE, normalize?=FALSE, name?=NULL, span?=NULL, ...) | Compose HRF with lag/block/normalization -> fn
@f gen_hrf_blocked|hrf_blocked (hrf?=hrf_gaussian, width?=5, precision?=0.1, half_life?=Inf, summate?=TRUE, normalize?=FALSE, ...) | Blocked HRF decorator -> fn
@f gen_hrf_lagged|hrf_lagged (hrf, lag?=2, normalize?=FALSE, ...) | Lagged HRF decorator -> fn
@f hrf_set (..., name?="hrf_set") | Combine HRFs into basis set -> fn
@f hrf_library (fun, pgrid, ...) | Build HRF library from param grid -> fn
@f as_hrf (f, name?=deparse(substitute(f)), nbasis?=1L, span?=24, params?=list()) | Wrap function as HRF object -> fn
@f lag_hrf (hrf, lag) | Lag HRF object -> fn
@f block_hrf (hrf, width, precision?=0.1, half_life?=Inf, summate?=TRUE, normalize?=FALSE) | Block HRF object -> fn
@f normalise_hrf (hrf) | Normalize HRF to peak 1 -> fn
@f empirical_hrf (t, y, name?="empirical_hrf") | Empirical HRF from data -> fn
@f hrf_basis_lwu (theta0, t, normalize_primary?="none") | LWU HRF Taylor basis -> mat
@f list_available_hrfs (details?=FALSE) | List available HRF types -> df

---

## 4. Event & Regressor Constructors
### Usage:
- Use `event_factor`, `event_variable`, `event_matrix`, `event_basis` to create event objects for event_term.
- Use `event_term()` to build event terms from lists of variables.
- Use `regressor()` or `Reg()` to create regressor objects for time series evaluation.
- Use `single_trial_regressor()` for trialwise modeling.

@f event_factor (fac, name, onsets, blockids?=1, durations?=0, subset?=NULL) | Categorical event sequence -> obj
@f event_variable (vec, name, onsets, blockids?=1, durations?=0, subset?=NULL) | Continuous event sequence -> obj
@f event_matrix (mat, name, onsets, blockids?=1, durations?=0, subset?=NULL) | Matrix event sequence -> obj
@f event_basis (basis, name?=NULL, onsets, blockids?=1, durations?=0, subset?=NULL) | Basis-modulated event sequence -> obj
@f event_term (evlist, onsets, blockids, durations?=0, subset?=NULL) | Build event term from variables -> obj
@f regressor|Reg (onsets, hrf?=HRF_SPMG1, duration?=0, amplitude?=1, span?=40, summate?=TRUE) | Create regressor object -> obj
@f single_trial_regressor (onsets, hrf?=HRF_SPMG1, duration?=0, amplitude?=1, span?=24) | Single-trial regressor -> obj

---

## 5. Model Fitting & Estimation
### Usage:
- Use `fmri_lm()` or `fmri_rlm()` to fit models to datasets.
- Use `estimate_betas()` for beta estimation (supports "mixed", "lss", "pls", "ols", etc.).
- Use `estimate_hrf()` for HRF estimation.
- Use `fmri_latent_lm()` for latent variable regression.

@f fmri_lm (formula, block, baseline_model?=NULL, dataset, durations?=0, drop_empty?=TRUE, robust?=FALSE, strategy?="runwise", nchunks?=10, use_fast_path?=FALSE, progress?=FALSE, ...) | Fit fMRI linear model -> obj
@f fmri_rlm (formula, block, baseline_model?=NULL, dataset, durations?=0, drop_empty?=TRUE, strategy?="runwise", nchunks?=10, ...) | Fit robust fMRI linear model -> obj
@f estimate_betas (x, fixed?=NULL, ran, block, method?="mixed", basemod?=NULL, ...) | Estimate betas (GLM, mixed, lss, pls, ols) -> lst
@f estimate_hrf (form, fixed?=NULL, block, dataset, bs?="tp", rsam?=seq(0,20,by=1), basemod?=NULL, k?=8, fx?=TRUE, progress?=TRUE) | Estimate HRF via GAM -> mat
@f fmri_latent_lm (formula, block, baseline_model?=NULL, dataset, durations, drop_empty?=TRUE, robust?=FALSE, autocor?="none", bootstrap?=FALSE, nboot?=1000, ...) | Latent variable regression -> obj

---

## 6. Design Matrix, Evaluation, and Visualization
### Usage:
- Use `design_matrix()` to extract model design matrices.
- Use `evaluate()` to compute regressor values over time grids.
- Use `autoplot()` or `plot()` for visualization.
- Use `design_map()`, `correlation_map()`, `design_plot()` for heatmaps and interactive plots.

@f design_matrix (x, ...) | Extract design matrix -> mat
@f evaluate (x, grid, ...) | Evaluate regressor or HRF over grid -> mat
@f autoplot (object, grid?=NULL, precision?=0.1, method?="fft", ...) | Plot regressor object -> obj
@f plot (x, ...) | Plot model or term -> obj
@f design_map (x, block_separators?=TRUE, rotate_x_text?=TRUE, fill_midpoint?=NULL, fill_limits?=NULL, ...) | Design matrix heatmap -> obj
@f correlation_map (x, method?="pearson", half_matrix?=FALSE, absolute_limits?=TRUE, ...) | Correlation heatmap -> obj
@f design_plot (fmrimod, term_name?=NULL, longnames?=FALSE, plot_title?=NULL, x_label?="Time (s)", y_label?="Amplitude", line_size?=1, color_palette?="viridis", facet_ncol?=2, theme_custom?=theme_minimal(), legend_threshold?=30, ...) | Interactive design plot (shiny) -> obj

---

## 7. AFNI/3dDeconvolve Integration
### Usage:
- Use `afni_hrf()` and `afni_trialwise()` to specify AFNI-style HRFs in model formulas.
- Use `gen_afni_lm()` and `afni_lm()` to generate and fit AFNI models.
- Use `build_afni_stims()` to generate AFNI stimulus files.
- Use `run()` on `afni_lm_spec` objects to execute AFNI 3dDeconvolve.

@f afni_hrf (..., basis, onsets?=NULL, durations?=NULL, prefix?=NULL, subset?=NULL, nbasis?=1, contrasts?=NULL, id?=NULL, lag?=0, precision?=0.3, summate?=TRUE, start?=NULL, stop?=NULL) | AFNI HRF spec for 3dDeconvolve -> obj
@f afni_trialwise (label, basis, onsets?=NULL, durations?=0, subset?=NULL, id?=NULL, start?=0, stop?=22, precision?=0.3, summate?=TRUE) | AFNI trialwise HRF spec -> obj
@f gen_afni_lm (x, ...) | Generate AFNI linear model from config -> obj
@f afni_lm (fmri_mod, dataset, working_dir?=".", polort?=-1, jobs?=1, censor?=NULL, options?=list()) | Prepare AFNI 3dDeconvolve model -> obj
@f build_afni_stims (x, ...) | Generate AFNI stimulus files -> lst
@f run[afni_lm_spec] (x, outdir, execute?=TRUE, execfun?=system, prepend?="", ...) | Run AFNI 3dDeconvolve -> NULL

---

## 8. Simulation & Benchmarking
### Usage:
- Use `simulate_bold_signal()`, `simulate_noise_vector()`, `simulate_simple_dataset()`, `simulate_fmri_matrix()` to generate synthetic data for testing and benchmarking.
- Use `evaluate_method_performance()` to compare estimation methods on benchmarks.

@f simulate_bold_signal (ncond, hrf?=HRF_SPMG1, nreps?=12, amps?=rep(1,ncond), isi?=c(3,6), ampsd?=0, TR?=1.5) | Simulate BOLD time series for conditions -> lst
@f simulate_noise_vector (n, TR?=1.5, ar?=c(0.3), ma?=c(0.5), sd?=1, drift_freq?=1/128, drift_amplitude?=2, physio?=TRUE, seed?=NULL) | Simulate fMRI noise vector -> vec
@f simulate_simple_dataset (ncond, nreps?=12, TR?=1.5, snr?=0.5, hrf?=HRF_SPMG1, seed?=NULL) | Simulate complete fMRI dataset -> lst
@f simulate_fmri_matrix (n?=1, total_time?=240, TR?=2, hrf?=HRF_SPMG1, n_events?=10, onsets?=NULL, isi_dist?="even", isi_min?=2, isi_max?=6, isi_rate?=0.25, durations?=0, duration_sd?=0, duration_dist?="lognormal", amplitudes?=1, amplitude_sd?=0, amplitude_dist?="lognormal", single_trial?=FALSE, noise_type?="none", noise_ar?=NULL, noise_sd?=1.0, random_seed?=NULL, verbose?=FALSE, buffer?=16) | Simulate multi-column fMRI time series -> lst
@f