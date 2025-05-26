## Micro-DSL (v2.5) & Output Format
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
`Sections+` (H2 `1.` Title; H2 Title (unnumbered ok))
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
* Format: ` - param_name : type_abbr Brief, essential clarification.`
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
* **Pure Markdown only.** Adhere strictly to the Micro-DSL v2.5.
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

## Formal Grammar "v2.5 Micro-EBNF"
Cheatsheet ::= Header Legend? Section+ Deps?
Legend ::= H2_TOKEN TEXT_CONTENT NEWLINE_TOKEN Block
Header ::= H1_TOKEN TEXT_CONTENT NEWLINE_TOKEN+ (H2Section)*
H2Section ::= H2_TOKEN TEXT_CONTENT NEWLINE_TOKEN Block
Section ::= H2_TOKEN (NUMBER_TOKEN PERIOD_TOKEN)? TEXT_CONTENT NEWLINE_TOKEN Entry+ HR_TOKEN?
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
ParamDesc ::= IDENT_TOKEN COLON_TOKEN TYPE_ABBR_TOKEN TEXT_CONTENT?
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

## 1. Core Constructors & Data Objects

@f event_model (x, data, block, sampling_frame, drop_empty?=TRUE, durations?=0, ...) | Build event-based fMRI model -> obj
@f baseline_model (basis?="constant", degree?=1, sframe, intercept?="runwise", nuisance_list?=NULL) | Construct baseline model -> obj
@f fmri_model (event_model, baseline_model) | Combine event and baseline models -> obj
@f fmri_dataset (scans, mask, TR, run_length, event_table?=df, base_path?=".", censor?=NULL, preload?=FALSE, mode?="normal") | File-backed fMRI dataset -> obj
@f fmri_mem_dataset (scans, mask, TR, run_length, event_table?=df, base_path?=".", censor?=NULL) | In-memory fMRI dataset -> obj
@f matrix_dataset (datamat, TR, run_length, event_table?=df) | Matrix-format fMRI dataset -> obj
@f latent_dataset (lvec, TR, run_length, event_table?=df) | Latent component dataset -> obj
@f sampling_frame (blocklens, TR, start_time?=TR/2, precision?=0.1) | Define scan/block timing -> obj

---

## 2. HRF & Basis Constructors

@f HRF (fun, name, nbasis?=1, span?=24, param_names?=NULL) | Create HRF object -> fn
@f as_hrf (f, name?=chr, nbasis?=1, span?=24, params?=lst) | Wrap function as HRF -> fn
@f bind_basis (...) | Combine HRFs into basis set -> fn
@f gen_hrf (hrf, lag?=0, width?=0, precision?=0.1, summate?=TRUE, normalize?=FALSE, name?=NULL, span?=NULL, ...) | Decorate HRF -> fn
@f hrf ( ... , basis?="spmg1", onsets?=NULL, durations?=NULL, prefix?=NULL, subset?=NULL, precision?=0.3, nbasis?=1, contrasts?=NULL, id?=NULL, lag?=0, summate?=TRUE) | HRF term spec for formulas -> obj
@f block_hrf (hrf, width, precision?=0.1, half_life?=Inf, summate?=TRUE, normalize?=FALSE) | Blocked HRF decorator -> fn
@f lag_hrf (hrf, lag) | Lag HRF by seconds -> fn
@f normalise_hrf (hrf) | Normalize HRF to peak 1 -> fn
@f empirical_hrf (t, y, name?="empirical_hrf") | Empirical HRF from data -> fn
@f hrf_set (...) | Combine HRFs into basis set -> fn
@f hrf_library (fun, pgrid, ...) | HRF library from param grid -> fn
@f list_available_hrfs (details?=FALSE) | List HRF types -> df

---

## 3. Event & Term Constructors

@f event_term (evlist, onsets, blockids, durations?=0, subset?=NULL) | Build event term -> obj
@f event_factor (fac, name, onsets, blockids?=1, durations?=0, subset?=NULL) | Categorical event -> obj
@f event_variable (vec, name, onsets, blockids?=1, durations?=0, subset?=NULL) | Continuous event -> obj
@f event_matrix (mat, name, onsets, blockids?=1, durations?=0, subset?=NULL) | Matrix event -> obj
@f event_basis (basis, name?=NULL, onsets, blockids?=1, durations?=0, subset?=NULL) | Basis-modulated event -> obj

---

## 4. Regressor Constructors & Evaluation

@f regressor (onsets, hrf?=HRF_SPMG1, duration?=0, amplitude?=1, span?=40, summate?=TRUE) | Create regressor object -> obj
@f single_trial_regressor (onsets, hrf?=HRF_SPMG1, duration?=0, amplitude?=1, span?=24) | Single-trial regressor -> obj
@f null_regressor (hrf?=HRF_SPMG1, span?=24) | Empty regressor -> obj
@f evaluate[Reg] (x, grid, precision?=0.33, method?="fft", sparse?=FALSE, ...) | Evaluate regressor at grid -> vec|mat
@f shift[Reg] (x, shift_amount, ...) | Shift regressor onsets -> obj

---

## 5. Model Fitting & Estimation

@f fmri_lm (formula, block, baseline_model?=NULL, dataset, durations?=0, drop_empty?=TRUE, robust?=FALSE, strategy?="runwise", nchunks?=10, use_fast_path?=FALSE, ...) | Fit fMRI linear model -> obj
@f fmri_rlm (formula, block, baseline_model?=NULL, dataset, durations?=0, drop_empty?=TRUE, strategy?="runwise", nchunks?=10, meta_weighting?="inv_var", ...) | Robust fMRI linear model -> obj
@f fmri_latent_lm (formula, block, baseline_model?=NULL, dataset, durations, drop_empty?=TRUE, robust?=FALSE, autocor?="none", bootstrap?=FALSE, nboot?=1000, ...) | Fit model to latent dataset -> obj
@f estimate_betas[fmri_dataset|matrix_dataset|latent_dataset] (x, fixed?=NULL, ran, block, method?="mixed", basemod?=NULL, hrf_basis?=NULL, hrf_ref?=NULL, ...) | Estimate betas -> obj
@f estimate_hrf (form, fixed?=NULL, block, dataset, bs?="tp", rsam?=seq(0,20,1), basemod?=NULL, k?=8, fx?=TRUE) | Estimate HRF via GAM -> mat
@f find_best_hrf[fmri_mem_dataset|matrix_dataset] (dataset, onset_var, hrflib?=NULL, rsam?=seq(0,24,1), ncomp?=5, nsplits?=3, block?=NULL, basemod?=NULL, ...) | Find best HRF(s) -> obj

---

## 6. Design Matrix & Visualization

@f design_matrix[event_term|event_model|baseline_model|fmri_model] (x, ...) | Extract design matrix -> tib
@f design_map[event_model|baseline_model|fmri_model] (x, block_separators?=TRUE, rotate_x_text?=TRUE, fill_midpoint?=NULL, fill_limits?=NULL, ...) | Design matrix heatmap -> obj
@f plot[event_term|baseline_model|fmri_model] (x, ...) | Plot design or baseline model -> obj
@f autoplot[Reg] (object, grid?=NULL, precision?=0.1, method?="fft", ...) | Plot regressor -> obj
@f correlation_map[event_model|baseline_model|fmri_model] (x, method?="pearson", half_matrix?=FALSE, absolute_limits?=TRUE, ...) | Correlation heatmap -> obj
@f design_plot (fmrimod, term_name?=NULL, longnames?=FALSE, plot_title?=NULL, x_label?="Time (s)", y_label?="Amplitude", line_size?=1, color_palette?="viridis", facet_ncol?=2, theme_custom?=theme_minimal(), legend_threshold?=30, ...) | Interactive design plot (shiny) -> obj

---

## 7. Contrast Specification & Weights

@f contrast (form, name, where?=NULL) | Define linear contrast -> obj
@f unit_contrast (A, name, where?=NULL) | Unit contrast vs baseline -> obj
@f pair_contrast (A, B, name, where?=NULL) | Pairwise contrast -> obj
@f one_against_all_contrast (levels, facname, where?=NULL) | One-vs-all contrasts -> obj
@f oneway_contrast (A, name, where?=NULL) | One-way F-contrast -> obj
@f interaction_contrast (A, name, where?=NULL) | Interaction F-contrast -> obj
@f poly_contrast (A, name, where?=NULL, degree?=1, value_map?=NULL) | Polynomial contrast -> obj
@f column_contrast (pattern_A, pattern_B?=NULL, name, where?=NULL) | Contrast by column regex -> obj
@f contrast_set (...) | Set of contrasts -> obj
@f pairwise_contrasts (levels, facname, where?=NULL, name_prefix?="con") | All pairwise contrasts -> obj
@f contrast_weights (x, ...) | Compute contrast weights -> lst
@f Fcontrasts (x, ...) | Generate F-contrasts -> lst
@f to