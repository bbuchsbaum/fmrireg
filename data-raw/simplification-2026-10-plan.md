# Simplification and legacy-retirement plan (2026-10)

How this was produced: a static reachability pass over `R/` (roots = NAMESPACE
exports, S3 registrations and `.onLoad`), then six read-only subsystem audits
(fitting pipeline; AR/robust/variance/contrasts; engines/plugins/betas; group
analysis; I/O/CLI/BIDS/plotting; repository, tests and docs). Claims marked
**verified** were reproduced in R or confirmed by reading every caller; the
rest are audit findings with file:line evidence that still need a second look
before acting.

Starting point (master `a4e3035`): 35.2k lines in `R/` (97 files), 1.7k lines
of C++, 41k lines of tests (265 files), 205 Rd pages, 1,196 tracked files /
14.9 MB.

---

## Status (2026-10-17)

- **Phase 1** and **Phase 2** are done on branch `claude/simplify-retire-cruft`
  (`R/` 35,176 -> 31,952 lines). All bugs B1-B12 are fixed with regression
  tests in `tests/testthat/test-simplify-2026-10-{A,B,C,D,E,merge}.R`, each of
  which failed on the code before its fix. Refactors were pinned field by field
  (137 fitting configurations, 26 engine fits, 24 `write_results()` layouts).
- B4 was fixed minimally (`format = "fmrilm"` now errors and names the working
  route); the conversion belongs to the group-analysis plan below.
- Skipped from Phase 2 because group analysis is moving to fmrigds: the
  `fmri_meta_gds.R` rewrite, the meta C++ kernel consolidation, the single
  contrast resolver and the `fmri_ttest_gds.R` cleanup.
- **Phase 4 item 1, step P1** (branch `claude/group-delegation-p1`): fmrigds
  pinned to e9863f6; reducer `"fmrireg:meta_reg"` registered at load; the gds
  `fmri_meta()` path rewritten as reduce -> `coef_contrast()` -> wrap (18
  configurations pinned field by field at 1e-10, identical) and fitting every
  subject-level contrast; `spatial_fdr()` wraps `fmrigds::spatial_fdr_sabha()`
  and `src/spatial_fdr.cpp` is gone; `"fmrireg:fdr:spatial"` is an alias of
  `"fdr:sabha"`. `R/` -470 lines, `src/` -347. No numeric change: no fmrireg
  path or test used fmrigds's `meta:re_reg`, so its DL-trace fix does not
  reach fmrireg. Left for P2: `fmri_ttest` delegation (its classic/Welch
  engines still use only the first subject-level contrast), the meta C++
  wrapper consolidation, the legacy stack (P3).
- **Step P2** (branch `claude/group-delegation-p2`): fmrigds required at
  `>= 0.1.0.9001`. `fmri_ttest()` normalises every input once
  (`.fmri_ttest_input()`) and its classic / Welch engines are fmrigds's
  `"ols:voxelwise"` / `"ttest:welch"` kernels on the subjects x features
  matrix; classic contrasts come from the reducer's per-feature `cov_tri`.
  `fmri_ttest_gds.R` is deleted; `fmri_meta.gds` added. 1,032 configurations
  (in-memory blocks with / without V and with missing subjects, `group_data`
  gds, plain fmrigds gds; every engine, formula, contrast, sign, mc, paired,
  mu0, voxelwise_cov, weights, combine) pinned against master: classic and
  meta results identical at 1e-10; differences only where intended (Welch SE
  and per-feature deletion, recorded engine/sign); gds and in-memory results
  agree in all 328 comparable configurations; 189 configurations that errored
  now run. `ols_t_cpp()` / `welch_t_cpp()` are deprecated (remove next
  minor, with their C++). fmrigds#33 batches `"ols:voxelwise"` over complete
  features (4x faster) so the classic engine stays fast on whole-brain data.
  Left: meta C++ wrapper consolidation; `ols_t_vcov_cpp` / `fmri_ols_fit`
  voxelwise covariates and `fmri_meta_fit_extended` stay in fmrireg until
  fmrigds has a per-feature covariate reducer; the legacy stack (P3).
- **Meta kernel consolidation** (branch `claude/meta-kernel-consolidation`):
  `src/meta_kernels.cpp` has one export, `meta_fit_cov_cpp()` (optional
  voxelwise covariate); `meta_fit_cpp`, `meta_fit_contrasts_cpp` and the C++
  `meta_fit_vcov_cpp` are gone (`src/` -300 lines). The exported R functions
  keep their names and outputs; 67 configurations (four methods x robust x
  threads, missing / non-positive variances, contrasts, voxelwise, relative
  weights) identical at 1e-12.
- **Step P3, routing** (branch `claude/group-delegation-p3`, with fmrigds#34
  adding t/df NIfTI input): `group_data(format = "nifti")` reads var and
  t-only maps through fmrigds, so no non-deprecated path builds a legacy
  `group_data_nifti`/`_h5`/`_csv` object any more. 57 configurations (beta+var
  with and without mask, t-only with scalar / per-subject df; fmri_meta
  methods, weights, robust, contrasts; fmri_ttest engines; three
  t-combinations) equal the legacy reader at 1e-8. Left for the next minor:
  delete the deprecated `group_data_from_*()` constructors, the legacy
  readers and `fmri_meta.group_data_{csv,h5,nifti}` methods, the five
  legacy-only accessors, and `ols_t_cpp()` / `welch_t_cpp()` with their C++.
  The internal demo `.demo_group_data_csv()` (doc_example_helpers.R) still
  builds a legacy csv object and must move to `group_data(format = "csv")`
  first, along with the tests that construct legacy objects by hand.
- Kept on purpose: `as.array.NeuroVec` (fmristore's `LatentNeuroVec` needs it),
  the `.estimate_ar` fallback (live when fmriAR fails, now with a warning),
  plotting's `.ar_whiten` (differs from `ar_whiten_transform` on the first p
  rows), `ar_order` in the fast contrast helpers (Phase 3 deprecation).
- Behaviour changes worth knowing: a zero standard error gives t = 0, p = 1 on
  every path; failed AR estimation warns; rrr engine AR failures now also fall
  back (with the warning) like lowrank; reference-backend joint AR equals the
  matrix backend; `print()` output for `fmri_lm` and `fmri_model` changed.

---

## Phase 1 — done on branch `claude/simplify-retire-cruft`

No exported function, S3 registration or fitted result changes. Full test
suite and `R CMD check` (code, Rd, examples) pass.

| | before | after |
|---|---|---|
| tracked files / size | 1,196 / 14.9 MB | 1,070 / 4.3 MB |
| `inst/` shipped in the installed package | 2.08 MB | 0.08 MB |
| lines in `R/` | 35,176 | 34,019 |
| Imports + LinkingTo + Suggests | 33 + 4 + 26 | 31 + 3 + 19 |

1. **Repository cruft** — `outtakes/`, `outtakess/`, `experimental/`,
   `scripts/`, `test_data/`, unused `inst/extdata` (images_study, auditory,
   design_aug.txt, config.R), `data-raw/hyper_roi_1006.RDS` (6.6 MB) and its
   reader, design notes for shipped refactors, eco-atlas leftovers,
   `.DS_Store`/`.vscode`/workspace file, committed benchmark outputs.
2. **Ignore files** — `.gitignore` anchored (it ignored
   `tests/testthat/test_convenience_functions.R`, `data-raw/*.md`, testthat
   `_snaps/*.md` and the tracked `vignettes/albers-header.html`).
   `.Rbuildignore` lost ~20 dead patterns, and `^.*\.md$` became
   `^CLAUDE\.md$` so NEWS.md/README.md ship again (`utils::news()` was empty).
3. **Duplicate definitions (verified)** — `create_fmri_model` ×3 (the copy in
   `fmri_model.R` called `event_model(x =)` and could not have run; the
   `fmri_model_utils.R` copy built a poly baseline; `fmrilm.R` won by
   collation), `get_formula.fmri_model` ×2, `reshape_coef` ×2.
4. **Unreachable internals (verified)** — bootstrap family, `compute_contrast`,
   both sandwich helpers and their doc page, `fast_rlm_run`,
   `multiresponse_rlm` (robust always routes to the integrated solver, so
   `robustbase::lmrob` never ran), `ar_whiten_inplace`, the `fast_preproject`
   alias, `pull_stat_revised`, stubs that only `stop()` (`prediction_matrix`,
   `design_env`, `lss_compute_r`), and eight orphan helpers.
5. **Dead C++** — `src/inst_correlation.cpp` (never called from R).
6. **Dependencies** — `robustbase`, `progress` (only dead code), `roptim`
   (LinkingTo, never included), 7 unused Suggests.

Note for review: two of the removed pieces had been fixed recently even though
nothing called them — the sandwich helpers (2026-08-07) and the
instantaneous-correlation kernel (review C9). Restore from history if either
was meant to be wired up.

---

## Bugs found by the audit (fix independently of the cleanup)

| # | Bug | Status | Where |
|---|---|---|---|
| B1 | Every gds `fmri_meta()`/`fmri_ttest()` call re-registers fmrigds's own reducers with `force = TRUE`. `fmrigds::get_reducer("meta:fe")$provides` changes from `beta_g,var_g,se_g,z_g,p_g,Q,I2,n_eff` to `beta,se,z,p,I2,Q,df`; `meta:re_reg` switches DL→PM; fmrigds's Simes `fdr:spatial` is replaced. Plain `fmrigds::reduce()` gives different results once fmrireg has run, and no fmrireg code path consumes these reducers. | verified | `R/fmrigds_registration.R:178-252`, `fmrigds_registration_helper.R` |
| B2 | `print.fmri_lm` is defined twice; the older crayon version in `fmrilm.R` wins by collation and the newer cli version in `fmri_lm_methods.R:400` (scope, df, config) never runs. | verified | `fmrilm.R` (end), `fmri_lm_methods.R:400` |
| B3 | `print.fmri_model` reads fields fmridesign objects do not have (`event_table`, `drift_term`, `block_term`, `nuisance_term`): "Events:" prints blank and the baseline section prints only headers. | verified (output seen in tests) | `fmri_model.R` print method |
| B4 | `group_data(list(fit1, fit2), format = "fmrilm")` fails with "invalid 'file' argument"; tests only check construction of the deprecated constructor. | verified | `group_data.R`, `group_data_gds.R` |
| B5 | `parse_contrast_formula` returns all-zero weights, silently, for anything but `~ a - b` with known names. | verified | `fmri_meta_methods.R` |
| B6 | `with_package(c("shiny","plotly","bslib"))` only checks the first package (`requireNamespace` ignores the rest). | verified | `all_generic.R:4`, `design_plot.R` |
| B7 | rrr engine passes raw `cfg$ar$censor` to AR estimation instead of `resolve_censor()`: `"auto"`/`"none"` strings are not resolved and a dataset censor column is ignored. No test covers rrr + censoring. | code reading | `fmri_lm_engine_rrr.R:439` |
| B8 | `standard_error(<fmri_latent_lm>, recon = TRUE)` recomputes from `result$residuals`/`result$qr`, which the latent AR path never stores. | audit | `fmri_latent_lm.R:296-336` |
| B9 | `bids_export` with `format = c("h5","nifti")` computes volumes twice and writes each JSON sidecar twice to the same path; `overwrite`/`h5_filepath` are unused. | audit | `bids_export.R:839-1185` |
| B10 | `beta_stats_matrix_voxelwise` sets t = 0 below an absolute SE threshold (`eps^0.5`), unlike the other three t-from-covariance paths (NA). | audit | `fmri_lm_internal.R:316` |
| B11 | Reference chunkwise AR estimates phi per voxel chunk from ≤64 voxels; the matrix backend estimates per run from all voxels, so reference results depend on `nchunks`. | audit | `.iterative_ar_gls_via_fmriAR` |
| B12 | Both red CI jobs are test bugs: `test-fmridesign-generics.R:27-30` greps deparsed `.onLoad` for `:::` (covr inserts it); `test-review-E-bids.R:48` compares paths with `identical()` (fails on Windows). | audit | tests |

---

## Phase 2 — internal consolidation (no API change)

Ordered by value × safety. Each needs the usual regression test against a
closed form or a reference implementation (CLAUDE.md standard).

### Fitting pipeline
- **One AR-order helper.** `switch(cfg$ar$struct, ar1=1, ar2=2, arp=p, iid=0)`
  appears ~12 times (strategies 295/545/630, runwise 140/314/376, chunkwise
  500, fmrilm 1310, rrr 397/402, lowrank 325, adapter 150/429/760); use
  `.target_ar_order(cfg$ar)` / `cfg$ar$p %||% 0L`. ~60 LOC.
- **Contrast-weight split** duplicated 4× (runwise 53-72, chunkwise 151-170 and
  448-468, `fmri_latent_lm.R:217-227`); `prepare_fmri_lm_contrasts()` already
  returns `$simple`/`$f`. ~45 LOC.
- **Per-run design build** (`term_matrices` → `list2env` →
  `model.matrix(get_formula)`) ×6 → one `.run_design(model, run)`. Also merge
  the global-AR and global-robust-scale OLS loops in `fmri_lm_fit` (one data
  pass instead of two). ~70 LOC.
- **`runwise_lm_voxelwise`**: parallel and sequential bodies are copies
  (`fmri_lm_runwise.R:428-489`, `506-563`) → one `fit_voxel(v)`. ~55 LOC.
- **Integrated solver shims** (`fmri_lm_integrated_solver.R`): dead
  `config$method`, `config$ar_options`, character `config$robust`, struct
  `"none"`, `theta`/`phi` aliases; the `ar_options` translation is built twice
  (59-73, 178-192) only to be mapped back by `.normalize_ar_options`;
  `compute_standard_errors`, `fitted`, `residuals`, `effective_df` outputs are
  unread (drops two n×V matrices per chunk). With them goes the second
  effective-df implementation (`compute_ar_effective_df`,
  `.compute_ar_effective_df_compat`, `.get_ar_order`), which contradicts the
  "AR order does not reduce df" convention. ~225 LOC.
- **Chunkwise**: `exists("precomp")` used as a flag 8× (353-415); unused `Vu`
  and a design rebuilt twice in the fast path; `ar_coef_list` always empty;
  `covariance_model_basis` always equals `cov.unscaled`. ~45 LOC.
- **Legacy front-ends**: `.fmri_lm_formula_legacy` and
  `.fmri_lm_model_legacy` are ~90% identical (884-996 / 999-1101); the
  "engines default to joint" rule is written 5×; `fmri_lm_fit` takes both
  `compute` and the loose `strategy/nchunks/use_fast_path/progress/parallel_*`;
  `engine_ar_options`/`engine_robust_options` are always NULL. ~150 LOC.
- **AR compat layer**: move `estimate_ar_parameters` and `ar_whiten_transform`
  into `fmriAR_adapter.R` and delete `compat-legacy-ar-impl.R`; make
  `ar_whiten_transform` call the adapter (it hard-codes `parallel = TRUE` and
  ignores `fmrireg.ar.parallel_whiten`); drop the test-only options
  `fmrireg.ar.fixed_order_legacy_fallback` / `nonvoxel_legacy_fallback` and the
  `fmrireg.debug.ar` blocks; `.shared_ar_correction_design` is the identity.
  Optional: drop the `stats::ar` Yule-Walker fallback (`.estimate_ar`), the
  only AR estimator that bypasses fmriAR. ~150 LOC.
- **Reference backend scope** (medium risk): limit it to OLS (+ runwise AR) as
  a parity oracle and error otherwise, as runwise already does for robust;
  removes a third AR/robust solver stack (B11). ~60 LOC.

### Inference
- Unify the four t-from-per-voxel-covariance implementations on `.safe_tstat`
  (fixes B10). ~60 LOC.
- `process_t_contrasts` / `process_f_contrasts` (con_stats 645-727) are
  near-copies; `compute_lm_contrasts_from_suffstats` repeats the second half of
  `compute_lm_contrasts`; `getFromNamespace("qr.lm")` → `lmfit$qr`. ~75 LOC.
- `ar_order` is threaded through ~10 inference functions but only clamps
  df ≥ 1; `.voxel_inference_df` duplicates `.effective_df_from_rank`. ~40 LOC
  (the exported argument is a Phase 3 deprecation).

### Engines and plugin API
- One `.new_external_fmri_lm()` result packager: the ~45-line block is repeated
  in `fit_glm_on_transformed_series`, `fit_glm_with_config`,
  `fit_glm_from_suffstats` (`aaa_plugin_api.R:627-930`), the rrr and lowrank
  engines and `.latent_ar_fit`. ~150 LOC.
- rrr `.rrr_whiten_data` re-implements shared-AR + whitening that the lowrank
  engine already wraps (and has B7). ~60 LOC.
- Engine plumbing: duplicate preflights (rrr runs its own twice), signature
  shims, hard-coded builtin aliases, unread capability keys
  (`rank_constrained`, `whitening`, `se_modes`, `contrast_scope`). ~50 LOC.
- `srht_apply`/`srht_adjoint` and `fmri_meta_fit` use hand-written `.Call`
  strings; call the generated wrappers.

### C++
- `meta_fit_cpp`, `meta_fit_contrasts_cpp`, `meta_fit_cov_cpp`,
  `meta_fit_vcov_cpp` wrap the same `fit_one` loop (~70 identical lines each,
  `meta_kernels.cpp:268-623`); keep `meta_fit_cov_cpp` (+ optional voxelwise
  covariate) and derive the rest in R. `inv_sympd_safe` is copied in
  `meta_kernels.cpp:30` and `ols_t.cpp:10`. ~240 LOC.
- `RcppParallel` (Imports + LinkingTo) is used by no C++ in this package and
  by none of fmriAR/fmrilss/fmrihrf; `.onLoad` only calls
  `setThreadOptions()`, which changes RcppParallel's pool for *other*
  packages. Decide what `fmrireg.num_threads`/`FMRIREG_NUM_THREADS` should
  control (the OpenMP `n_threads` already reads the option), update
  vignette `a_09`, then drop the dependency.

### Group analysis
- `fmri_meta_gds.R` (551 lines): with beta+se present the code always fits
  directly, so the reducer branch, `coef:` extraction, `compute_manual_fit`,
  the `cov:`/`con_est` assay extraction and ~100 lines of dimnames patching are
  dead or fallback-only; `contrasts=` triggers a second fit (P3-2);
  `first_contrast_matrix` is defined twice. Rewrite as materialise → fit once
  → contrasts from `cov_tri` (~150 lines). ~400 LOC.
- `fmri_ttest_gds.R`: the final plan-building branch (203-266) is unreachable
  (`engine` is `match.arg`'d); `weights_gds` unused; dispatch via `get0`.
  ~95 LOC.
- Five contrast-weight builders → one resolver that errors on unknown names
  and unparseable formulas (fixes B5). ~70 LOC.
- Stop overriding fmrigds reducers (fixes B1); keep fmrireg's spatial-FDR
  posthoc under a namespaced name. ~180 LOC.

### I/O and glue
- `bids_export.R`: three h5/NIfTI writer pairs differ only in the final write
  call; mask/space validation ×6; naming loops repeated in
  `.predict_output_files`; `.extract_contrast_type` guesses from names
  although every builder sets `type`. Reuse the pattern of
  `.write_fmri_meta_map` (plan list → one writer). Fixes B9; file names must
  not change. ~300 LOC.
- `print.fmri_model` → header + fmridesign's print methods (fixes B3);
  keep the cli `print.fmri_lm` (fixes B2).
- `all_generic.R`: ~155 of 492 lines are "now imported" comments and blanks.
- `%||%` defined 3× (keep `utils_operators.R`); `as.array.NeuroVec` S3 method
  for a neuroim2 class that already has S4 methods.
- ~40 unused `importFrom` lines (17 of 18 ggplot2, 16 fmrihrf, fmridataset
  `as_sampling_frame`/`temporal_schema`, fmrilss `lss`/`mixed_solve`, dplyr
  `starts_with`, duplicate assertthat import); stale `globalVariables` in
  `imports.R:18-39` duplicate roxyglobals' `globals.R`. Verify with
  `R CMD check` ("no visible global function").
- Dead `requireNamespace()` guards on hard Imports (crayon, fmrigds ×5,
  future.apply ×2, fmrilss, jsonlite, tibble).

---

## Phase 3 — API retirements (deprecate now, remove next minor)

| Item | Why | Cost |
|---|---|---|
| Legacy `fmri_lm(...)` arguments (`strategy`, `cor_struct`, `robust`, `nchunks`, `use_fast_path`, `ar_options`, ...) | Accepted with **no** deprecation warning although `control`/`compute` replaced them; the adapter is documented as "deleted as one unit after the deprecation window". Removal deletes ~500 LOC (front-ends, `*_from_legacy`, integrated-solver shims). | Start warning now; ~100 test files and several roxygen examples (`fmri_lm` uses `strategy="chunkwise", nchunks=1`) need migrating first. |
| Exported `.fmri_lm_formula_legacy`, `.resample_param` | dot-names exported by accident | none known |
| `estimate()` | exported stub that only errors | trivial |
| `fmrireg_cli`, `install_cli`, `exec/fmrireg`, `inst/cli` (~635 LOC) | the CLI only wraps `help` and three benchmark lookups that are already exported functions | two exports |
| `design_plot()` (Shiny/plotly) | only user of plotly + bslib (hard Imports); broken example; deprecated `geom_line(size=)`/`aes_string`; B6 | drop 2 Imports |
| Raw C++ exports `ols_t_cpp`, `welch_t_cpp`, `ols_t_vcov_cpp`, `meta_fit_vcov_cpp`, `mixed_solve_cpp` | low-level kernels; `mixed_solve_cpp` is a misnamed pass-through to `fmrilss::mixed_solve` | un-export |
| `soft_subspace_options` | duplicates `projection_spec`, no internal caller | un-export |
| `lowrank_control(time_sketch = list(method = "ihs"))` | deprecated in 0.3.0, runs exact OLS | ~80 LOC |
| `fmri_latent_lm(bootstrap=, nboot=)`; `lowrank_control(ncomp=, noise_pcs=)`; `glm_lss(use_cpp=)` (deprecated since 0.1.1); unused `estimate_betas` args (`maxit`, `fracs`, `block`); `return_cov`, `chunk_size`, `verbose` in group fits; `fmri_ttest(mask=)` | accepted but never read | deprecate |
| `ar_order` in `compute_lm_contrasts(_from_suffstats)` | only clamps df ≥ 1 | deprecate |
| Engine-fit top-level `betas_fixed`, `sigma2`, `vcov_inv` | duplicate `result` fields; vignette `sketched-ar` reads `betas_fixed` | point to `coef()` |
| `estimate_hrf(bs=, fx=)` | deprecated arguments | bump fmrihrf pin first (its new `hrf_bspline` replaces fmrireg's private basis, ~60 LOC) |
| Imports `crayon`, `furrr`, `purrr`, `magrittr`, `tidyr`, `plotly`, `bslib` | each used in 1-2 places (or only by `design_plot`); `crayon` → `cli`, `furrr` → `future.apply` or `fmrilss::mixed_solve(Y = matrix)`, `purrr` → base, `%>%` → `|>` | 7 Imports |

---

## Phase 4 — structural

1. **Group analysis moves to fmrigds (proposed division of labour).**
   fmrigds owns the group data model, I/O and every across-subject estimator
   as a registered reducer (FE/DL, regression, OLS t, permutation, z/p
   combination, post-hoc FDR); fmrireg owns subject-level fitting, turning
   `fmri_lm` fits into group inputs (`as_group_data(fits, contrast =)`), thin
   `fmri_meta()`/`fmri_ttest()`/`spatial_fdr()` front-ends with fmrireg's
   argument names and classes, run pooling, and BIDS/image output.
   Estimator correctness tests live with the estimator. fmrireg never
   overrides a reducer name; anything it registers is namespaced.
   `cov_tri` (coefficient covariance for contrasts): fmrigds already
   transports a reducer's `cov_tri` into `cov:<a>:<b>` assays
   (`coef_cov_tri()`); the cost is the estimator, not the covariance. Now:
   register `fmrireg:meta_reg` (PM/DL/FE, Huber, weights) wrapping the
   existing `meta_fit_cov_cpp` (~80 LOC), which also fixes fit-only-the-first-
   contrast. Later, only if another consumer needs PM/Huber: port `fit_one`
   to fmrigds and delete fmrireg's kernel. Steps: (a) fmrigds: fix the bogus
   `cov_tri` assay created when a reducer lists it in `provides`, add
   `coef_contrast()`, extra per-voxel assays for covariates, `fdr:sabha`
   (fmrireg's empirical-null spatial FDR) and `ttest:welch`; (b) fmrireg P1:
   bump the fmrigds pin (fixes upstream DL-regression tau2; `a_09` uses the
   now-unexported `plan()`), register `fmrireg:meta_reg` at load, rewrite the
   gds `fmri_meta` path (~-600 LOC); (c) P2: `fmri_ttest` delegates, delete
   `fmri_ttest_gds.R`, `ols_t.cpp` (except vcov) and the moved FDR code
   (~-1,000 LOC); (d) P3, next minor: remove the legacy csv/h5/nifti stack and
   accessors (needs an fmrigds t/df NIfTI adapter) (~-2,500 LOC).
   The earlier, fmrireg-internal version of this item follows for reference.
   **One group-analysis stack.** gds is the future (`group_data()` always
   builds gds; the `group_data_from_*` constructors are `.Deprecated`). The
   legacy stack is still live only for var/t-only NIfTI. Plan:
   (a) Phase 2 items above; add `.fmri_meta_from_matrices(Y, V | T, DF, X)`
   used by both stacks. (b) `.as_group_gds(x)` converter at the boundary;
   legacy `fmri_meta` methods become shims; route var/t NIfTI through it
   (≈ −700 LOC); give non-voxel fits class `fmri_meta_roi`. (c) One
   `fmri_ttest` input normaliser; delete `fmri_ttest_gds.R` (≈ −300 LOC) and
   make `paired`/`mu0`/`voxelwise_cov` work on gds input (today they only work
   on hand-built `list(blocks=)` objects, which the main correctness tests
   use). (d) Next minor: remove deprecated constructors, legacy readers
   (~1,000 LOC of `group_data_{csv,h5,nifti}.R`) and the five legacy-only
   accessors (`extract_csv_data`, `get_rois`, `get_contrasts`,
   `read_h5_full`, `read_nifti_full`).
2. **Benchmark data on demand.** `data/fmri_benchmark_datasets.rda` (853 KB)
   is deterministic and small to generate; freezing it already forced one
   regeneration and a ~200-line HRF recipe/version layer. Generate on first use
   with a fixed seed and cache; drop `data/`, the recipe layer and the
   data-raw generators (~850 KB, ~850 lines). Or retire the feature with the
   CLI.
3. **fmri_latent_lm** → let `fmri_lm` tag latent fits itself and compute the
   reconstruction from stored `cov.unscaled` (fixes B8).

---

## Legacy group-data deletion checklist (next minor release)

Everything below is marked in the code with `LEGACY-GROUP-DELETE (next
minor)`; `grep -rn "LEGACY-GROUP-DELETE" R src tests vignettes` lists every
site, and after the deletion it must list nothing. Users have been warned
for a release: every public entry point emits a deprecation warning (the
constructors and `fmri_meta.group_data_csv/_h5` since 0.3.0; the accessors,
`fmri_meta.group_data_nifti` and legacy input to `fmri_ttest()` since the
2026-10 prep, via `.legacy_group_deprecated()`). The list was derived from the
call graph: functions reachable only from the deprecated entry points.

**Delete whole files**
- `R/group_data_csv.R`: `group_data_from_csv()`, `.group_data_from_csv_impl()`,
  `validate_effect_cols()`, `wide_to_long_format()`, `extract_csv_data()`,
  `get_rois()`, `get_contrasts()`, `validate_group_data_csv()`.
- `R/group_data_h5.R`: `group_data_from_h5()`, `group_data_from_fmrilm()`,
  `read_h5_metadata()`, `detect_h5_file_type()`,
  `extract_subject_ids_from_paths()`, `read_h5_chunk()`, `.close_h5_handle()`,
  `.read_h5_subject()`, `read_h5_full()`, `validate_group_data_h5()`.
- `R/group_data_nifti.R`: `group_data_from_nifti()`,
  `.group_data_from_nifti_impl()`, `read_nifti_header()`, `read_nifti_data()`,
  `validate_nifti_dimensions()`, `read_nifti_chunk()`, `read_nifti_full()`,
  `validate_group_data_nifti()`. (`t_to_beta_se()` already moved to
  `R/metafuns.R`.)

**Delete functions in kept files**
- `R/fmri_meta.R`: `fmri_meta.group_data_h5()`, `fmri_meta.group_data_nifti()`,
  `fmri_meta.group_data_csv()`, `fit_meta_chunk()`, `fit_meta_single()`,
  `parse_meta_formula()`, `.meta_slice_weights()`.
- `R/group_data.R`: `validate_group_data()`, `print.group_data()`,
  `summary.group_data()`, `n_subjects.group_data()`,
  `get_subjects.group_data()`, `get_covariates.group_data()`,
  `.legacy_group_deprecated()`, `.legacy_group_quiet()`.
- `R/fmri_ttest_helpers.R`: `ols_t_cpp()`, `welch_t_cpp()` (and the
  `ttest-kernels-deprecated` roxygen block); `src/ols_t.cpp`: `ols_t_cpp()`,
  `welch_t_cpp()` (keep `ols_t_vcov_cpp()`, used by `fmri_ols_fit()`); rerun
  `Rcpp::compileAttributes()`.
- `R/doc_example_helpers.R`: `.demo_group_data_csv()` (live examples use
  `.demo_group_data()`).

**Delete branches in kept functions** (each marked in place)
- `R/fmri_meta.R`, `fmri_meta()`: the `validate_group_data()` call.
- `R/fmri_meta_methods.R`: `use_exact_roi` and its branch (csv contrast SEs);
  the h5/nifti tests in `coef_image.fmri_meta()`; the h5/nifti branches of
  `reconstruct_image()`, of the meta mask helper and of
  `.extract_fmri_meta_source_files()`.
- `R/fmri_ttest_helpers.R`: the nifti/h5 branches of
  `.fmri_ttest_materialize_effects()`; the no-blocks `else` branch of
  `.fmri_ttest_input()` (legacy file-backed objects). The in-memory
  `list(blocks =)` form and `paired_diff_block()` stay.
- `R/group_data_gds.R`: the conversion of `group_data_nifti` objects in
  `.group_data_gds()`.

**Behaviour to settle when deleting**
- (Done in the prep: `summary.group_data_gds()` exists, so gds objects no
  longer dispatch to the legacy `summary.group_data()`; gds fits on non-voxel
  samples are `fmri_meta_roi`; a csv table without an ROI column is one ROI,
  "overall".) The `summary.group_data` roxygen block also documents
  `summary.group_data_gds` (shared `@rdname`): keep that block when deleting
  the legacy method.
- `fmri_meta()` / `fmri_ttest()` on a legacy object will error with "no
  method"; say so in NEWS.

**Generated files**: `roxygen2::roxygenise()` removes the exports
(`group_data_from_csv/_h5/_nifti/_fmrilm`, `extract_csv_data`, `get_rois`,
`get_contrasts`, `read_h5_full`, `read_nifti_full`, `ols_t_cpp`,
`welch_t_cpp`, `S3method(fmri_meta, group_data_csv/_h5/_nifti)`, the
`group_data` print/summary/accessor methods) and their Rd files.

**Tests**: each test that exercises only the legacy stack carries the marker
on the line above its `test_that()` (or at the top of the file when the whole
file goes): delete the marked blocks and files. Every test of a live feature
that used legacy objects as input now builds gds group data instead.
- whole files: `test-coverage-h5-meta-mock.R`, `test-gds-h5-parity.R`,
  `test-gds-nifti-parity.R`, `test-review-B-group-io.R`;
- marked blocks in: `test-coverage-cli-group-soft.R`,
  `test-coverage-design-h5-spatial.R`, `test-coverage-group-nifti-io.R` (keep
  its `t_to_beta_se` block), `test-coverage-group-ttest-meta.R`,
  `test-coverage-meta-recon-lowrank-ar.R`, `test-coverage-meta-robust-tau2.R`,
  `test-coverage-model-meta-contrasts.R`,
  `test-coverage-preflight-group-buffer.R`, `test-coverage-ttest-engines-mc.R`
  (four lines inside a kept block), `test-gds-integration.R`,
  `test-gds-parity-methods.R`, `test-gds-ttest-parity.R`,
  `test-group_analysis.R`, `test-review-B-meta-weights.R`,
  `test-review-B-ttest-classic.R`, `test-simplify-2026-10-P2.R`,
  `test-simplify-2026-10-P3.R`.
- `tests/testthat/setup.R` sets `fmrireg.suppress_deprecation = TRUE`; drop it
  if nothing deprecated remains.

**Docs**: NEWS (breaking: removed functions, replacements); the `a_09`
vignette note on legacy constructors.

## Tests, docs and CI

- **Coverage-wave tests.** The 93 `test-coverage-*` files (12.4k lines, one PR)
  have 42% numeric assertions vs 89% for `test-review-*`; 28 files have none
  and they hold 123 of the suite's 130 always-true assertions (e.g.
  `expect_true(inherits(x,"error") || !is.null(x) || is.null(x))`). Compare
  covr before/after, delete the 28, fold the rest into thematic files with one
  expected outcome per check.
- **Merge clusters**: spatial FDR (3 files), effective df (3), tidy (3),
  robust (5), solver/QR (14), AR (26), CLI (4, each with its own
  `capture_cli()`), group/meta (37). Retire the `test_` prefix while merging;
  move the 50 file-local `matrix_frame` builders into `helper-frame.R`.
- 268 `skip_if_not_installed()` calls on hard Imports (9 on "fmrireg"
  itself); tests that need undeclared pryr/devtools never run; 13 files still
  use `context()`; add `Config/testthat/edition: 3`.
- **Workflows**: delete `pkgcheck` (fails structurally every push and rewrites
  issue #179) and `update-codemeta` (bot commits to master; codemeta only
  serves rOpenSci) plus `codemeta.json`; replace `validate-readme`'s grep with
  "knit README, `git diff --exit-code`"; README coverage badge points at
  `main` (default branch is `master`); `codecov.yml` filters on a flag the
  upload never sets; bump `checkout@v3`.
- **Docs**: `_pkgdown.yml` has no `reference:` index (205 topics alphabetical),
  an empty "tutorials" menu, a news link to a non-existent v0.1.0 tag, and a
  home description about HRF construction (now fmrihrf). `a_05_contrasts`
  re-teaches fmridesign (~330 of 533 lines) and `a_10_dataset` re-teaches
  fmridataset — trim and link. NEWS 0.3.0 is 660 lines; condense to user-facing
  bullets. CLAUDE.md's architecture section names files that no longer exist
  (`event-classes.R`, `hrf.R`, `fmri_dataset.R`, `contrast.R`) and the
  `*_dataset` classes now called `*_frame`.
- `make clean` runs `rm -rf man/*.Rd`, deleting the two hand-written Rd files
  (`fmrireg_cli.Rd`, `install_cli.Rd`).

---

## Looks like cruft, keep

`.demo_*` helpers (roxygen examples use them via `:::`); `.preproject_pinv`
(reached through `$.fmri_preproject`); `solve_integrated_glm` (latent path,
plugin API, reference backend); `cpp_srht_*` and `meta_fit_cpp` (string
`.Call`s); `estimate_ar_parameters`/`ar_whiten_transform` (live, just
misfiled); `.iterative_ar_gls_via_fmriAR` (ARMA); `robust_iterative_fitter`
(the single IRLS); `fmri_lm_variance.R`; `resolve_basis` (injected into
formulas as `fmrireg::resolve_basis`); the engine registry and plugin
vignette; the batch pipeline (`fmri_template` → `instantiate` → `run_jobs` →
`reduce_*` → `export_jobs`/`collect_results`), which is one coherent design;
`fmri_meta.group_data_nifti` and the NIfTI readers until the converter lands;
`combine_t_statistics`; `meta_stouffer`/`do_fixef` (run pooling, different
zero-SE semantics from the C++ kernel); `spatial_fdr(fit, p = coef)`;
`%||%` in `utils_operators.R` (base R only has it from 4.4).
