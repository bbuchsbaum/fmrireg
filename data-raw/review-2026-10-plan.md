# Review 2026-10: correctness and performance fix plan

Branch: `fix/review-2026-10`. Baseline before any change: 6,897 tests in 239
files, 0 failures, 2 skips (R 4.3.3, dependencies at the commits pinned in
`DESCRIPTION`).

## Ground rules

- Every fix lands with a regression test that **fails on the baseline code and
  passes after the fix**. Tests must check numbers, not shapes: compare against
  a closed form, `stats::lm`, `sandwich`, `metafor`, or a simulation with known
  truth.
- Statistical-calibration tests use pure-noise simulations with a fixed seed
  and assert the false-positive rate sits inside a binomial tolerance band.
- Performance changes must be result-preserving: each one gets an
  equivalence test (old path vs new path, or new path vs reference) and, where
  the claim is "does not load the full dataset" or "does not build an n x n
  matrix", a test that counts calls or inspects allocations directly.
- No skipped, disabled or loosened existing tests. If an existing test encoded
  a wrong number, the change is called out in the commit message.
- Each stream commits separately; the full suite runs after merging.

## Stream A: spatial FDR (`R/spatial_fdr.R`, `src/spatial_fdr.cpp`)

- [x] A1. Empirical null: scale estimated from the central 60% of |z| is
      ~0.57 under the null; correct for truncation (or use all finite z).
      Test: 1e5 pure-noise z, default args, rejections at alpha=.05 ≈ 0;
      sigma0 within 0.95–1.05.
- [x] A2. Group-ID mismatch between R (`factor` order) and C++ (first-
      appearance order). Make one numbering authoritative; vectorise the
      weight expansion. Test: result invariant to relabelling/permuting group
      labels; neighbor smoothing uses the right groups.

## Stream B: group analysis and meta-analysis

- [x] B1. `meta_stouffer` (run-pooled F-contrasts) returns a weighted mean,
      not a Stouffer z. Use Σwz/√Σw², upper-tail `qnorm`/`pnorm`. Test:
      k runs at identical p give z = √k·qnorm(1-p); null calibration.
- [x] B2. `weights = "equal"/"custom"` treated as absolute variances. Treat as
      relative weights with an estimated residual scale. Custom weights must
      be sliced per chunk. Test: equal-weights FE equals one-sample t
      (`t.test`) for any data scale; custom-weights equals weighted LS
      (`lm(weights=)`); chunk_size does not change results.
- [x] B3. t→z conversion overflows to Inf for p < ~1e-16 (7 copies). One
      helper using log-scale tails. Test: t = 40, df = 20 gives finite z
      matching `qnorm(pt(..., log.p=TRUE), log.p=TRUE)`.
- [x] B4. Stouffer/Lancaster t-combination counts missing subjects. Test:
      result with NA subjects equals result with those subjects removed.
- [x] B5. Meta estimators verified against `metafor::rma` (FE, DL, PM) for
      intercept-only and meta-regression; label "reml" honestly (either
      implement or map with a message). Test: estimates/SE/tau2/Q within 1e-6.
- [x] B6. Contrast SE ignores coefficient covariance by default. Test against
      `metafor` contrast SE (`anova(rma, X=)`).
- [x] B7. Classic t-test engine: NA subject wipes the voxel. Drop non-finite
      subjects per voxel (as the meta engine does). Test vs `t.test` on the
      complete cases.
- [x] B8. Paired test with `voxelwise_cov` errors on rownames. Test that it
      runs and matches `lm` per voxel.
- [x] B9. HDF5 reader leaks handles / silences failures. Test handle closes
      and that a bad file is an error, not NA.
- [x] P-B1. NIfTI group reader re-reads every subject per chunk. Read once.
      Test: number of reads == number of subjects (mock the reader); result
      identical across chunk sizes.
- [x] P-B2. `schedule(dynamic)` in meta kernels; no result change (existing
      tests + B5 references).

## Stream C: betas, contrasts, tidy, design utilities

- [x] C1. `estimate_betas(method = "ols"/"mixed")` residualises Y but not X.
      Test: noiseless simulation recovers true betas exactly (vs `lm`).
- [x] C2. `tidy(type = "contrasts")` returns only the first contrast. Test:
      two t-contrasts + one F-contrast, all present with correct values.
- [x] C3. `generate_main_effect_contrast` / `generate_interaction_contrast` row
      order vs `expand.grid`. Test: contrast applied to cell means from a
      known factorial design yields the textbook main effect / interaction.
- [x] C4. Post-hoc `fit_contrasts.fmri_lm`: silent positional fallback and
      dropped specs. Map by full column names, error when unmatched, support
      every spec via `contrast_weights()`. Test: post-hoc contrasts equal the
      ones computed during the fit, for pair, oneway and unit specs, with
      several terms.
- [x] C5. Post-hoc contrasts on AR/robust chunkwise fits mix unwhitened
      `(X'X)^-1` with whitened sigma. Store the matching covariance, or refuse.
      Test: post-hoc equals in-fit contrast for an ar1 fit.
- [x] C6. `hrf_smoothing_kernel` normalisation (asymmetric, NaN). Test:
      symmetric, unit diagonal, finite with buffer_scans = 0.
- [x] C7. Absolute 1.5e-8 / sqrt(eps) thresholds zero out small-scale data.
      Test: t/F invariant to rescaling Y by 1e-4.
- [x] C8. F p-value `1 - pf()` loses precision. Use `lower.tail = FALSE`.
- [x] P-C1. LSS fixed effects and mixed path: one multi-response solve instead
      of per-voxel `lm.fit`/`lsfit`; drop identity `map_voxels`. Equivalence
      test vs previous per-voxel results.
- [x] P-C2. Reference F-contrast path: vectorise; compute `cov.unscaled` once.
      Equivalence test vs `.fast_F_contrast` and `anova(lm)`.

## Stream D: GLM fitting core

- [x] D1. `reestimate_phi = TRUE` estimates AR from whitened residuals.
      Use unwhitened residuals. Test: AR(1) phi = 0.5 simulation, robust +
      reestimate_phi recovers phi ≈ 0.5 and matches non-reestimated fit.
- [x] D2. Robust df counts zero-weight volumes; chunkwise vs runwise use
      different sigma; returned weights differ from the weights used for the
      final betas. Test: runwise == chunkwise for robust fits on 1 and 2
      runs; betas reproducible from returned weights by WLS.
- [x] D3. Runwise volume weights break on multi-run data. Test: two runs,
      full-length weight vector, matches WLS `lm(weights=)` on stacked data
      with run intercepts; zero weights reduce df.
- [x] D4. Effective-df: one implementation, the one CLAUDE.md documents; delete
      unreachable copies. Test against hand computation.
- [x] D5. Rank-deficient designs: `n <- dfres + p` wrong. Test with an aliased
      column: df matches `lm`.
- [x] P-D1. `.fast_preproject` builds an n×n identity on every call. Compute
      the pseudo-inverse lazily. Test: results identical; `Pinv` never built on
      the OLS/voxelwise-AR path (trace calls).
- [x] P-D2. Chunking loads the full dataset (several times). Lazy chunk
      descriptors, run rows from the temporal schema, no `Y_orig`. Test: count
      data-access calls on a mock dataset; results identical across
      `nchunks`.
- [x] P-D3. HAC/Satterthwaite: skip covariance when unused; hoist bread and
      `XtXinv` out of the voxel loop. Equivalence test vs `sandwich::NeweyWest`
      / `vcovHAC` for a few voxels.
- [x] P-D4. `solve_glm_core`: single `qr.qty` pass for betas and exact RSS.
      Equivalence vs `lm`.

## Stream E: engines, latent, I/O, CRAN hygiene

- [x] E1. `rrr_gls` SEs ~√(V/k) too small. Use per-voxel residual variance.
      Test: pure-noise calibration (|t|>1.96 rate ≈ 5%); full-rank result
      equals OLS SEs.
- [x] E2. `fmri_latent_lm` silently ignores `robust` / `autocor`. Read the
      real config fields; error for unsupported combinations; remove the
      undefined-function globals. Test: AR(1) latent fit differs from OLS and
      matches dense AR fit on the reconstructed data.
- [x] E3. RNG: no `set.seed` without restoring state (`rrr_gls`, simulators,
      `bids_export` staging dir). Test: `.Random.seed` unchanged after calls.
- [x] E4. `write_results(format = "nifti")` must not require fmristore.
      Test with `requireNamespace` mocked to FALSE.
- [x] E5. Benchmark data must not hard-fail on fmrihrf version mismatch.
      Test with the version mocked.
- [x] E6. BIDS sidecar DegreesOfFreedom scalar (or range), contrast-label
      collisions detected, `desc` sanitised and kept. Tests on file names and
      JSON content.
- [x] P-E1. `rrr_gls` bootstrap: hoist QR of the fixed design. Equivalence
      test with fixed seed.

## Final pass

- [x] F1. Merge streams; full test suite green.
- [x] F2. Remove dead code proven unreachable (list in commit), keeping tests
      that exercised live behaviour.
- [x] F3. `R CMD build` + `R CMD check --no-manual` with no new warnings.
- [x] F4. Update CLAUDE.md with what changed and the new reference tests.
- [x] F5. Results table below: test counts before/after, timings of the
      performance items.

## Results

### Test suite (R 4.3.3, `NOT_CRAN=true`)

| | Files | Expectations | Failures | Errors | Skips |
|---|---|---|---|---|---|
| Before (5bd3ad1) | 239 | 6,897 | 0 | 0 | 2 |
| After (all files but `test_betas.R`) | 262 | 8,089 | 0 | 0 | 2 |
| After, `test_betas.R` alone | 1 | 9 | 0 | 0 | 0 |

Every new `test-review-*` file was run against the pre-fix code and failed
there (except pure equivalence guards, which pass on both sides by design
and are paired with a call-count or allocation assertion that fails before).

`R CMD check --no-manual --no-tests --ignore-vignettes`: examples, Rd,
code/documentation, R code and compiled-code checks all OK. Remaining
1 WARNING / 2 NOTEs are environmental (missing `en_US.UTF-8` locale,
Suggests not installable offline, installed size).

### Defects found during the merge (final pass, `test-review-F-integration.R`)

- Reference backend read `result$sigma`, which partially matched `sigma2`:
  AR standard errors ~3x too large.
- Multi-run runwise fits stored one scalar sigma/rss/resvar for all voxels
  (`as.matrix()` on per-voxel vectors before `colSums`); pre-dates the review.
  Post-hoc contrast SEs on those fits were ~3x too large.
- `fmri_ttest` group contrasts: t/z flipped twice (meta, classic) or the
  estimate not flipped (Welch, voxelwise covariate, GDS); 11 failing checks.
- Soft subspace projection: residual df now tr((I-H)P^2) on runwise iid fits.
- `fit_contrasts()` refuses multi-run runwise fits that whiten/weight per run.

### Performance (stream reports; elapsed seconds, noisy shared machine)

| Item | Before | After |
|---|---|---|
| P-D1 `.fast_preproject`, n=1200, x50 | 4.66 s | 0.05 s |
| P-D2 joint AR1, data passes | 3 | 2 |
| P-D2 joint AR1 + HAC, data passes | 4 | 2 |
| P-D3 HAC Parzen lag 4 (200 x 400) | 0.83 s | 0.70 s |
| P-E1 `rrr_gls` bootstrap, 100 resamples | 5.75 s | 1.74 s |
| P-B1 NIfTI group reads, S subjects, k chunks | 2·S·k | 2·S |
| P-C2 `(X'X)^-1` per model, 3 contrasts | 4 | 1 |
| `test_betas.R` (mixed-model betas) | 30.0 min | 34.6 min |

`test_betas.R` was already slow (each mixed-model test ~10 min); the mixed
path now carries the baseline as fixed effects (C1), and the after-run
overlapped other test runs.

### Not done / known limitations

- Satterthwaite df still forms an n x n residual covariance and its
  eigendecomposition per voxel (O(n^3)); P-D3 removed only the redundant work.
- Reference chunkwise backend still estimates AR per voxel chunk, so its
  results depend on `nchunks` (documented in stream D; pre-existing).
- Soft-projection df correction covers runwise iid fits only.
- `ols_t_vcov_cpp` (voxelwise covariate OLS) has no missing-data handling;
  `robust = "t"` is treated as none by the C++ meta kernels; the meta writer's
  `.stat_desc()` still drops `desc`.
- `create_fmri_model()` is still defined in three files (the `fmrilm.R`
  definition wins by collation order); consolidating it changes documented
  defaults and is left for a deliberate decision.
- `test_betas.R` runtime (mixed model per voxel) is unchanged in kind.

---

# Performance pass 2 (remaining hotspots)

Same ground rules as above, plus:

- Each change is result-preserving and gets an **equivalence test against an
  independent reference** (`lm`, `sandwich`, `fmrilss::mixed_solve`,
  `metafor`, a direct dense computation) at a stated tolerance.
- Each change also gets a **deterministic guard that fails on the old code**:
  a call count, a "never forms an n x n / T x V matrix" assertion, or a read
  count. No wall-clock assertions in tests.
- Timings are measured back to back (old vs new in the same R process, one
  thread) and recorded in `bench/review-2026-10-perf2.R` and the results
  table below; the bench script is build-ignored.
- Streams run targeted test files only; the full suite runs once after merge.

## Stream P1: fitting core

- [x] P1-1. Satterthwaite df: no per-voxel n x n residual covariance or
      `eigen()`; use the banded/Toeplitz structure (traces from lag sums).
      Reference: the current dense formula on small n.
- [x] P1-2. Voxelwise AR statistics: vectorise `beta_stats_matrix_voxelwise` /
      `fit_lm_contrasts_voxelwise` (`R/fmri_lm_internal.R`) over voxels;
      `future_lapply` in `R/fmri_lm_runwise.R` must not ship the full run
      matrix to every worker.
- [x] P1-3. Reference chunkwise backend: estimate AR once per run over all
      voxels and pass it to every chunk (results independent of `nchunks`).
- [x] P1-4. Robust IRLS: avoid the T x V copies per iteration
      (`crossprod(X * w, Y)` instead of forming `Y * sqrt(w)`).
- [x] P1-5. AR integration: `whiten_glm_context` stops computing unused OLS
      residuals; legacy whitening honours `fmrireg.ar.parallel_whiten` instead
      of `parallel = TRUE` per voxel.

## Stream P2: trial betas and outputs

- [x] P2-1. Mixed-model betas (`estimate_betas(method = "mixed")`): one
      spectral decomposition per design, then a 1-D REML search per voxel,
      instead of a full `mixed_solve` per voxel. Reference:
      `fmrilss::mixed_solve` per voxel. Target: `test_betas.R` from ~30 min to
      well under a minute.
- [x] P2-2. `coef_images()`: build the design names / mask space once, not per
      coefficient.
- [x] P2-3. `write_results(format = c("h5", "nifti"))`: compute each volume
      and sidecar once and hand it to every writer.
- [x] P2-4. Plot helpers: no silent fallback to loading the whole dataset to
      plot one voxel.

## Stream P3: group analysis and engines

- [x] P3-1. Meta kernels: compute the fixed-effects fit once per voxel and
      reuse it for DL / PM / Q (currently up to three times); vectorise the
      per-row DL loop.
- [x] P3-2. GDS meta path: one fit per call (no repeated `compute()` /
      manual refit / contrast refit); chunked.
- [x] P3-3. GDS t-combination: one vectorised call instead of a per-feature
      `vapply`.
- [x] P3-4. `rrr_gls` on latent data: stay in latent space instead of forming
      the dense T x V `Z %*% t(L)`.
- [x] P3-5. Soft projection: closed-form GCV (O(k) per lambda) and no unused
      right singular vectors.

## Performance pass 2 results

Test suite after pass 2 (`NOT_CRAN=true`, all files including
`test_betas.R`): 276 files, 10,047 expectations, 0 failures, 0 errors,
2 skips (same skips as the baseline); whole suite 5.5 min (`test_betas.R`
alone was 30 min). `R CMD check`: unchanged (same environmental 1 WARNING /
2 NOTEs; examples, Rd, code and compiled-code checks OK).

Timings: `bench/perf2-P{1,2,3}.R`, one thread, old and new back to back in
one process (old code copied in or sourced from 206600c); shared machine,
so treat small ratios as noise. "Max diff" is old vs new output.

| Item | Case | Old (s) | New (s) | Speed-up | Max diff |
|---|---|---|---|---|---|
| P1-1 Satterthwaite df | n=400, p=12, V=60, Parzen lag 8 | 5.42 | 0.006 | ~900x | ≤1e-10 |
| P1-1 | truncated window lag 8 (12/60 voxels need PSD clipping) | 5.45 | 0.44 | 12x | ≤1e-10 |
| P1-1 | truncated window lag 30 (all clipped) | 5.59 | 2.05 | 2.7x | ≤1e-10 |
| P1-2 voxelwise AR stats | 3 t + 2 F contrasts, V=5000 | 1.51 | 0.039 | 39x | ≤1e-13 |
| P1-3 reference backend | AR1, T=300, V=2000, 10 chunks | 0.40 | 0.08 | 5x | (now chunk-invariant) |
| P1-4 robust IRLS | Huber, T=400, V=5000 | 1.05-1.31 | 0.77-0.79 | 1.3-1.7x | ≤1e-13 |
| P1-5 AR integration | `whiten_glm_context`, fixed phi, V=20000 | 0.24 | 0.093 | 2.6x | ≤1e-13 |
| P2-1 mixed betas | n=300, 40 trials, 200 voxels | 11.80 | 0.122 | 97x | 7e-15 rel |
| P2-1 | n=900, 90 trials, 30 voxels | 45.87 | 0.965 | 48x | 2e-14 rel |
| P2-1 | `test_betas.R` | 1705 | 25 | 68x | 9/9 pass |
| P2-2 `coef_images` | 4 coefficients, 6,436 voxels | 0.019 | 0.006 | 3x | identical |
| P2-3 `write_results` h5+nifti | by_stat | 4.41 | 4.41 | 1x (I/O bound) | identical files |
| P3-1 meta kernels | DL, regression, P=20k | 0.180 | 0.120 | 1.5x | 0 |
| P3-2 GDS meta | PM ~group + contrast, P=3000 | 0.646 | 0.329 | 2x | 0 |
| P3-3 GDS t-combination | Lancaster, P=3000 | 0.571 | 0.084 | 6.8x | 0 |
| P3-4 `rrr_gls` latent | T=400, k=20, V=60k | 3.954 | 0.114 | 35x | 9e-14 |
| P3-4 | bootstrap(20) | 5.745 | 0.369 | 16x | 2e-13 |
| P3-5 soft projection GCV | n=300, p=400, V=20k | 42.8 | 1.28 | 33x | lambda 1e-9 |

Deterministic guards added (all failed on the old code): no dense n x n
covariance or `eigen()` on the Satterthwaite fast path; contrast functions
called once per contrast, not per voxel; `future_lapply` never ships the
run matrix; AR and IRLS estimated once per fit on the reference backend;
one `solve_glm_core` per IRLS fit and <= 3 T x V allocations per iteration
(Rprofmem); no pseudo-inverse for fixed phi; mixed-model design decomposed
once and `mixed_solve` never called per voxel; `coef_images` resolves names
and mask space once; `write_results` builds each volume/sidecar once; plot
readers never load the full dataset; one FE solve per voxel in the meta
kernel (test-only counter `meta_fit_diag_cpp`); one GDS fit per call,
`ceiling(P / chunk_size)` chunks; one t-combination call; latent `rrr_gls`
never forms T x V; GCV touches Y once and calls `svd(nv = 0)`.

Behaviour changes:
- Reference backend (`backend = "reference"`) results no longer depend on
  `voxel_chunks` (AR plan and robust weights are estimated once per fit). A
  robust model reads the full dataset once to form its weights.
- Mixed betas: a design-level failure warns once (not once per voxel);
  `progress` has no effect on the mixed path.

Not done:
- `rrr_gls` on latent data with AR noise still forms the dense response (AR
  is estimated from voxel residuals, a non-linear step).
- `write_results` time is dominated by fmristore's `gc()` on every H5 close
  (outside this package).
- Satterthwaite with a truncated (non-tapered) HAC window still needs an
  exact PSD projection, per run block, for voxels whose Toeplitz symbol is
  not provably positive; Parzen/Bartlett/Tukey stay on the fast path.

---

# Pass 3: remaining issues

Everything still open from the review, pass 1 "Not done" and pass 2 "Not
done". Same ground rules: a regression test that fails on the old code for
every behaviour change; reference-based tests for numerics; guards for
performance; no loosened tests; streams run targeted tests only.

Not addressed here (documented as inherent): Satterthwaite with an
untapered HAC window must keep the exact PSD projection for voxels whose
Toeplitz symbol is not positive (dropping it changes results);
`write_results` spends ~half its time in fmristore's `gc()` on H5 close
(upstream package).

## Stream R1: fitting and engines

- [x] R1-1. Soft-projection residual df beyond runwise iid: AR, robust and
      chunkwise paths use tr((I-H)P^2)-consistent df, or refuse the
      combination clearly. Test vs Monte Carlo E[RSS] / null calibration.
- [x] R1-2. `rrr_gls` on latent data with run/global-pooled AR: estimate phi
      from latent scores using the loadings Gram (as `fmri_latent_lm` does)
      and stay in latent space; voxelwise AR keeps the dense path or errors.
      Test vs dense fit; guard no T x V.
- [x] R1-3. `create_fmri_model()` defined in three files: keep the definition
      that runs today (collation winner), delete the others, make the docs
      match the real defaults. Test pins the defaults.
- [x] R1-4. Dead/inconsistent AR helpers: `fast_rlm_run` (MAD-scale, one
      caller), `.compute_ar_effective_df_compat` (double counts AR df),
      `exact_first %||% TRUE` vs FALSE elsewhere, repeated `ar_order`
      switch blocks -> one helper.
- [x] R1-5. Landmark extension (`fmri_lm_engine_lowrank`): stored
      `resvar`/`sigma` understate voxel residual variance; compute the real
      residual variance or store it under an honest name. Test vs dense.

## Stream R2: group analysis

- [x] R2-1. `ols_t_vcov_cpp` (voxelwise covariate OLS): drop non-finite
      subjects per voxel like `ols_t_cpp`. Test vs `lm` on complete cases.
- [x] R2-2. `robust = "t"`: accepted by `fmri_meta` and silently treated as
      none by the C++ kernels. Implement or reject everywhere with a clear
      error. `combine` with beta+SE input passes `method = "combine:..."`
      into the fitter and fails.
- [x] R2-3. `fmri_meta`: `n_threads` default in the generic not forwarded;
      `print.fmri_meta` labels I^2 (a fraction) as %; `parse_contrast_formula`
      silently returns all-zero weights for unmatched names (error instead).
- [x] R2-4. Meta writer `.stat_desc()` drops `desc` (same fix as E6 for
      fmri_lm outputs).
- [x] R2-5. GDS reducers declare `provides` fields some branches never
      return; registration errors are swallowed; registration is re-run on
      every call. Return what is declared; register once; surface failures.

## Stream R3: package hygiene and API

- [x] R3-1. DESCRIPTION: fmrigds -> Suggests (all uses already guarded);
      plotly/bslib -> Suggests with `design_plot()` guarded; drop unused
      `aes_string` import; `geom_line(size=)` -> `linewidth`; `design_plot`
      docs match its defaults.
- [x] R3-2. Exports: `.resample_param`, `.fmri_lm_formula_legacy` and the
      placeholder `estimate()` stub; `tidy` should be `generics::tidy` so
      `broom::tidy()` dispatches to fmrireg methods; one `%||%`.
- [x] R3-3. `export_jobs`: job-id uniqueness checked after sanitising (two ids
      mapping to one result file).
- [x] R3-4. `contrasts_api` `drop_failed = TRUE` drops misspelled contrasts
      silently (warn); `fmri_betas` stores the user's `basemod` (often NULL)
      instead of the baseline used; misleading "Cholesky failed" error on
      non-finite betas.
- [x] R3-5. `bootstrap_glm_inference` (internal, unwired, test-only; block
      bootstrap ignores runs, wild bootstrap destroys autocorrelation, BCa
      NaN): remove with its tests, or fix and test; prefer removal (one way
      to do one thing).

## Pass 3 results

Every item has a test in `tests/testthat/test-pass3-*.R` that fails on the
pre-pass-3 code.

Full suite after pass 3 (`NOT_CRAN=true`, all files incl. `test_betas.R`):

| Library | Files | Expectations | Failures | Errors | Skips |
|---|---|---|---|---|---|
| with fmrigds | 291 | 11,150 | 0 | 0 | 2 |
| without fmrigds | 291 | 9,876 | 0 | 0 | 46 |

`R CMD check --no-manual --no-tests --ignore-vignettes`: examples, Rd,
code/documentation, R code and compiled-code checks OK; only the known
environmental WARNING/NOTEs (locale, unavailable Suggests, installed size).

Follow-ups found while merging (each with a failing-first test):
- R3's fmrigds guard sat at the top of `group_data()` and also blocked the
  legacy NIfTI routes; it now guards only the gds routes. Found by running
  the whole suite against a library *without* fmrigds
  (`test-pass3-fmrigds-optional.R` mocks it unavailable). One test that
  silently needed fmrigds now skips; the `group_analysis` and
  `multisubject_fanout` vignettes skip their gds chunks without it (both
  knit with and without fmrigds).
- `export_jobs` ids that differ only in case collide on macOS/Windows file
  systems; now rejected.
- `latent_sketch` ignored `landmarks` for iid noise and for by-cluster AR;
  now applied for iid, refused for by-cluster.
- `projection_spec()` documents which fits support soft projection.
- Six mocks in two AR routing tests had fixed formals and rejected R1's new
  `acvf_correction` argument; they now take `...` (no assertion changed).

User-visible behaviour changes in pass 3:
- Soft projection: exact df/covariance and projection-corrected AR for
  runwise iid (incl. volume weights) and run/global-pooled AR; robust,
  joint/chunkwise, reference backend, MA, and AR + volume weights now error.
- Integrated-solver AR `effective_df` is n - rank (was reduced a second time
  for autocorrelation after whitening); missing `exact_first` means FALSE
  everywhere; an unknown AR structure errors.
- Landmark fits report each voxel's own residual variance (larger sigma/SE).
- `create_fmri_model()` has one definition; its docs describe the real
  defaults.
- `robust = "t"` errors; `combine` on beta+SE-only data errors; formula
  contrasts with unknown names error; voxelwise-covariate t-tests accept
  missing data (per-voxel df); a non-default `desc` appears in meta
  filenames; GDS reducers return exactly what they declare and register once.
- fmrigds, plotly and bslib are Suggests; `design_plot()`/`group_data()`
  error clearly without them. `tidy` is `generics::tidy`. `estimate()`,
  `.resample_param` and `.fmri_lm_formula_legacy` are no longer exported.
  `bootstrap_glm_inference` and `compute_contrast` are removed (unwired).
  Contrasts dropped by `drop_failed = TRUE` now warn.

Still open (inherent or upstream): Satterthwaite with an untapered HAC
window keeps the exact PSD projection for voxels whose Toeplitz symbol is
not provably positive; `write_results` time is dominated by fmristore's
`gc()` on H5 close.

---

# Pass 4: the two remaining limitations

## Satterthwaite df with an untapered HAC window — resolved

Voxels whose PSD projection cannot be certified no longer take a dense
eigendecomposition per run block. `src/satterthwaite_psd.cpp`
(`satterthwaite_psd_traces_cpp`, called once per fit for all such voxels):
banded Cholesky certificate that nothing is clipped; otherwise exact
centrosymmetric split, banded-to-tridiagonal reduction (dsbtrd-style Givens),
`dstebz` bisection for the extreme eigenvalues (tolerance over the whole V,
as the dense rule) and for every eigenvalue below it, `dstein` eigenvectors,
and low-rank corrections of tr(RV) and tr(RVRV). Never forms an n x n or
m x m matrix. Dense helper kept only as a per-voxel fallback on a LAPACK
failure (never triggered in tests or bench; tested by forcing it).

- Equivalence (`test-satterthwaite-psd-equivalence.R`): df within 1e-10
  relative of the dense reference across tapers x lags 0-40 x one/two runs/
  censoring, 0 / few / >= 20 clipped eigenvalues, tiny scales, and
  eigenvalues placed ~200 eps ||V|| either side of the tolerance; max
  observed difference 3.3e-15.
- Guard (`test-satterthwaite-psd-guard.R`): no `eigen()` calls (old: 12).
- Timing (`bench/satterthwaite-psd.R`, one thread, n = 400, p = 12): per
  projected voxel 110 ms -> 0.8-6.3 ms (one run), 34 ms -> 0.4-4.8 ms (two
  runs); 7-138x.

## H5 output closing — fmrireg side done, fmristore side pending

`.h5_close()` closes a result file directly when it is the only open HDF5
object, else falls back to `close_all()` (whose full `gc()` was ~half of
`write_results()` time). Tests: `test-h5-close.R`. With the pinned fmristore
the fallback is taken (no change). The upstream fix, which makes fmristore
writers close every handle they create and fixes nested-group creation, is
in `data-raw/fmristore-close-h5-handles.patch` (with its own failing-first
test). With it installed, `write_results(format = "h5")` on a 40x40x20 two-run
fit drops from 2.6-2.9 s to 1.0-1.3 s with identical voxel values
(`bench/h5-close.R`).

To finish: open the fmristore PR with that patch, then bump the fmristore pin
in DESCRIPTION `Remotes:` to the merged commit.

## Pass 4 results

Full suite (`NOT_CRAN=true`): with fmrigds 294 files, 11,298 expectations,
0 failures, 0 errors, 2 skips; without fmrigds 294 files, 10,024
expectations, 0 failures, 0 errors, 46 skips. `R CMD check`: examples, Rd,
code and compiled-code checks OK; only the known environmental
WARNING/NOTEs.
