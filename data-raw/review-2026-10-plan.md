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
