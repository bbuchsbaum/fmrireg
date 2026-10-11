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

- [ ] A1. Empirical null: scale estimated from the central 60% of |z| is
      ~0.57 under the null; correct for truncation (or use all finite z).
      Test: 1e5 pure-noise z, default args, rejections at alpha=.05 ≈ 0;
      sigma0 within 0.95–1.05.
- [ ] A2. Group-ID mismatch between R (`factor` order) and C++ (first-
      appearance order). Make one numbering authoritative; vectorise the
      weight expansion. Test: result invariant to relabelling/permuting group
      labels; neighbor smoothing uses the right groups.

## Stream B: group analysis and meta-analysis

- [ ] B1. `meta_stouffer` (run-pooled F-contrasts) returns a weighted mean,
      not a Stouffer z. Use Σwz/√Σw², upper-tail `qnorm`/`pnorm`. Test:
      k runs at identical p give z = √k·qnorm(1-p); null calibration.
- [ ] B2. `weights = "equal"/"custom"` treated as absolute variances. Treat as
      relative weights with an estimated residual scale. Custom weights must
      be sliced per chunk. Test: equal-weights FE equals one-sample t
      (`t.test`) for any data scale; custom-weights equals weighted LS
      (`lm(weights=)`); chunk_size does not change results.
- [ ] B3. t→z conversion overflows to Inf for p < ~1e-16 (7 copies). One
      helper using log-scale tails. Test: t = 40, df = 20 gives finite z
      matching `qnorm(pt(..., log.p=TRUE), log.p=TRUE)`.
- [ ] B4. Stouffer/Lancaster t-combination counts missing subjects. Test:
      result with NA subjects equals result with those subjects removed.
- [ ] B5. Meta estimators verified against `metafor::rma` (FE, DL, PM) for
      intercept-only and meta-regression; label "reml" honestly (either
      implement or map with a message). Test: estimates/SE/tau2/Q within 1e-6.
- [ ] B6. Contrast SE ignores coefficient covariance by default. Test against
      `metafor` contrast SE (`anova(rma, X=)`).
- [ ] B7. Classic t-test engine: NA subject wipes the voxel. Drop non-finite
      subjects per voxel (as the meta engine does). Test vs `t.test` on the
      complete cases.
- [ ] B8. Paired test with `voxelwise_cov` errors on rownames. Test that it
      runs and matches `lm` per voxel.
- [ ] B9. HDF5 reader leaks handles / silences failures. Test handle closes
      and that a bad file is an error, not NA.
- [ ] P-B1. NIfTI group reader re-reads every subject per chunk. Read once.
      Test: number of reads == number of subjects (mock the reader); result
      identical across chunk sizes.
- [ ] P-B2. `schedule(dynamic)` in meta kernels; no result change (existing
      tests + B5 references).

## Stream C: betas, contrasts, tidy, design utilities

- [ ] C1. `estimate_betas(method = "ols"/"mixed")` residualises Y but not X.
      Test: noiseless simulation recovers true betas exactly (vs `lm`).
- [ ] C2. `tidy(type = "contrasts")` returns only the first contrast. Test:
      two t-contrasts + one F-contrast, all present with correct values.
- [ ] C3. `generate_main_effect_contrast` / `generate_interaction_contrast` row
      order vs `expand.grid`. Test: contrast applied to cell means from a
      known factorial design yields the textbook main effect / interaction.
- [ ] C4. Post-hoc `fit_contrasts.fmri_lm`: silent positional fallback and
      dropped specs. Map by full column names, error when unmatched, support
      every spec via `contrast_weights()`. Test: post-hoc contrasts equal the
      ones computed during the fit, for pair, oneway and unit specs, with
      several terms.
- [ ] C5. Post-hoc contrasts on AR/robust chunkwise fits mix unwhitened
      `(X'X)^-1` with whitened sigma. Store the matching covariance, or refuse.
      Test: post-hoc equals in-fit contrast for an ar1 fit.
- [ ] C6. `hrf_smoothing_kernel` normalisation (asymmetric, NaN). Test:
      symmetric, unit diagonal, finite with buffer_scans = 0.
- [ ] C7. Absolute 1.5e-8 / sqrt(eps) thresholds zero out small-scale data.
      Test: t/F invariant to rescaling Y by 1e-4.
- [ ] C8. F p-value `1 - pf()` loses precision. Use `lower.tail = FALSE`.
- [ ] P-C1. LSS fixed effects and mixed path: one multi-response solve instead
      of per-voxel `lm.fit`/`lsfit`; drop identity `map_voxels`. Equivalence
      test vs previous per-voxel results.
- [ ] P-C2. Reference F-contrast path: vectorise; compute `cov.unscaled` once.
      Equivalence test vs `.fast_F_contrast` and `anova(lm)`.

## Stream D: GLM fitting core

- [ ] D1. `reestimate_phi = TRUE` estimates AR from whitened residuals.
      Use unwhitened residuals. Test: AR(1) phi = 0.5 simulation, robust +
      reestimate_phi recovers phi ≈ 0.5 and matches non-reestimated fit.
- [ ] D2. Robust df counts zero-weight volumes; chunkwise vs runwise use
      different sigma; returned weights differ from the weights used for the
      final betas. Test: runwise == chunkwise for robust fits on 1 and 2
      runs; betas reproducible from returned weights by WLS.
- [ ] D3. Runwise volume weights break on multi-run data. Test: two runs,
      full-length weight vector, matches WLS `lm(weights=)` on stacked data
      with run intercepts; zero weights reduce df.
- [ ] D4. Effective-df: one implementation, the one CLAUDE.md documents; delete
      unreachable copies. Test against hand computation.
- [ ] D5. Rank-deficient designs: `n <- dfres + p` wrong. Test with an aliased
      column: df matches `lm`.
- [ ] P-D1. `.fast_preproject` builds an n×n identity on every call. Compute
      the pseudo-inverse lazily. Test: results identical; `Pinv` never built on
      the OLS/voxelwise-AR path (trace calls).
- [ ] P-D2. Chunking loads the full dataset (several times). Lazy chunk
      descriptors, run rows from the temporal schema, no `Y_orig`. Test: count
      data-access calls on a mock dataset; results identical across
      `nchunks`.
- [ ] P-D3. HAC/Satterthwaite: skip covariance when unused; hoist bread and
      `XtXinv` out of the voxel loop. Equivalence test vs `sandwich::NeweyWest`
      / `vcovHAC` for a few voxels.
- [ ] P-D4. `solve_glm_core`: single `qr.qty` pass for betas and exact RSS.
      Equivalence vs `lm`.

## Stream E: engines, latent, I/O, CRAN hygiene

- [ ] E1. `rrr_gls` SEs ~√(V/k) too small. Use per-voxel residual variance.
      Test: pure-noise calibration (|t|>1.96 rate ≈ 5%); full-rank result
      equals OLS SEs.
- [ ] E2. `fmri_latent_lm` silently ignores `robust` / `autocor`. Read the
      real config fields; error for unsupported combinations; remove the
      undefined-function globals. Test: AR(1) latent fit differs from OLS and
      matches dense AR fit on the reconstructed data.
- [ ] E3. RNG: no `set.seed` without restoring state (`rrr_gls`, simulators,
      `bids_export` staging dir). Test: `.Random.seed` unchanged after calls.
- [ ] E4. `write_results(format = "nifti")` must not require fmristore.
      Test with `requireNamespace` mocked to FALSE.
- [ ] E5. Benchmark data must not hard-fail on fmrihrf version mismatch.
      Test with the version mocked.
- [ ] E6. BIDS sidecar DegreesOfFreedom scalar (or range), contrast-label
      collisions detected, `desc` sanitised and kept. Tests on file names and
      JSON content.
- [ ] P-E1. `rrr_gls` bootstrap: hoist QR of the fixed design. Equivalence
      test with fixed seed.

## Final pass

- [ ] F1. Merge streams; full test suite green.
- [ ] F2. Remove dead code proven unreachable (list in commit), keeping tests
      that exercised live behaviour.
- [ ] F3. `R CMD build` + `R CMD check --no-manual` with no new warnings.
- [ ] F4. Update CLAUDE.md with what changed and the new reference tests.
- [ ] F5. Results table below: test counts before/after, timings of the
      performance items.

## Results

(filled in at the end)
