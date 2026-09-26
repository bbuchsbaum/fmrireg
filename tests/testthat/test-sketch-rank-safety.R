# Rank safety of the time-sketched GLM (engine = "latent_sketch").
#
# Keep-but-aliased: estimability is judged on the full (whitened) design with
# the same pivoted-QR rule as exact OLS; aliased coefficients are NA, the
# estimable ones are the sketch solve of the reduced design, and contrasts
# that load on an aliased column are NA with a named warning. Before this was
# enforced, chol() of a numerically singular sketched Gram matrix "succeeded":
# with a duplicated nuisance column CountSketch reported sigma2 ~ 221 (truth 1)
# and event standard errors of 38-166 (truth ~3.6), and IHS reported aliased
# variances of ~2e13.

.rank_case_data <- function(case, Tlen = 160L, nvox = 150L, seed = 1L) {
  TR <- 2
  on <- seq(6, Tlen * TR - 30, by = 16)
  cond <- rep(c("A", "B"), length.out = length(on))
  ev_full <- data.frame(onsets = on, condition = cond, run = 1L)
  if (identical(case, "empty")) {
    # Condition C is declared but every onset lies after the last scan, so
    # its convolved regressor is identically zero.
    ev_alias <- rbind(ev_full, data.frame(onsets = Tlen * TR + c(20, 60),
                                          condition = "C", run = 1L))
  } else {
    ev_alias <- ev_full
  }
  ev_full$condition <- factor(ev_full$condition)
  ev_alias$condition <- factor(ev_alias$condition)
  Y <- withr::with_seed(seed, matrix(stats::rnorm(Tlen * nvox), Tlen) + 100)
  sf <- fmridesign::sampling_frame(blocklens = Tlen, TR = TR)
  tt <- seq_len(Tlen)
  n12 <- cbind(n1 = sin(tt / 7), n2 = cos(tt / 11))
  nuis_alias <- if (identical(case, "duplicate")) cbind(n12, n3 = n12[, "n1"]) else n12
  list(
    alias = list(
      ds = matrix_frame(Y, TR = TR, run_length = Tlen, event_table = ev_alias),
      bm = suppressWarnings(fmridesign::baseline_model(
        basis = "poly", degree = 1, sframe = sf, nuisance_list = list(nuis_alias)))
    ),
    reduced = list(
      ds = matrix_frame(Y, TR = TR, run_length = Tlen, event_table = ev_full),
      bm = fmridesign::baseline_model(
        basis = "poly", degree = 1, sframe = sf, nuisance_list = list(n12))
    ),
    nvox = nvox, Tlen = Tlen
  )
}

.rank_fit <- function(d, method, noise, m = 40L, contrasts = TRUE, seed = 7L) {
  parcels <- rep(1:3, length.out = ncol(fmridataset::collect_assay(d$ds)))
  ctl <- switch(noise,
    iid = fmri_lm_control(),
    ar1 = fmri_lm_control(noise = noise_spec("ar1")),
    by_cluster = fmri_lm_control(noise = noise_spec("ar1", pooling = "parcel",
                                                    parcels = parcels))
  )
  lr <- lowrank_control(
    parcels = if (identical(noise, "by_cluster")) parcels else NULL,
    time_sketch = list(method = method, m = m)
  )
  cons <- list(pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"))
  if (contrasts) {
    cons <- c(cons, list(pair_contrast(~ condition == "A", ~ condition == "C",
                                       name = "A_vs_C")))
  }
  con <- do.call(contrast_set, cons)
  form <- onsets ~ hrf(condition, contrasts = con)
  warns <- character(0)
  fit <- withCallingHandlers(
    withr::with_seed(seed, fmri_lm(
      form, block = ~run, dataset = d$ds, baseline_model = d$bm,
      engine = "latent_sketch", control = ctl, lowrank = lr,
      drop_empty = FALSE
    )),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = warns)
}

.rank_payload <- function(fit) {
  b <- fit$result$betas$data[[1]]
  list(est = as.matrix(b$estimate[[1]]), se = as.matrix(b$se[[1]]))
}

.rank_contrast <- function(fit, name) {
  ct <- fit$result$contrasts
  hit <- grep(name, ct$name, fixed = TRUE)
  stopifnot(length(hit) == 1L)
  ct$data[[hit]]
}

.rank_covs <- function(fit) {
  fit$result$covariance_by_cluster %||% list(fit$result$cov.unscaled)
}

for (case in c("duplicate", "empty")) {
  test_that(sprintf("sketched fits alias the %s column and solve the rest exactly", case), {
    d <- .rank_case_data(case)
    ols <- fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = d$reduced$ds,
                   baseline_model = d$reduced$bm)
    ols_se <- as.matrix(.rank_payload(ols)$se)[, 1:2]
    for (noise in c("iid", "ar1", "by_cluster")) {
      for (method in c("gaussian", "countsketch", "srht", "ihs")) {
        lab <- paste(case, noise, method)
        fa <- .rank_fit(d$alias, method, noise, contrasts = identical(case, "empty"))
        fr <- .rank_fit(d$reduced, method, noise, contrasts = FALSE)
        fit <- fa$fit; ref <- fr$fit
        vn <- colnames(design_matrix(fit$model))
        al <- if (identical(case, "empty")) {
          which(vn == "condition_condition.C")
        } else {
          which(vn == "nuis_n3_block_1")
        }
        expect_length(al, 1L)
        keep <- setdiff(seq_along(vn), al)

        # The aliased column is marked exactly as exact OLS marks it.
        for (cv in .rank_covs(fit)) {
          expect_identical(as.integer(attr(cv, "aliased")), al, label = lab)
          expect_false(attr(cv, "is_full_rank"), label = lab)
          expect_true(all(is.finite(cv[keep, keep])), label = lab)
        }
        expect_true(any(grepl("rank deficient", fa$warnings)), label = lab)

        # Aliased coefficient and SE are NA; everything else is the sketch
        # solve of the reduced design with the same sketch.
        b <- .rank_payload(fit); r <- .rank_payload(ref)
        expect_true(all(is.na(b$est[, al])) && all(is.na(b$se[, al])), label = lab)
        expect_true(all(is.na(fit$betas_fixed[al, ])), label = lab)
        expect_equal(unname(b$est[, keep]), unname(r$est), tolerance = 1e-7, label = lab)
        expect_equal(unname(b$se[, keep]), unname(r$se), tolerance = 1e-7, label = lab)
        expect_equal(fit$sigma2, ref$sigma2, tolerance = 1e-7, label = lab)
        expect_equal(fit$result$rdf, ref$result$rdf, tolerance = 1e-7, label = lab)
        expect_lt(abs(mean(fit$sigma2) - 1), 0.15, label = lab)
        # Event SEs stay on the scale of the reduced exact OLS fit: about
        # sqrt(T / m) = 2 times it for sketch-and-solve, equal for IHS. The
        # pre-fix CountSketch SEs were 20-100 times OLS.
        se_ratio <- mean(b$se[, 1:2]) / mean(ols_se)
        expect_gt(se_ratio, 0.9, label = lab)
        expect_lt(se_ratio, 3 * sqrt(d$Tlen / 40), label = lab)

        # An estimable contrast is unaffected by the aliased column.
        ab <- .rank_contrast(fit, "A_vs_B"); ab_ref <- .rank_contrast(ref, "A_vs_B")
        expect_equal(ab$estimate, ab_ref$estimate, tolerance = 1e-7, label = lab)
        expect_equal(ab$se, ab_ref$se, tolerance = 1e-7, label = lab)
        expect_true(all(is.finite(ab$se)), label = lab)

        if (identical(case, "empty")) {
          # A contrast loading on the aliased condition is NA with a warning
          # naming the contrast and the column.
          ac <- .rank_contrast(fit, "A_vs_C")
          expect_true(all(is.na(ac$estimate)) && all(is.na(ac$se)), label = lab)
          # ... once, not once per cluster.
          expect_equal(sum(grepl("A_vs_C.*non-estimable.*condition_condition.C",
                                 fa$warnings)), 1L, label = lab)
          if (!identical(noise, "by_cluster")) {
            expect_warning(
              pc <- fit_contrasts(fit, list(C_only = structure(1, colind = al))),
              "non-estimable"
            )
            expect_true(all(is.na(pc$C_only$estimate)), label = lab)
          }
        }
      }
    }
  })
}

test_that("a sketch too small for the estimable design is an error naming m", {
  # Row-sampling sketch that keeps the first m scans: the third column is
  # nonzero only after scan m, so the sketched design loses rank although
  # the full design is well conditioned. No ridge fallback may hide this.
  Tlen <- 60L; m <- 20L
  X <- cbind(1, seq_len(Tlen) / Tlen, c(rep(0, 40L), rep(1, 20L)))
  Z <- matrix(stats::rnorm(Tlen * 3), Tlen)
  op <- list(method = "rows", m = m,
             apply = function(A) as.matrix(A)[seq_len(m), , drop = FALSE],
             gram = function() diag(m))
  expect_error(
    fmrireg:::.lowrank_time_solve(X, Z, list(method = "rows", m = m), op),
    "time_sketch\\$m = 20 is too small"
  )
  # The same design with a sketch that sees every column solves.
  op_ok <- op
  op_ok$apply <- function(A) as.matrix(A)[seq(1L, Tlen, by = 3L), , drop = FALSE]
  sol <- fmrireg:::.lowrank_time_solve(X, Z, list(method = "rows", m = m), op_ok)
  expect_true(all(is.finite(sol$M)))
  expect_length(attr(sol$cov_unscaled, "aliased"), 0L)
})

test_that("IHS with fewer sketch rows than estimable columns is an error naming m", {
  # Every IHS Hessian sketch is then singular. The kernel used to fall back
  # to pinv() and returned coefficients hundreds of OLS SEs off while
  # reporting the exact OLS covariance.
  dset <- sketch_matrix_dataset(run_length = 60L, nvox = 10L, seed = 3L)
  for (m in c(3L, 6L)) {
    expect_error(
      suppressWarnings(fmri_lm(
        onsets ~ hrf(condition), block = ~run, dataset = dset,
        engine = "latent_sketch",
        lowrank = lowrank_control(time_sketch = list(method = "ihs", m = m,
                                                     iters = 4L, tol = 0))
      )),
      sprintf("time_sketch\\$m = %d is too small", m)
    )
  }
  # Aliased columns do not count towards the m requirement, and the fit
  # reproduces OLS on the reduced design.
  Tlen <- 80L
  X <- cbind(1, sin(seq_len(Tlen) / 5), cos(seq_len(Tlen) / 9))
  X <- cbind(X, X[, 2])
  Z <- withr::with_seed(4, matrix(stats::rnorm(Tlen * 4), Tlen))
  expect_no_error(withr::with_seed(5, fmrireg:::.lowrank_time_solve(
    X, Z, list(method = "ihs", m = 3L, iters = 2L, tol = 0), NULL)))
  sol <- withr::with_seed(5, fmrireg:::.lowrank_time_solve(
    X, Z, list(method = "ihs", m = 20L, iters = 100L, tol = 1e-8), NULL))
  expect_identical(as.integer(attr(sol$cov_unscaled, "aliased")), 4L)
  expect_equal(sol$M[1:3, ], qr.solve(X[, 1:3], Z), tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_equal(sol$cov_unscaled[1:3, 1:3], solve(crossprod(X[, 1:3])),
               tolerance = 1e-10, ignore_attr = TRUE)
})
