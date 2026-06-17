# Regression tests for oneway / cell-oriented F-contrasts.
#
# oneway_contrast() emits weights in coefficient/cell orientation
# (k cells x (k-1) hypotheses) with a colind of length k. Three bugs used to
# make these unusable end to end:
#   1. fit_Fcontrasts() sized `cmat` as nrow(conmat) x ncoef and set the
#      numerator rank r = nrow(conmat), so the slow path crashed ("number of
#      items to replace ...") and would have used the wrong df even if it had
#      not.
#   2. create_full_contrast_matrix() rejected the k x (k-1) matrix outright, so
#      the fast path silently dropped the contrast.
#   3. .fast_F_contrast() computed the per-voxel quadratic form with colSums
#      (length r) instead of rowSums (length V), so even a reachable fast F
#      contrast produced wrong-length / wrong-value results.
# Additionally an F-contrast (which also inherits "contrast") was routed through
# the t-contrast engine, emitting a spurious "Failed to compute t-contrast" warning.

test_that(".orient_fcontrast normalises to hypotheses-in-rows", {
  W <- rbind(c(1, 0), c(0, 1), c(-1, -1)) # 3 cells x 2 hyps (cell-oriented)
  expect_equal(.orient_fcontrast(W, colind = 1:3), t(W))      # transposed -> 2 x 3
  expect_equal(.orient_fcontrast(t(W), colind = 1:3), t(W))   # already 2 x 3 -> identity
  expect_error(.orient_fcontrast(W, colind = 1:4), "do not match length")
})

test_that("fit_Fcontrasts matches anova reference and is orientation-invariant", {
  set.seed(42)
  n <- 60
  g <- factor(rep(c("A", "B", "C"), length.out = n))
  y <- rnorm(n)
  fit <- lm(y ~ 0 + g) # cell-means coding, 3 coefficients
  W <- rbind(c(1, 0), c(0, 1), c(-1, -1)) # oneway cell-oriented layout

  res_cell <- fit_Fcontrasts(fit, W, colind = 1:3)
  ref_F <- anova(lm(y ~ g))[1, "F value"] # factor F, numerator df = 2

  expect_equal(res_cell$stat, ref_F, tolerance = 1e-6)
  # Wrong r (= 3) would rescale the statistic; matching anova confirms r = 2.
  expect_equal(res_cell$prob, 1 - pf(res_cell$stat, 2, fit$df.residual), tolerance = 1e-9)

  res_col <- fit_Fcontrasts(fit, t(W), colind = 1:3) # column-oriented
  expect_equal(res_col$stat, res_cell$stat, tolerance = 1e-10)
})

test_that(".fast_F_contrast quadratic form is per-voxel (rowSums), r != V", {
  set.seed(11)
  p <- 3
  V <- 5
  r <- 2
  X <- matrix(rnorm(40 * p), 40, p)
  XtXinv <- solve(crossprod(X))
  B <- matrix(rnorm(p * V), p, V)
  sigma2 <- runif(V, 0.5, 2)
  L <- rbind(c(1, 0, -1), c(0, 1, -1)) # r x p

  out <- .fast_F_contrast(B, sigma2, XtXinv, L, df = 37)
  expect_length(out$stat, V) # colSums would give length r = 2

  M <- L %*% XtXinv %*% t(L)
  qf_ref <- vapply(seq_len(V), function(v) {
    u <- L %*% B[, v]
    drop(t(u) %*% solve(M) %*% u)
  }, numeric(1))
  F_ref <- (qf_ref / r) / sigma2
  expect_equal(out$stat, F_ref, tolerance = 1e-8)
})

make_oneway_dataset <- function(seed = 7, aliased = FALSE) {
  set.seed(seed)
  ev <- data.frame(
    onset = sort(runif(30, 0, 290)), run = 1,
    cond = factor(rep(c("A", "B", "C"), 10))
  )
  if (aliased) {
    ev$m <- rnorm(30)
    ev$m2 <- ev$m # exact duplicate -> rank-deficient nuisance
  }
  fmridataset::matrix_dataset(
    matrix(rnorm(150 * 4), 150, 4),
    TR = 2, run_length = 150, event_table = ev
  )
}

test_that("oneway F-contrast works end to end on both paths (full rank)", {
  dset <- make_oneway_dataset()
  con <- contrast_set(oneway_contrast(~cond, name = "cond_F"))
  form <- onset ~ hrf(cond, contrasts = con)

  warns <- character(0)
  fit_path <- function(fast) {
    withCallingHandlers(
      fmri_lm(form,
        block = ~run, dataset = dset, strategy = "runwise",
        use_fast_path = fast, progress = FALSE
      ),
      warning = function(w) {
        warns[[length(warns) + 1]] <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )
  }

  f_slow <- fit_path(FALSE)
  f_fast <- fit_path(TRUE)

  v_slow <- as.numeric(unlist(stats(f_slow, type = "F")[["cond_F"]]))
  v_fast <- as.numeric(unlist(stats(f_fast, type = "F")[["cond_F"]]))

  expect_true(length(v_slow) > 0 && all(is.finite(v_slow)))
  expect_true(length(v_fast) > 0 && all(is.finite(v_fast)))
  expect_equal(v_slow, v_fast, tolerance = 1e-6)
  # F-contrast must not be misrouted through the t-contrast engine.
  expect_false(any(grepl("t-contrast 'cond_F'", warns)))
  expect_false(any(grepl("Dimension mismatch", warns)))
})

test_that("estimable oneway F-contrast is not NA-poisoned under rank deficiency", {
  dset <- make_oneway_dataset(aliased = TRUE)
  con <- contrast_set(oneway_contrast(~cond, name = "cond_F"))
  form <- onset ~ hrf(cond, contrasts = con) + hrf(m) + hrf(m2)

  f_slow <- suppressWarnings(fmri_lm(form,
    block = ~run, dataset = dset,
    strategy = "runwise", use_fast_path = FALSE, progress = FALSE
  ))
  f_fast <- suppressWarnings(fmri_lm(form,
    block = ~run, dataset = dset,
    strategy = "runwise", use_fast_path = TRUE, progress = FALSE
  ))

  v_slow <- as.numeric(unlist(stats(f_slow, type = "F")[["cond_F"]]))
  v_fast <- as.numeric(unlist(stats(f_fast, type = "F")[["cond_F"]]))

  expect_true(length(v_slow) > 0 && all(is.finite(v_slow)))
  expect_true(length(v_fast) > 0 && all(is.finite(v_fast)))
  expect_equal(v_slow, v_fast, tolerance = 1e-6)
})

test_that("post-hoc fit_contrasts.fmri_lm computes a cell-oriented oneway F", {
  dset <- make_oneway_dataset()
  fit <- suppressWarnings(fmri_lm(onset ~ hrf(cond),
    block = ~run, dataset = dset,
    strategy = "runwise", use_fast_path = TRUE, progress = FALSE
  ))
  W <- rbind(c(1, 0), c(0, 1), c(-1, -1)) # cell-oriented, 3 x 2
  attr(W, "colind") <- fit$result$event_indices[1:3]

  ph <- fit_contrasts(fit, list(cond_F = W))
  expect_false(is.null(ph$cond_F)) # was silently dropped before the fix
  expect_true(all(is.finite(ph$cond_F$stat)))

  con <- contrast_set(oneway_contrast(~cond, name = "cond_F"))
  fit2 <- suppressWarnings(fmri_lm(onset ~ hrf(cond, contrasts = con),
    block = ~run, dataset = dset, strategy = "runwise",
    use_fast_path = TRUE, progress = FALSE
  ))
  v_ref <- as.numeric(unlist(stats(fit2, type = "F")[["cond_F"]]))
  expect_equal(as.numeric(ph$cond_F$stat), v_ref, tolerance = 1e-6)
})

test_that("voxelwise-AR path computes a cell-oriented oneway F", {
  dset <- make_oneway_dataset()
  con <- contrast_set(oneway_contrast(~cond, name = "cond_F"))
  form <- onset ~ hrf(cond, contrasts = con)
  f <- suppressWarnings(fmri_lm(form,
    block = ~run, dataset = dset, strategy = "runwise",
    use_fast_path = FALSE, cor_struct = "ar1", ar_voxelwise = TRUE, progress = FALSE
  ))
  v <- as.numeric(unlist(stats(f, type = "F")[["cond_F"]]))
  expect_true(length(v) > 0 && all(is.finite(v)))
})

test_that(".fast_F_contrast square case (r == V) uses rowSums, not colSums", {
  # When r == V the colSums/rowSums bug produced equal-length but WRONG values,
  # so finiteness checks alone would not catch it -- assert exact values.
  set.seed(99)
  p <- 4
  V <- 2
  r <- 2
  X <- matrix(rnorm(50 * p), 50, p)
  XtXinv <- solve(crossprod(X))
  B <- matrix(rnorm(p * V), p, V)
  sigma2 <- c(1.3, 0.7)
  L <- rbind(c(1, -1, 0, 0), c(0, 0, 1, -1)) # r x p

  out <- .fast_F_contrast(B, sigma2, XtXinv, L, df = 46)
  M <- L %*% XtXinv %*% t(L)
  qf_ref <- vapply(seq_len(V), function(v) {
    u <- L %*% B[, v]
    drop(t(u) %*% solve(M) %*% u)
  }, numeric(1))
  expect_equal(out$stat, (qf_ref / r) / sigma2, tolerance = 1e-8)
})
