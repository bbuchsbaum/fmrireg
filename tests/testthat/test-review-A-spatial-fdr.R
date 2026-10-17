# Review 2026-10, stream A: spatial FDR correctness
#
# A1: the default empirical null must be consistent under N(0, 1) and robust
#     to a minority of true signals.
# A2: one group numbering (sorted labels / factor levels) is authoritative
#     for pi0, smoothing, weights and the returned `groups`.
# Plus: the weighted-BH step is checked against a direct R implementation.

# Reference weighted BH (Genovese, Roeder & Wasserman 2006): weights are
# normalised to sum to the number of valid tests, then BH is run on p / w.
wbh_reference <- function(p, w, alpha) {
  valid <- is.finite(p)
  n <- sum(valid)
  wn <- numeric(length(p))
  wn[valid] <- w[valid] * n / sum(w[valid])
  ps <- rep(NA_real_, length(p))
  ps[valid] <- p[valid] / wn[valid]
  srt <- sort(ps[valid])
  ok <- which(srt <= alpha * seq_len(n) / n)
  k <- if (length(ok)) max(ok) else 0L
  thr <- if (k > 0) alpha * k / n else 0
  reject <- valid & !is.na(ps) & ps <= thr & k > 0
  q <- rep(NA_real_, length(p))
  q[valid] <- stats::p.adjust(pmin(ps[valid], 1), method = "BH")
  list(reject = reject, q = q, w_norm = wn, k = k, threshold = thr)
}

# Direct R version of the groupwise Storey pi0 estimate, indexed by
# sorted group label.
pi0_reference <- function(p, group, tau) {
  f <- if (is.factor(group)) group else factor(group)
  tail <- tapply(is.finite(p) & p > tau, f, sum)
  m <- tapply(is.finite(p), f, sum)
  est <- as.numeric(tail) / ((1 - tau) * as.numeric(m))
  est[is.na(est)] <- 1
  pmin(pmax(est, 0), 1)
}

# ---------------------------------------------------------------------------
# A1. Empirical null calibration
# ---------------------------------------------------------------------------

test_that("A1: empirical null is consistent under a pure N(0,1) null", {
  set.seed(20261011)
  m <- 1e5
  z <- rnorm(m)
  grp <- rep(seq_len(100), each = m / 100)

  res <- spatial_fdr(z = z, group = grp)

  expect_gte(res$sigma0, 0.95)
  expect_lte(res$sigma0, 1.05)
  expect_lt(abs(res$mu0), 0.05)
  # Under the global null BH makes any rejection with probability <= alpha;
  # the old estimator (sigma0 ~ 0.57) rejected ~11,700 of these.
  expect_lte(sum(res$reject), 5)
})

test_that("A1: null calibration holds across replicates (binomial band)", {
  set.seed(7)
  nrep <- 40
  m <- 1e4
  grp <- rep(seq_len(20), each = m / 20)
  any_rej <- logical(nrep)
  sig <- numeric(nrep)
  for (r in seq_len(nrep)) {
    res <- spatial_fdr(z = rnorm(m), group = grp)
    any_rej[r] <- any(res$reject)
    sig[r] <- res$sigma0
  }
  # P(any rejection) <= alpha = .05 for BH under the global null; allow the
  # upper 99.9% binomial quantile.
  expect_lte(sum(any_rej), stats::qbinom(0.999, nrep, 0.05))
  expect_lt(abs(mean(sig) - 1), 0.02)
})

test_that("A1: empirical null stays near N(0,1) with 10% true signals and keeps power", {
  set.seed(42)
  m <- 1e5
  G <- 100
  grp <- rep(seq_len(G), each = m / G)
  is_sig <- grp <= 10                      # 10% of features, clustered
  z <- rnorm(m)
  z[is_sig] <- rnorm(sum(is_sig), mean = 4)

  res <- spatial_fdr(z = z, group = grp)
  oracle <- spatial_fdr(z = z, group = grp, empirical_null = FALSE)

  expect_gte(res$sigma0, 0.95)
  expect_lte(res$sigma0, 1.05)
  expect_lt(abs(res$mu0), 0.05)

  power <- mean(res$reject[is_sig])
  power_oracle <- mean(oracle$reject[is_sig])
  expect_gt(power, power_oracle - 0.03)
  fdp <- sum(res$reject & !is_sig) / max(1, sum(res$reject))
  expect_lt(fdp, 0.08)
})

test_that("A1: shifted/scaled null is recovered", {
  set.seed(111)
  z <- rnorm(5e4, mean = 0.5, sd = 1.2)
  res <- spatial_fdr(z = z, group = rep(1:50, each = 1000))
  expect_equal(res$mu0, 0.5, tolerance = 0.05 / 0.5)
  expect_equal(res$sigma0, 1.2, tolerance = 0.05)
})

test_that("A1: small-sample fallback (robust MAD on all finite z) is consistent", {
  set.seed(3)
  sig <- replicate(200, spatial_fdr(z = rnorm(200), group = rep(1:4, 50))$sigma0)
  expect_lt(abs(mean(sig) - 1), 0.03)
})

# ---------------------------------------------------------------------------
# A2. Group numbering
# ---------------------------------------------------------------------------

test_that("A2: pi0 and weights follow sorted group labels, not first appearance", {
  set.seed(1)
  n <- 200
  # Group 2 appears first; the signal is in group 1.
  group <- c(rep(2, n), rep(1, n))
  p <- c(runif(n), rbeta(n, 0.1, 5))

  res <- spatial_fdr(p = p, group = group, lambda = 0, empirical_null = FALSE)

  expect_identical(as.integer(res$groups), as.integer(factor(group)))
  expect_equal(res$pi0_raw, pi0_reference(p, group, 0.5), tolerance = 1e-12)
  expect_equal(res$m_g, c(n, n))
  # Signal block (group 1) must get the larger weight
  expect_gt(mean(res$weights[group == 1]), mean(res$weights[group == 2]))
  # Null block must not be up-weighted
  expect_lte(max(res$weights[group == 2]), 1)
  # Weights equal 1/pi0_smooth of each feature's own group, normalised
  a <- 1 / res$pi0_smooth[as.integer(factor(group))]
  expect_equal(res$weights, a * length(p) / sum(a), tolerance = 1e-12)
})

test_that("A2: results are invariant to relabelling groups and permuting features", {
  set.seed(99)
  G <- 12
  per <- 60
  group <- rep(seq_len(G), each = per)
  z <- rnorm(G * per)
  z[group %in% c(3, 4, 5)] <- z[group %in% c(3, 4, 5)] + 3
  neighbors <- lapply(seq_len(G), function(g) setdiff(c(g - 1, g + 1), c(0, G + 1)))

  base <- spatial_fdr(z = z, group = group, neighbors = neighbors, lambda = 1)

  # Relabel: group g -> perm[g]; neighbors relabelled consistently
  perm <- sample(G)
  group2 <- perm[group]
  nb2 <- vector("list", G)
  for (g in seq_len(G)) nb2[[perm[g]]] <- perm[neighbors[[g]]]
  rel <- spatial_fdr(z = z, group = group2, neighbors = nb2, lambda = 1)
  expect_identical(rel$reject, base$reject)
  expect_equal(rel$q, base$q, tolerance = 1e-12)
  expect_equal(rel$weights, base$weights, tolerance = 1e-12)
  expect_equal(rel$pi0_smooth[perm], base$pi0_smooth, tolerance = 1e-12)

  # Character labels map through sorted levels the same way
  lab <- sprintf("g%02d", group)
  chr <- spatial_fdr(z = z, group = lab, neighbors = neighbors, lambda = 1)
  expect_identical(chr$reject, base$reject)
  expect_equal(chr$q, base$q, tolerance = 1e-12)

  # Permute feature order, then un-permute the outputs
  o <- sample(length(z))
  pf <- spatial_fdr(z = z[o], group = group[o], neighbors = neighbors, lambda = 1)
  inv <- order(o)
  expect_identical(pf$reject[inv], base$reject)
  expect_equal(pf$q[inv], base$q, tolerance = 1e-12)
  expect_equal(pf$weights[inv], base$weights, tolerance = 1e-12)
  expect_equal(pf$pi0_smooth, base$pi0_smooth, tolerance = 1e-12)
  expect_identical(pf$groups[inv], base$groups)
})

test_that("A2: neighbor smoothing uses the groups the neighbors list names", {
  # Chain 1 - 2 - 3. Data order puts group 3 first, then 1, then 2, so a
  # first-appearance numbering would smooth the wrong pairs.
  p <- c(rep(0.9, 10),            # group 3: pi0_raw = 1
         rep(0.001, 10),          # group 1: pi0_raw = 0
         c(rep(0.9, 5), rep(0.001, 5)))  # group 2: pi0_raw = 1 (5/(0.5*10))
  group <- c(rep(3, 10), rep(1, 10), rep(2, 10))
  neighbors <- list(2L, c(1L, 3L), 2L)

  res <- spatial_fdr(p = p, group = group, neighbors = neighbors, lambda = 1,
                     min_pi0 = 0.05, empirical_null = FALSE)
  raw <- pi0_reference(p, group, 0.5)
  expect_equal(raw, c(0, 1, 1))
  expect_equal(res$pi0_raw, raw)
  expected <- c((raw[1] + raw[2]) / 2,
                (raw[2] + mean(raw[c(1, 3)])) / 2,
                (raw[3] + raw[2]) / 2)
  expect_equal(res$pi0_smooth, pmax(expected, 0.05), tolerance = 1e-12)
  # Group 3 (all p = 0.9) is adjacent only to group 2, so its smoothed pi0
  # stays at 1; under the first-appearance numbering it would be 0.5.
  expect_equal(res$pi0_smooth[3], 1)
})

test_that("A2: factor levels (including unused ones) define the numbering", {
  set.seed(5)
  f <- factor(rep(c("a", "c"), each = 50), levels = c("a", "b", "c"))
  p <- c(rbeta(50, 0.2, 4), runif(50))
  neighbors <- list(2L, c(1L, 3L), 2L)
  res <- spatial_fdr(p = p, group = f, neighbors = neighbors, lambda = 1,
                     empirical_null = FALSE)
  expect_equal(res$G, 3L)
  expect_equal(res$m_g, c(50L, 0L, 50L))
  expect_identical(as.integer(res$groups), as.integer(f))
  expect_equal(res$pi0_raw, pi0_reference(p, f, 0.5), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Weighted BH against a direct R implementation, with NA handling
# ---------------------------------------------------------------------------

test_that("weighted BH matches a direct R implementation (with NAs)", {
  set.seed(2026)
  m <- 2000
  p <- c(rbeta(300, 0.2, 6), runif(m - 300))
  w <- rexp(m) + 0.1
  p[sample(m, 40)] <- NA

  for (alpha in c(0.01, 0.05, 0.2)) {
    got <- fmrigds:::sabha_weighted_bh_cpp(p, w, alpha = alpha)
    ref <- wbh_reference(p, w, alpha)
    expect_identical(as.logical(got$reject), ref$reject)
    expect_equal(got$k, ref$k)
    expect_equal(got$threshold, ref$threshold, tolerance = 1e-12)
    expect_equal(got$w_norm, ref$w_norm, tolerance = 1e-12)
    expect_equal(sum(got$w_norm[is.finite(p)]), sum(is.finite(p)))

    q <- fmrigds:::sabha_bh_qvalues_scaled_cpp(p / got$w_norm)
    expect_equal(q, ref$q, tolerance = 1e-12)
    expect_identical(ref$reject, !is.na(q) & q <= alpha)
  }
})

test_that("spatial_fdr end-to-end equals weighted BH with 1/pi0 group weights; NA handling", {
  set.seed(77)
  G <- 8
  group <- sample(rep(seq_len(G), each = 100))
  p <- runif(length(group))
  sig <- group %in% c(2, 7)
  p[sig] <- rbeta(sum(sig), 0.15, 5)
  p[c(5, 17, 230)] <- NA
  group[c(40, 41)] <- NA

  res <- spatial_fdr(p = p, group = group, lambda = 0, min_pi0 = 0.05,
                     empirical_null = FALSE, alpha = 0.05)

  pi0 <- pmin(pmax(pi0_reference(p[!is.na(group)], group[!is.na(group)], 0.5),
                   0.05), 1)
  expect_equal(res$pi0_smooth, pi0, tolerance = 1e-12)
  w <- ifelse(is.na(group), 1, 1 / pi0[group])
  ref <- wbh_reference(p, w, 0.05)
  expect_identical(res$reject, ref$reject)
  expect_equal(res$q, ref$q, tolerance = 1e-12)
  expect_equal(res$weights[is.finite(p)], ref$w_norm[is.finite(p)], tolerance = 1e-12)

  # NA p: never rejected, NA q; NA group: excluded from pi0, weight 1 pre-normalisation
  expect_false(any(res$reject[is.na(p)]))
  expect_true(all(is.na(res$q[is.na(p)])))
  expect_true(all(is.na(res$groups[c(40, 41)])))
  expect_equal(sum(res$m_g), sum(!is.na(group) & is.finite(p)))

  # NA z behaves like NA p
  z <- c(rnorm(500), NA, Inf - Inf)
  rz <- spatial_fdr(z = z, group = rep(1:2, length.out = length(z)))
  expect_true(all(is.na(rz$p[501:502])))
  expect_false(any(rz$reject[501:502]))
  expect_true(all(is.na(rz$q[501:502])))
})
