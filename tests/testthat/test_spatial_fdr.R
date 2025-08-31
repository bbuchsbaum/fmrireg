test_that("spatial_fdr basic pipeline computes expected group pi0 and weights", {
  skip_on_cran()

  # Two groups, 5 hypotheses each
  p <- c(0.001, 0.005, 0.02, 0.03, 0.10,  0.2, 0.5, 0.6, 0.8, 0.9)
  group <- c(rep(1, 5), rep(2, 5))

  # tau = 0.5 => group 1: tail=0 of 5 -> pi0_raw = 0
  #              group 2: tail=3 (0.6,0.8,0.9) of 5 -> 3/(0.5*5)=1.2 -> clip to 1
  res <- fmrireg:::spatial_fdr(p = p, group = group, alpha = 0.1, tau = 0.5,
                               neighbors = NULL, lambda = 0, empirical_null = FALSE,
                               min_pi0 = 0.05, verbose = FALSE)

  expect_s3_class(res, "spatial_fdr_result")
  expect_equal(as.integer(res$groups), as.integer(factor(group)))
  expect_equal(res$G, 2)
  expect_equal(res$m_g, c(5L, 5L))
  expect_length(res$pi0_raw, 2)
  expect_equal(res$pi0_raw, c(0, 1), tolerance = 1e-12)
  expect_equal(res$pi0_smooth, pmax(res$pi0_raw, 0.05))  # min_pi0 applied

  # Weights proportional to 1/pi0_smooth, but normalized to sum = #valid p
  # pi0_smooth: c(0.05, 1) -> a_g = c(20, 1)
  # Normalized weights sum to 10, so group weights become roughly c(1.90476, 0.09524)
  w <- res$weights
  expect_equal(sum(w[is.finite(p)]), sum(is.finite(p)))
  expect_true(all(w[1:5] > w[6:10]))

  # Ratio of group weights equals ratio of inverse-pi0 (up to common scale)
  w1 <- mean(w[1:5]); w2 <- mean(w[6:10])
  expect_equal(w1 / w2, (1/0.05) / (1/1.0), tolerance = 1e-6)

  # Rejection set matches rule: reject if p_i <= threshold * w_i
  thr <- res$threshold
  expect_identical(res$reject, p <= thr * w)

  # q-values equal BH q-values computed on scaled p' = p / w
  # Implement BH q-values in R for verification
  q_scaled <- p / w
  isv <- is.finite(q_scaled)
  qs <- sort(q_scaled[isv], index.return = TRUE)
  n <- sum(isv)
  qv <- rep(NA_real_, length(p))
  minv <- 1
  for (k in seq_len(n)) {
    idx <- n - k + 1
    val <- n * qs$x[idx] / idx
    if (!is.finite(val)) val <- 1
    minv <- min(minv, val)
    qv[which(isv)[qs$ix[idx]]] <- min(1, max(0, minv))
  }
  expect_equal(res$q, qv, tolerance = 1e-12)
})

test_that("spatial_fdr smoothing across neighbors matches one-step formula", {
  skip_on_cran()

  # Three groups with chain adjacency: 1 - 2 - 3
  # Make group 2 very small pi0, others large, to test smoothing pulls toward neighbors.
  p <- c(0.001, 0.01, 0.02,  0.9, 0.95, 0.99,  0.9, 0.95, 0.99)
  group <- c(2,2,2, 1,1,1, 3,3,3)  # note shuffled order to test compression

  # tau = 0.5 -> pi0_raw: g1=1, g2=0, g3=1 (after clipping to [0,1])
  neighbors <- list(
    c(2L),      # neighbors of group 1
    c(1L, 3L),  # group 2
    c(2L)       # group 3
  )

  res0 <- fmrireg:::spatial_fdr(p = p, group = group, alpha = 0.1, tau = 0.5,
                                neighbors = neighbors, lambda = 0, min_pi0 = 0.05,
                                empirical_null = FALSE, verbose = FALSE)
  res1 <- fmrireg:::spatial_fdr(p = p, group = group, alpha = 0.1, tau = 0.5,
                                neighbors = neighbors, lambda = 1.0, min_pi0 = 0.05,
                                empirical_null = FALSE, verbose = FALSE)

  # Without smoothing, pi0_smooth equals pi0_raw with min_pi0 applied
  expect_equal(res0$pi0_smooth, pmax(res0$pi0_raw, 0.05), tolerance = 1e-12)

  # With lambda=1 and single iteration, formula (before min clamp):
  # nxt[g] = (pi0_raw[g] + mean(pi0_raw[neighbors[g]])) / 2
  cur_raw <- res0$pi0_raw
  expected <- numeric(3)
  expected[1] <- (cur_raw[1] + mean(cur_raw[2])) / 2
  expected[2] <- (cur_raw[2] + mean(c(cur_raw[1], cur_raw[3]))) / 2
  expected[3] <- (cur_raw[3] + mean(cur_raw[2])) / 2
  # spatial_fdr applies min_pi0 clamp after smoothing
  expected <- pmax(expected, 0.05)
  expect_equal(as.numeric(res1$pi0_smooth), expected, tolerance = 1e-12)
})

test_that("spatial_fdr handles NA inputs and invalid parameters", {
  skip_on_cran()

  p <- c(0.01, NA, 0.2, 0.8)
  group <- c(1, 1, 2, NA)

  # NA p and NA group should not cause errors; NA p yields NA q and no rejection
  res <- fmrireg:::spatial_fdr(p = p, group = group, alpha = 0.1, tau = 0.5,
                               neighbors = NULL, lambda = 0, empirical_null = FALSE,
                               min_pi0 = 0.05, verbose = FALSE)
  expect_false(res$reject[2])
  expect_true(is.na(res$q[2]))

  # Parameter validation
  expect_error(fmrireg:::spatial_fdr(p = p, group = group[-1], alpha = 0.1),
               "group length must match")
  expect_error(fmrireg:::spatial_fdr(p = p, group = group, alpha = 1.0),
               "alpha must be in")
  expect_error(fmrireg:::spatial_fdr(p = p, group = group, alpha = 0.1, tau = 0),
               "tau must be in")
  expect_error(fmrireg:::spatial_fdr(p = p, group = group, alpha = 0.1, tau = 0.5,
                                     neighbors = list(1L)),
               "neighbors must be a list of length")
})

test_that("spatial_fdr with z input and empirical_null fallback uses theoretical null", {
  skip_on_cran()

  # Fewer than 10 valid z-values triggers theoretical null (mu0=0, sigma0=1)
  z <- c(-0.5, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 2.5)
  group <- c(1,1,1,1,2,2,2,2)

  res <- fmrireg:::spatial_fdr(z = z, group = group, alpha = 0.1, tau = 0.5,
                               neighbors = NULL, empirical_null = TRUE,
                               verbose = FALSE)

  expect_equal(res$mu0, 0)
  expect_equal(res$sigma0, 1)
  # p computed from two-sided normal with mu=0, sigma=1
  expect_equal(res$p, 2 * pnorm(abs(z), lower.tail = FALSE), tolerance = 1e-12)
})
