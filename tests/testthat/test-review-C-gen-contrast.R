# Review 2026-10, stream C, C3: generate_main_effect_contrast() and
# generate_interaction_contrast() rows must follow the rows of `des`
# (expand.grid order: first factor varies fastest), which is also the column
# order fmridesign uses for crossed terms such as hrf(Time, Cond).

.cell_means <- function(des, a, b, ab) {
  # mu = a[Time] + b[Cond] + ab[Time, Cond]
  ti <- as.integer(des$Time)
  ci <- as.integer(des$Cond)
  a[ti] + b[ci] + ab[cbind(ti, ci)]
}

test_that("C3: main-effect contrast applied to expand.grid cell means isolates that factor", {
  des <- expand.grid(Time = factor(1:4), Cond = factor(c("face", "scene")))
  a <- c(0, 1, 3, 7)
  b <- c(0, 10)
  mu <- .cell_means(des, a, b, matrix(0, 4, 2))

  M <- generate_main_effect_contrast(des, "Time")
  # Successive Time differences, summed over the two Cond levels.
  expect_equal(drop(crossprod(M, mu)) / 2, -diff(a))
  # Cond main effect does not leak into the Time contrast.
  mu_cond_only <- .cell_means(des, rep(0, 4), b, matrix(0, 4, 2))
  expect_equal(drop(crossprod(M, mu_cond_only)), rep(0, 3))

  MC <- generate_main_effect_contrast(des, "Cond")
  expect_equal(drop(crossprod(MC, mu)) / 4, -diff(b))
  expect_equal(drop(crossprod(MC, .cell_means(des, a, c(0, 0), matrix(0, 4, 2)))), 0)
})

test_that("C3: interaction contrast is zero for additive means and recovers the interaction", {
  des <- expand.grid(Time = factor(1:4), Cond = factor(c("face", "scene")))
  a <- c(0, 1, 3, 7)
  b <- c(0, 10)
  I <- generate_interaction_contrast(des, c("Time", "Cond"))
  expect_equal(drop(crossprod(I, .cell_means(des, a, b, matrix(0, 4, 2)))), rep(0, 3))

  ab <- matrix(0, 4, 2)
  ab[3, 2] <- 5
  # Difference-of-differences: (d_{t+1} - d_t) where d_t = mu[t, face] - mu[t, scene]
  d <- ab[, 1] - ab[, 2]
  expect_equal(drop(crossprod(I, .cell_means(des, a, b, ab))), -diff(d))
})

test_that("C3: contrasts follow the rows of a permuted design", {
  des <- expand.grid(Time = factor(1:3), Cond = factor(c("x", "y")))
  perm <- c(4, 1, 6, 2, 5, 3)
  des_p <- des[perm, ]
  M <- generate_main_effect_contrast(des, "Time")
  Mp <- generate_main_effect_contrast(des_p, "Time")
  expect_equal(unname(Mp), unname(M[perm, , drop = FALSE]))
})

test_that("C3: contrast rows match fmridesign's crossed-term column order", {
  set.seed(3)
  levs <- expand.grid(Time = factor(1:3), Cond = factor(c("face", "scene")))
  etab <- levs[rep(seq_len(nrow(levs)), 3), ]
  etab$onset <- seq(5, by = 12, length.out = nrow(etab))
  etab$run <- 1L
  n <- max(etab$onset) + 30
  a <- c(0, 2, 5)
  b <- c(0, 4)
  ab <- matrix(0, 3, 2); ab[2, 2] <- 3
  cell_beta <- .cell_means(levs, a, b, ab)

  sf <- fmrihrf::sampling_frame(blocklens = n, TR = 1)
  emod <- event_model(onset ~ hrf(Time, Cond), data = etab, block = ~ run,
                      sampling_frame = sf)
  X <- as.matrix(design_matrix(emod))
  # Column j of the design corresponds to row j of expand.grid(levels).
  expect_equal(ncol(X), nrow(levs))
  expect_equal(colnames(X),
               paste0("Time_Cond_Time.", levs$Time, "_Cond.", levs$Cond))

  y <- drop(X %*% cell_beta) + 1
  beta_hat <- coef(lm(y ~ X))[-1]
  M <- generate_main_effect_contrast(levs, "Time")
  I <- generate_interaction_contrast(levs, c("Time", "Cond"))
  expect_equal(unname(drop(crossprod(M, beta_hat))) / 2,
               -diff(a) - diff(rowMeans(ab)), tolerance = 1e-8)
  expect_equal(unname(drop(crossprod(I, beta_hat))), -diff(ab[, 1] - ab[, 2]),
               tolerance = 1e-8)
})
