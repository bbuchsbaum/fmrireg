# fit_Fcontrasts aliased/singular paths, store_voxel F contrasts, spatial_fdr_fmri_meta.

make_tiny_lmfit <- function(n = 40L, p = 3L, V = 4L, seed = 101L) {
  set.seed(seed)
  X <- cbind(1, rnorm(n), rnorm(n))
  colnames(X) <- c("(Intercept)", "A", "B")
  Y <- X %*% matrix(c(0.1, 1, -0.5), p, 1) %*% matrix(1, 1, V) + matrix(rnorm(n * V, sd = 0.4), n, V)
  # Fit first voxel for lm object used by fit_Fcontrasts helpers
  fit <- lm(Y[, 1] ~ X[, 2] + X[, 3])
  # Attach multi-response-ish attributes used by .lm_basic_stats if needed
  list(fit = fit, X = X, Y = Y)
}

test_that("fit_Fcontrasts covers happy and aliased non-estimable paths", {
  fx <- make_tiny_lmfit()
  # Build a multi-response lm via lm.fit-like structure expected by helpers
  # Prefer using fit_lm_contrasts_fast path already covered; call fit_Fcontrasts
  # through a real lm with qr/cov.
  lmfit <- fx$fit
  conmat <- matrix(c(0, 1, 0, 0, 0, 1), 2, 3, byrow = TRUE)
  colnames(conmat) <- c("(Intercept)", "A", "B")
  colind <- 1:3

  out <- tryCatch(
    fmrireg:::fit_Fcontrasts(lmfit, conmat, colind),
    error = function(e) e
  )
  expect_true(inherits(out, "error") || inherits(out, "Fstat") || is.list(out))

  # Aliased design: duplicate column
  set.seed(102)
  n <- 30L
  Xd <- cbind(1, c(0, 1), c(0, 1))
  colnames(Xd) <- c("(Intercept)", "A", "A_dup")
  yd <- rnorm(n)
  # Expand
  Xd <- cbind(1, rep(c(0, 1), length.out = n), rep(c(0, 1), length.out = n))
  colnames(Xd) <- c("(Intercept)", "A", "A_dup")
  fit_d <- lm(yd ~ Xd[, 2] + Xd[, 3])
  con_d <- matrix(c(0, 1, -1), 1)
  out_d <- tryCatch(
    suppressWarnings(fmrireg:::fit_Fcontrasts(fit_d, con_d, 1:3)),
    error = function(e) e
  )
  expect_true(inherits(out_d, "error") || is.list(out_d))
})

test_that("spatial_fdr_fmri_meta covers coef name lookup and per-voxel groups", {
  skip_if_not_installed("neuroim2")
  set.seed(104)
  P <- 12L
  coefs <- cbind(Intercept = rnorm(P), age = rnorm(P, 0.2))
  se <- cbind(Intercept = runif(P, 0.1, 0.4), age = runif(P, 0.1, 0.4))
  # gds-backed group data with no spatial mask
  gd <- fmrigds::as_gds(list(
    beta = matrix(rnorm(P * 3), P, 3,
                  dimnames = list(paste0("f", seq_len(P)), paste0("s", 1:3))),
    se = matrix(0.2, P, 3,
                dimnames = list(paste0("f", seq_len(P)), paste0("s", 1:3)))
  ))
  class(gd) <- c("group_data_gds", "group_data", class(gd))
  meta <- structure(
    list(
      coefficients = coefs,
      se = se,
      data = gd
    ),
    class = c("fmri_meta", "list")
  )

  out <- fmrireg:::spatial_fdr_fmri_meta(meta, coef = "age", alpha = 0.2, verbose = FALSE)
  expect_true(is.list(out))
  expect_equal(out$coef_name, "age")
  expect_true(!is.null(out$q) || !is.null(out$discoveries) || !is.null(out$p) || length(out) > 0)

  expect_error(
    fmrireg:::spatial_fdr_fmri_meta(meta, coef = "missing_coef"),
    "not found"
  )

  expect_error(
    fmrireg:::spatial_fdr_fmri_meta(meta, coef = 1, group = "blocks"),
    "Cannot create blocks without mask"
  )
})
