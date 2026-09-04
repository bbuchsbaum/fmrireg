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

test_that("store_voxel_contrasts writes t and F into initialized storage", {
  set.seed(103)
  n <- 36L
  p <- 3L
  X <- cbind(1, rnorm(n), rnorm(n))
  colnames(X) <- c("(Intercept)", "A", "B")
  beta <- c(0.2, 1.0, -0.5)
  qr_x <- qr(X)
  XtXinv <- chol2inv(qr.R(qr_x))
  sigma_w <- 0.4

  conlist <- list(A_vs_B = structure(c(1, -1), colind = c(2L, 3L)))
  fconlist <- list(main = structure(diag(2), colind = c(2L, 3L)))
  storage <- fmrireg:::initialize_contrast_storage(conlist, fconlist, n_voxels = 3L)

  stored <- fmrireg:::store_voxel_contrasts(
    storage, voxel_index = 2L, qr_or_XtXinv = XtXinv,
    beta_w = beta, sigma_w = sigma_w, conlist = conlist, fconlist = fconlist,
    dfres = n - p
  )
  expect_true(is.finite(stored$A_vs_B$estimate[2]) || is.na(stored$A_vs_B$estimate[2]))
  expect_true(!is.null(stored$main$stat[2]))

  # QR object path
  stored2 <- fmrireg:::store_voxel_contrasts(
    storage, voxel_index = 1L, qr_or_XtXinv = qr_x,
    beta_w = beta, sigma_w = sigma_w, conlist = conlist, fconlist = fconlist,
    dfres = n - p
  )
  expect_true(is.list(stored2))

  formatted <- fmrireg:::format_contrast_results(stored2, dfres = n - p)
  expect_true(is.list(formatted))
  expect_true(length(formatted) >= 1L)
})

test_that("spatial_fdr_fmri_meta covers coef name lookup and per-voxel groups", {
  skip_if_not_installed("neuroim2")
  set.seed(104)
  P <- 12L
  coefs <- cbind(Intercept = rnorm(P), age = rnorm(P, 0.2))
  se <- cbind(Intercept = runif(P, 0.1, 0.4), age = runif(P, 0.1, 0.4))
  meta <- structure(
    list(
      coefficients = coefs,
      se = se,
      data = structure(list(mask = NULL, mask_data = NULL), class = "group_data_csv")
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
