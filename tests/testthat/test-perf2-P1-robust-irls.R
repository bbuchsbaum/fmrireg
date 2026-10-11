# Performance pass 2, item P1-4: robust IRLS without T x V weighted copies.
#
# Reference 1: the IRLS implementation as it stood before the change, copied
# verbatim below (evaluated in the package namespace).
# Reference 2: stats::lm(weights = ) with the final robust weights.

.p14_old_fitter <- function(initial_glm_ctx,
                                    cfg_robust_options,
                                    X_orig_for_resid,
                                    sigma_fixed = NULL,
                                    obs_weights = NULL) {
  
  # Validate inputs
  stopifnot(is.glm_context(initial_glm_ctx))
  stopifnot(is.list(cfg_robust_options))
  stopifnot(!is.null(X_orig_for_resid))
  
  # Extract options
  psi_type <- cfg_robust_options$type
  if (isFALSE(psi_type)) {
    stop("robust_iterative_fitter called with robust type = FALSE")
  }
  
  k_huber <- cfg_robust_options$k_huber
  c_tukey <- cfg_robust_options$c_tukey
  max_iter <- max(1L, as.integer(cfg_robust_options$max_iter %||% 2L))
  tol <- as.numeric(cfg_robust_options$tol %||% 1e-4)
  min_iter <- as.integer(cfg_robust_options$min_iter %||% 1L)
  if (!is.finite(tol) || tol < 0) tol <- 1e-4
  if (!is.finite(min_iter) || min_iter < 1L) min_iter <- 1L
  scope_opt <- cfg_robust_options$scale_scope
  if (is.null(scope_opt)) scope_opt <- "run"
  scale_scope <- tolower(as.character(scope_opt))
  if (identical(scale_scope, "local")) scale_scope <- "voxel"
  if (!(scale_scope %in% c("run", "global", "voxel"))) {
    scale_scope <- "run"
  }
  
  # Ensure matrices
  if (!is.matrix(X_orig_for_resid)) X_orig_for_resid <- as.matrix(X_orig_for_resid)
  Y <- initial_glm_ctx$Y
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  
  # Check for NAs
  if (anyNA(Y) || anyNA(X_orig_for_resid)) {
    stop("NA values detected in 'X_orig_for_resid' or 'Y' for robust_iterative_fitter")
  }
  
  # Initial fit from context (could be OLS or whitened OLS)
  initial_fit <- solve_glm_core(initial_glm_ctx, return_fitted = FALSE)
  betas <- initial_fit$betas
  
  # Initialize final values
  XtWXi_final <- initial_glm_ctx$proj$XtXinv
  dfres <- initial_glm_ctx$proj$dfres

  compute_scale <- function(resid_mat, scope, sigma_fixed_value = NULL) {
    eps <- .Machine$double.eps
    abs_resid <- abs(resid_mat)

    if (!is.null(sigma_fixed_value)) {
      sigma_vec <- as.numeric(sigma_fixed_value)
      if (length(sigma_vec) == 1L) {
        sigma_vec <- rep(sigma_vec, ncol(abs_resid))
      } else if (length(sigma_vec) != ncol(abs_resid)) {
        stop("sigma_fixed must have length 1 or match number of response columns")
      }
      sigma_vec[!is.finite(sigma_vec) | sigma_vec <= eps] <- eps
    } else if (identical(scope, "voxel")) {
      sigma_vec <- 1.4826 * matrixStats::colMedians(abs_resid)
      sigma_vec[!is.finite(sigma_vec) | sigma_vec <= eps] <- eps
    } else {
      row_med <- matrixStats::rowMedians(abs_resid)
      sigma_scalar <- 1.4826 * stats::median(row_med)
      if (!is.finite(sigma_scalar) || sigma_scalar <= eps) {
        sigma_scalar <- eps
      }
      sigma_vec <- rep(sigma_scalar, ncol(abs_resid))
    }

    u <- matrixStats::rowMedians(sweep(abs_resid, 2, sigma_vec, "/"))
    sigma_scalar <- stats::median(sigma_vec)
    if (!is.finite(sigma_scalar) || sigma_scalar <= eps) {
      sigma_scalar <- eps
    }

    list(u = u, sigma_scalar = sigma_scalar, sigma_components = sigma_vec)
  }
  
  psi_weights <- function(u) {
    switch(psi_type,
           huber = pmin(1, k_huber / abs(u)),
           bisquare = ifelse(abs(u) <= c_tukey, (1 - (u / c_tukey)^2)^2, 0))
  }

  n_obs <- nrow(Y)
  active_prior <- if (is.null(obs_weights)) {
    rep(TRUE, n_obs)
  } else {
    as.numeric(obs_weights) > 0
  }

  # IRLS iterations. Everything returned below comes from the last weighted
  # solve: the weights it used, its betas, its (X'WX)^-1 and its residuals.
  w <- rep(1, n_obs)
  scale_info <- NULL
  proj_w <- initial_glm_ctx$proj
  rss_w <- initial_fit$rss
  for (it in seq_len(max_iter)) {
    betas_prev <- betas

    # Calculate residuals using original (unweighted) X
    resid <- Y - X_orig_for_resid %*% betas

    scale_info <- compute_scale(resid, scale_scope, sigma_fixed)
    w <- psi_weights(scale_info$u)

    # Apply sqrt of weights for weighted least squares
    sqrtw <- sqrt(w)
    Xw <- X_orig_for_resid * sqrtw
    Yw <- Y * sqrtw

    proj_w <- .fast_preproject(Xw)
    fit_w <- solve_glm_core(
      glm_context(X = Xw, Y = Yw, proj = proj_w, robust_weights = w),
      return_fitted = FALSE
    )

    betas <- fit_w$betas
    XtWXi_final <- proj_w$XtXinv
    rss_w <- fit_w$rss

    # Convergence check on coefficient updates.
    if (it >= min_iter) {
      scale <- max(1, max(abs(betas_prev)))
      delta_beta <- max(abs(betas - betas_prev)) / scale
      if (is.finite(delta_beta) && delta_beta <= tol) {
        break
      }
    }
  }

  # Residual df of the WLS fit: a zero-weight row carries no information, as
  # in stats::lm.wfit().
  dfres <- sum(w > 0 & active_prior) - proj_w$rank
  sigma2 <- if (dfres > 0) rss_w / dfres else rep(NA_real_, ncol(Y))

  list(
    betas_robust = betas,
    XtWXi_final = XtWXi_final,
    sigma2 = sigma2,
    rss = rss_w,
    sigma_robust_scale_final = scale_info$sigma_scalar,
    sigma_robust_scale_components = scale_info$sigma_components,
    robust_weights_final = w,
    dfres = dfres,
    rank = proj_w$rank,
    scale_scope = scale_scope
  )
}

environment(.p14_old_fitter) <- asNamespace("fmrireg")

.p14_ctx <- function(X, Y) {
  fmrireg:::glm_context(X = X, Y = Y, proj = fmrireg:::.fast_preproject(X))
}

.p14_problem <- function(n = 80L, V = 7L, seed = 2L) {
  set.seed(seed)
  X <- cbind(Intercept = 1, a = rnorm(n), b = sin(seq_len(n) / 5))
  Y <- X %*% matrix(rnorm(3 * V), 3, V) + matrix(rnorm(n * V), n, V)
  # gross outliers in a few volumes
  Y[c(5, 17, 40), ] <- Y[c(5, 17, 40), ] + 8
  Y[60, 2] <- Y[60, 2] - 12
  list(X = X, Y = Y)
}

.p14_configs <- list(
  list(type = "huber", k_huber = 1.345, c_tukey = 4.685, max_iter = 6L,
       scale_scope = "voxel"),
  list(type = "huber", k_huber = 1.345, c_tukey = 4.685, max_iter = 3L,
       scale_scope = "run"),
  list(type = "bisquare", k_huber = 1.345, c_tukey = 4.685, max_iter = 8L,
       scale_scope = "global"),
  list(type = "bisquare", k_huber = 1.345, c_tukey = 4.685, max_iter = 5L,
       scale_scope = "voxel")
)

test_that("robust IRLS output equals the previous implementation", {
  pr <- .p14_problem()
  for (cfg in .p14_configs) {
    for (sf in list(NULL, 0.9)) {
      new <- fmrireg:::robust_iterative_fitter(.p14_ctx(pr$X, pr$Y), cfg, pr$X,
                                               sigma_fixed = sf)
      old <- .p14_old_fitter(.p14_ctx(pr$X, pr$Y), cfg, pr$X, sigma_fixed = sf)
      info <- paste(cfg$type, cfg$scale_scope, is.null(sf))
      expect_equal(new$robust_weights_final, old$robust_weights_final,
                   tolerance = 1e-10, info = info)
      expect_equal(new$betas_robust, old$betas_robust, tolerance = 1e-10,
                   info = info)
      expect_equal(new$XtWXi_final, old$XtWXi_final, tolerance = 1e-10,
                   info = info)
      expect_equal(new$rss, old$rss, tolerance = 1e-10, info = info)
      expect_equal(new$sigma2, old$sigma2, tolerance = 1e-10, info = info)
      expect_identical(new$dfres, old$dfres)
      expect_identical(new$rank, old$rank)
      expect_equal(new$sigma_robust_scale_final, old$sigma_robust_scale_final,
                   tolerance = 1e-10)
    }
  }
  # prior observation weights (zero rows) and a rank-deficient design
  X2 <- cbind(pr$X, pr$X[, 2] * 2)
  colnames(X2)[4] <- "alias"
  ow <- rep(1, nrow(pr$X))
  ow[c(3, 9)] <- 0
  cfg <- .p14_configs[[4]]
  new <- suppressWarnings(fmrireg:::robust_iterative_fitter(
    .p14_ctx(X2, pr$Y), cfg, X2, obs_weights = ow))
  old <- suppressWarnings(.p14_old_fitter(.p14_ctx(X2, pr$Y), cfg, X2,
                                          obs_weights = ow))
  expect_equal(new$betas_robust, old$betas_robust, tolerance = 1e-10)
  expect_identical(attributes(new$betas_robust), attributes(old$betas_robust))
  expect_equal(new$sigma2, old$sigma2, tolerance = 1e-10)
  expect_identical(new$dfres, old$dfres)
})

test_that("robust IRLS betas and scale are the WLS fit with the final weights", {
  pr <- .p14_problem(V = 4L, seed = 8L)
  cfg <- .p14_configs[[1]]
  fit <- fmrireg:::robust_iterative_fitter(.p14_ctx(pr$X, pr$Y), cfg, pr$X)
  w <- fit$robust_weights_final
  for (v in seq_len(ncol(pr$Y))) {
    ref <- lm(pr$Y[, v] ~ pr$X - 1, weights = w)
    expect_equal(unname(fit$betas_robust[, v]), unname(coef(ref)),
                 tolerance = 1e-10)
    expect_equal(fit$sigma2[v], summary(ref)$sigma^2, tolerance = 1e-10)
    expect_equal(as.numeric(fit$XtWXi_final * fit$sigma2[v]), as.numeric(vcov(ref)),
                 tolerance = 1e-10)
  }
  expect_identical(rownames(fit$betas_robust), colnames(pr$X))
})

test_that("robust IRLS never forms the weighted response", {
  # The weighted solve is crossprod-style: one solve_glm_core() call (the
  # initial fit) and no T x V weighted copy of Y per iteration.
  pr <- .p14_problem(n = 200L, V = 50L)
  cfg <- .p14_configs[[1]]
  cfg$max_iter <- 4L
  cfg$tol <- 0
  calls <- 0L
  orig_solve <- fmrireg:::solve_glm_core
  testthat::local_mocked_bindings(
    solve_glm_core = function(...) { calls <<- calls + 1L; orig_solve(...) }
  )
  fit <- fmrireg:::robust_iterative_fitter(.p14_ctx(pr$X, pr$Y), cfg, pr$X)
  expect_equal(calls, 1L)

  # Allocation count (deterministic): T x V double blocks allocated per
  # iteration (including the initial fit and final RSS, amortised). The
  # previous implementation allocated 8.75; this one allocates 2.
  skip_if_not(capabilities("profmem"))
  ctx <- .p14_ctx(pr$X, pr$Y)
  big <- 8 * nrow(pr$Y) * ncol(pr$Y)
  logfile <- tempfile()
  utils::Rprofmem(logfile, threshold = big)
  fit <- fmrireg:::robust_iterative_fitter(ctx, cfg, pr$X)
  utils::Rprofmem(NULL)
  lines <- readLines(logfile)
  sizes <- suppressWarnings(as.numeric(sub(":.*", "", lines)))
  n_big <- sum(sizes >= big, na.rm = TRUE)
  unlink(logfile)
  expect_lte(n_big / cfg$max_iter, 3)
})
