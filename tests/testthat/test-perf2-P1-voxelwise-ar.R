# Performance pass 2, item P1-2: voxelwise AR statistics are computed for all
# voxels at once (not one .fast_t_contrast / .fast_F_contrast call per voxel),
# and the parallel voxel path ships each worker its own columns only.

# ---- reference: the per-voxel loops as they stood before the change ---------

.p12_old_beta_stats <- function(Betas, XtXinv_list, sigma, dfres, varnames,
                                robust_weights_list = NULL) {
  V <- ncol(Betas)
  p <- nrow(Betas)
  est <- se <- tt <- pr <- matrix(NA_real_, V, p)
  for (v in seq_len(V)) {
    M <- XtXinv_list[[v]]
    if (is.null(M) || !is.matrix(M) || any(dim(M) != c(p, p)) ||
        !is.finite(sigma[v])) next
    rw <- if (!is.null(robust_weights_list)) robust_weights_list[[v]] else NULL
    df_eff <- if (is.null(rw)) dfres else
      fmrireg:::.robust_effective_df(rw, attr(M, "rank") %||% p)
    se_vec <- sqrt(diag(M)) * sigma[v]
    t_vec <- ifelse(abs(se_vec) < .Machine$double.eps^0.5, 0, Betas[, v] / se_vec)
    est[v, ] <- Betas[, v]
    se[v, ] <- se_vec
    tt[v, ] <- t_vec
    pr[v, ] <- 2 * pt(-abs(t_vec), df_eff)
  }
  list(estimate = est, se = se, stat = tt, prob = pr)
}

.p12_old_contrasts <- function(Betas, sigma2, XtXinv_list, conlist, fconlist,
                               dfres, robust_weights_list = NULL, ar_order = 0) {
  p <- nrow(Betas)
  V <- ncol(Betas)
  out <- list()
  for (nm in names(conlist)) {
    l <- conlist[[nm]]
    full_l <- matrix(0, 1, p)
    full_l[, attr(l, "colind")] <- l
    res <- matrix(NA_real_, V, 5)
    for (v in seq_len(V)) {
      if (is.null(XtXinv_list[[v]]) || !is.finite(sigma2[v])) next
      rw <- if (!is.null(robust_weights_list)) robust_weights_list[[v]] else NULL
      r <- fmrireg:::.fast_t_contrast(Betas[, v, drop = FALSE], sigma2[v],
                                      XtXinv_list[[v]], full_l, dfres, rw, ar_order)
      res[v, ] <- c(r$estimate, r$se, r$stat, r$prob, r$sigma)
    }
    out[[nm]] <- res
  }
  for (nm in names(fconlist)) {
    L <- fconlist[[nm]]
    colind <- attr(L, "colind")
    Cw <- fmrireg:::.orient_fcontrast(L, colind)
    full_L <- matrix(0, nrow(Cw), p)
    full_L[, colind] <- Cw
    res <- matrix(NA_real_, V, 4)
    for (v in seq_len(V)) {
      if (is.null(XtXinv_list[[v]]) || !is.finite(sigma2[v])) next
      rw <- if (!is.null(robust_weights_list)) robust_weights_list[[v]] else NULL
      r <- suppressWarnings(fmrireg:::.fast_F_contrast(
        Betas[, v, drop = FALSE], sigma2[v], XtXinv_list[[v]], full_L, dfres,
        rw, ar_order))
      res[v, ] <- c(r$estimate, r$se, r$stat, r$prob)
    }
    out[[nm]] <- res
  }
  out
}

.p12_problem <- function(V = 12L, n = 40L, seed = 5L) {
  set.seed(seed)
  X <- cbind(1, rnorm(n), rnorm(n), rnorm(n))
  colnames(X) <- c("Intercept", "A", "B", "C")
  p <- ncol(X)
  Betas <- matrix(rnorm(p * V), p, V)
  XtXinv_list <- lapply(seq_len(V), function(v) {
    Xv <- X * (1 + 0.2 * v / V) + 0.05 * v * sin(seq_len(n))
    M <- solve(crossprod(Xv))
    attr(M, "rank") <- p
    M
  })
  sigma <- runif(V, 0.5, 2)
  # degenerate voxels: missing covariance, non-finite scale
  XtXinv_list[2] <- list(NULL)
  sigma[3] <- NA_real_
  rw <- lapply(seq_len(V), function(v) runif(n, 0.2, 1))
  tcon <- list(
    A = structure(1, colind = 2L),
    AvB = structure(c(1, -1), colind = c(2L, 3L))
  )
  fcon <- list(
    AB = structure(diag(2), colind = c(2L, 3L)),
    one = structure(matrix(c(0, 1, -1), 1), colind = c(2L, 3L, 4L)),
    ABC = structure(cbind(c(1, -1, 0), c(0, 1, -1)), colind = c(2L, 3L, 4L)),
    # rank-deficient hypothesis: L (X'X)^-1 L' is singular for every voxel
    sing = structure(rbind(c(1, 0), c(2, 0)), colind = c(2L, 3L))
  )
  list(X = X, Betas = Betas, XtXinv_list = XtXinv_list, sigma = sigma,
       rw = rw, tcon = tcon, fcon = fcon, n = n, p = p, V = V)
}

test_that("voxelwise beta statistics equal the per-voxel reference", {
  pr <- .p12_problem()
  for (rw in list(NULL, pr$rw)) {
    got <- fmrireg:::beta_stats_matrix_voxelwise(
      pr$Betas, pr$XtXinv_list, pr$sigma, pr$n - pr$p, colnames(pr$X),
      robust_weights_list = rw, ar_order = 1L
    )$data[[1]]
    ref <- .p12_old_beta_stats(pr$Betas, pr$XtXinv_list, pr$sigma,
                               pr$n - pr$p, colnames(pr$X), rw)
    for (fld in c("estimate", "se", "stat", "prob")) {
      expect_equal(unname(got[[fld]][[1]]), ref[[fld]], tolerance = 1e-12,
                   info = fld)
    }
    expect_equal(colnames(got$estimate[[1]]), colnames(pr$X))
  }
})

test_that("voxelwise t and F contrasts equal the per-voxel reference", {
  pr <- .p12_problem()
  for (rw in list(NULL, pr$rw)) {
    for (ar_order in c(0L, 1L)) {
      got <- suppressWarnings(fmrireg:::fit_lm_contrasts_voxelwise(
        pr$Betas, pr$sigma^2, pr$XtXinv_list, pr$tcon, pr$fcon,
        pr$n - pr$p, robust_weights_list = rw, ar_order = ar_order
      ))
      ref <- .p12_old_contrasts(pr$Betas, pr$sigma^2, pr$XtXinv_list,
                                pr$tcon, pr$fcon, pr$n - pr$p, rw, ar_order)
      for (nm in names(pr$tcon)) {
        d <- got[[nm]]$data[[1]]
        expect_equal(cbind(d$estimate, d$se, d$stat, d$prob, d$sigma),
                     ref[[nm]], tolerance = 1e-10, info = nm)
        expect_equal(got[[nm]]$type, "contrast")
      }
      for (nm in names(pr$fcon)) {
        d <- got[[nm]]$data[[1]]
        expect_equal(cbind(d$estimate, d$se, d$stat, d$prob),
                     ref[[nm]], tolerance = 1e-10, info = nm)
        expect_equal(got[[nm]]$type, "Fcontrast")
      }
    }
  }
  # A singular hypothesis still warns, as the per-voxel solve() did.
  expect_warning(
    fmrireg:::fit_lm_contrasts_voxelwise(
      pr$Betas, pr$sigma^2, pr$XtXinv_list, list(), pr$fcon["sing"],
      pr$n - pr$p),
    "Singular"
  )
})

test_that("voxelwise AR fit equals lm on each voxel's whitened data", {
  set.seed(21)
  n <- 60L
  V <- 5L
  etab <- data.frame(onset = seq(2, 110, by = 12),
                     cond = factor(rep(c("A", "B", "C"), length.out = 10)),
                     run = 1)
  Y <- vapply(seq_len(V), function(v) {
    as.numeric(arima.sim(list(ar = 0.2 + 0.1 * v), n)) + 5
  }, numeric(n))
  dset <- matrix_frame(Y, TR = 2, run_length = n, event_table = etab)
  con <- contrast_set(
    pair_contrast(~ cond == "A", ~ cond == "B", name = "AvB"),
    oneway_contrast(~ cond, name = "main")
  )
  for (par in c(FALSE, TRUE)) {
    fit <- fmri_lm(onset ~ hrf(cond, contrasts = con), block = ~run,
                   dataset = dset, use_fast_path = FALSE, ar_voxelwise = TRUE,
                   cor_struct = "ar1", parallel_voxels = par)
    X <- as.matrix(design_matrix(fit$model))
    b <- fit$result$betas$data[[1]]
    tab_c <- fit$result$contrasts
    avb <- tab_c$data[[which(tab_c$name == "AvB")]]
    main <- tab_c$data[[which(tab_c$name == "main")]]
    for (v in seq_len(V)) {
      e <- residuals(lm(Y[, v] ~ X - 1))
      phi <- fmrireg:::estimate_ar_parameters(e, 1L)
      w <- fmrireg:::ar_whiten_transform(X, matrix(Y[, v]), phi, FALSE)
      fv <- lm(w$Y[, 1] ~ w$X - 1)
      sm <- summary(fv)$coefficients
      expect_equal(unname(b$estimate[[1]][v, ]), unname(coef(fv)),
                   tolerance = 1e-8)
      expect_equal(unname(b$se[[1]][v, ]), unname(sm[, 2]), tolerance = 1e-8)
      expect_equal(unname(b$stat[[1]][v, ]), unname(sm[, 3]), tolerance = 1e-8)
      # AvB t-contrast
      cn <- colnames(X)
      l <- numeric(ncol(X))
      l[grep("cond.A$|condA|cond_A", cn)] <- 1
      l[grep("cond.B$|condB|cond_B", cn)] <- -1
      est <- sum(l * coef(fv))
      se <- sqrt(drop(t(l) %*% vcov(fv) %*% l))
      expect_equal(avb$estimate[v], est, tolerance = 1e-8)
      expect_equal(avb$stat[v], est / se, tolerance = 1e-8)
      # main-effect F (equal condition effects) equals the nested-model F
      cols <- grep("^cond", cn)
      reduced <- lm(w$Y[, 1] ~ cbind(rowSums(w$X[, cols]), w$X[, -cols]) - 1)
      expect_equal(main$stat[v], anova(reduced, fv)$F[2], tolerance = 1e-8)
    }
  }
})

test_that("voxelwise contrasts make O(contrasts) calls, not O(voxels)", {
  pr <- .p12_problem(V = 40L)
  calls <- 0L
  orig_t <- fmrireg:::.fast_t_contrast
  orig_f <- fmrireg:::.fast_F_contrast
  testthat::local_mocked_bindings(
    .fast_t_contrast = function(...) { calls <<- calls + 1L; orig_t(...) },
    .fast_F_contrast = function(...) { calls <<- calls + 1L; orig_f(...) }
  )
  fmrireg:::fit_lm_contrasts_voxelwise(
    pr$Betas, pr$sigma^2, pr$XtXinv_list, pr$tcon, pr$fcon[1:3],
    pr$n - pr$p, robust_weights_list = pr$rw, ar_order = 1L
  )
  expect_lte(calls, length(pr$tcon) + 3L)
})

test_that("parallel voxelwise AR ships each worker its own columns only", {
  skip_if_not_installed("future.apply")
  set.seed(3)
  n <- 40L
  V <- 25L
  etab <- data.frame(onset = c(2, 20, 40, 60), cond = factor(c("A", "B", "A", "B")),
                     run = 1)
  Y <- matrix(rnorm(n * V), n, V)
  dset <- matrix_frame(Y, TR = 2, run_length = n, event_table = etab)
  seen <- list()
  testthat::local_mocked_bindings(
    future_lapply = function(X, FUN, ..., future.seed = FALSE) {
      seen$X <<- X
      seen$FUN <<- FUN
      seen$dots <<- list(...)
      lapply(X, FUN, ...)
    },
    .package = "future.apply"
  )
  fit <- fmri_lm(onset ~ hrf(cond), block = ~run, dataset = dset,
                 use_fast_path = FALSE, ar_voxelwise = TRUE,
                 cor_struct = "ar1", parallel_voxels = TRUE)
  expect_false(is.null(seen$FUN))
  # The worker function does not close over the run's data matrix ...
  env <- environment(seen$FUN)
  expect_false(exists("Y_run", envir = env, inherits = FALSE))
  big <- function(obj) is.matrix(obj) && ncol(obj) == V && nrow(obj) == n
  expect_false(any(vapply(as.list(env), big, logical(1))) &&
                 !isNamespace(env))
  # ... no shared argument carries it ...
  expect_false(any(vapply(seen$dots, big, logical(1))))
  # ... and each element holds exactly one voxel's series.
  expect_length(seen$X, V)
  expect_true(all(vapply(seen$X, function(x) length(x) == n, logical(1))))
  expect_equal(nrow(coef(fit)), V)
})
