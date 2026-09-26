# Regression tests for the inference reported by the time-sketched GLM
# (engine = "latent_sketch"): data-scale fields, calibrated sketch-and-solve
# standard errors, OLS-exact IHS, and per-cluster (by_cluster) solves.

.sketch_fit <- function(dset, method, m, ...) {
  suppressWarnings(fmri_lm(
    onsets ~ hrf(condition), block = ~run, dataset = dset,
    engine = "latent_sketch",
    lowrank = lowrank_control(time_sketch = list(method = method, m = m, ...))
  ))
}

.beta_payload <- function(fit) {
  b <- fit$result$betas$data[[1]]
  list(est = as.matrix(b$estimate[[1]]), se = as.matrix(b$se[[1]]),
       prob = as.matrix(b$prob[[1]]))
}

test_that("SRHT is a normalised sketch with an exact adjoint", {
  # The Hadamard transform is unnormalised; the plan's 1 / sqrt(m) scale must
  # make E||Sr||^2 = ||r||^2. Under the former sqrt(T / m) scale ||Sr||^2 was
  # about T ||r||^2, which put sigma2, rss and cov.unscaled T-fold off.
  set.seed(301)
  Tlen <- 128L; m <- 32L
  plan <- make_srht_plan(Tlen, m)
  S <- srht_apply(diag(Tlen), plan)
  # T a power of two: S S' = (T / m) I exactly.
  expect_equal(tcrossprod(S), diag(Tlen / m, m), tolerance = 1e-12)
  B <- matrix(rnorm(m * 3), m)
  expect_equal(fmrireg:::srht_adjoint(B, plan), crossprod(S, B), tolerance = 1e-12)

  # Zero-padded length (T = 150 -> 256): still normalised in expectation.
  Tlen <- 150L
  r <- rnorm(Tlen)
  ratio <- vapply(1:400, function(i) {
    sum(srht_apply(matrix(r), make_srht_plan(Tlen, 30L))^2) / sum(r^2)
  }, numeric(1))
  expect_lt(abs(mean(ratio) - 1), 0.05)
})

test_that("sketched fits report residual variance and covariance on the data scale", {
  dset <- sketch_matrix_dataset(run_length = 60L, nvox = 300L, seed = 11L)
  for (method in c("srht", "gaussian", "countsketch", "ihs")) {
    set.seed(12)
    fit <- .sketch_fit(dset, method, m = 40L)
    # iid N(0, 1) noise around a baseline of 100: sigma^2 = 1.
    expect_lt(abs(mean(fit$sigma2) - 1), 0.1, label = method)
    expect_equal(fit$result$resvar, fit$sigma2, label = method)
    expect_equal(fit$result$sigma, sqrt(fit$sigma2), label = method)
    # rss is the residual sum of squares of the fitted rows, whose
    # expectation is sigma^2 * kappa (kappa = tr(PSS'), the residual count
    # for IHS); it is not sigma2 * rdf, since rdf is a Satterthwaite df.
    expect_equal(fit$result$rss, fit$sigma2 * fit$sketch$kappa, label = method)
    if (identical(method, "ihs")) {
      expect_equal(fit$sketch$kappa, fit$result$rdf, label = method)
    }
  }
  # Direct check against the residuals of the sketched regression.
  Tlen <- 90L
  X <- cbind(1, rnorm(Tlen), cos(seq_len(Tlen) / 5))
  Z <- matrix(rnorm(Tlen * 5), Tlen)
  for (method in c("gaussian", "countsketch", "srht")) {
    set.seed(13)
    op <- fmrireg:::.lowrank_sketch_operator(Tlen, list(method = method, m = 30L))
    sol <- fmrireg:::.lowrank_time_solve(X, Z, list(method = method, m = 30L), op)
    Xs <- op$apply(X); Zs <- op$apply(Z)
    rs <- Zs - Xs %*% qr.solve(Xs, Zs)
    s2 <- fmrireg:::.lowrank_sigma2(sol, NULL, TRUE)
    expect_equal(s2 * sol$kappa, colSums(rs^2), tolerance = 1e-10, label = method)
  }
})

test_that("sketch-and-solve standard errors match the estimator's sampling SD", {
  # The sketch-and-solve estimate uses only m sketched rows, so its SD is about
  # sqrt(T / m) times the OLS SD. Reporting OLS-sized standard errors gave a
  # null type-I rate of 0.42 (SRHT) and 0.45 (Gaussian) at m = 40, T = 200.
  dset <- sketch_matrix_dataset(run_length = 60L, nvox = 3000L, seed = 21L)
  for (method in c("srht", "gaussian", "countsketch")) {
    set.seed(22)
    fit <- .sketch_fit(dset, method, m = 40L)
    b <- .beta_payload(fit)
    cols <- 1:2  # condition regressors, true effect 0
    ratio <- colMeans(b$se[, cols]) / sqrt(colMeans(b$est[, cols]^2))
    expect_true(all(abs(ratio - 1) < 0.1), label = paste(method, toString(round(ratio, 3))))
    type1 <- mean(b$prob[, cols] < 0.05)
    expect_lt(abs(type1 - 0.05), 0.02, label = paste(method, type1))
    # Residual df follow the sketch, not the full series.
    expect_lt(fit$result$rdf, 40)
    expect_gt(fit$result$rdf, 20)
  }
})

test_that("conditional sketch covariance is the sandwich for a non-orthogonal sketch", {
  # Gaussian and CountSketch rows are not orthonormal (S S' != c I), so the
  # covariance is (Xs'Xs)^-1 Xs' S S' Xs (Xs'Xs)^-1 and sigma^2 is estimated
  # from ||r_s||^2 / tr(P S S'). Check both against a Monte Carlo over noise
  # with the sketch held fixed.
  set.seed(31)
  Tlen <- 120L
  X <- cbind(1, rnorm(Tlen), sin(seq_len(Tlen) / 6))
  E <- matrix(rnorm(Tlen * 4000), Tlen)
  for (method in c("gaussian", "countsketch", "srht")) {
    op <- fmrireg:::.lowrank_sketch_operator(Tlen, list(method = method, m = 30L))
    sol <- fmrireg:::.lowrank_sketch_solve(X, E, op)
    emp <- apply(sol$M, 1, var)
    expect_true(all(abs(emp / diag(sol$cov_unscaled) - 1) < 0.1), label = method)
    s2 <- colSums(sol$residuals^2) / sol$kappa
    expect_lt(abs(mean(s2) - 1), 0.03, label = method)
    # Satterthwaite df: var(s2) = 2 / df.
    expect_lt(abs(var(s2) * sol$df / 2 - 1), 0.15, label = method)
  }
})

test_that("IHS fits reproduce OLS coefficients, standard errors and df", {
  dset <- sketch_matrix_dataset(run_length = 60L, nvox = 100L, signal = 0.5,
                                seed = 41L)
  set.seed(42)
  fit <- expect_no_warning(fmri_lm(
    onsets ~ hrf(condition), block = ~run, dataset = dset,
    engine = "latent_sketch",
    lowrank = lowrank_control(time_sketch = list(method = "ihs"))
  ))
  # Joint (all-runs) OLS on the same design.
  X <- as.matrix(design_matrix(fit$model))
  Y <- fmridataset::collect_assay(dset)
  XtXinv <- solve(crossprod(X))
  B_ols <- XtXinv %*% crossprod(X, Y)
  df <- nrow(X) - ncol(X)
  s2 <- colSums((Y - X %*% B_ols)^2) / df
  SE_ols <- sqrt(diag(XtXinv)) %o% sqrt(s2)

  b <- .beta_payload(fit)
  expect_lt(max(abs(t(b$est) - B_ols) / SE_ols), 5e-3)
  expect_equal(unname(t(b$se)), unname(SE_ols), tolerance = 1e-4)
  expect_equal(fit$sigma2, s2, tolerance = 1e-6)
  expect_equal(fit$result$rdf, df)
  # cov.unscaled carries the same rank attributes as exact fits' XtXinv.
  cu <- fit$result$cov.unscaled
  expect_true(attr(cu, "is_full_rank"))
  expect_length(attr(cu, "aliased"), 0L)
  expect_equal(unname(matrix(cu, nrow(cu))), unname(XtXinv), tolerance = 1e-8)
  expect_true(all(fit$sketch$converged))
})

test_that("by_cluster fits solve each cluster with its own whitened design", {
  # Each cluster has its own whitened design. The former code summed the
  # sketched Gram matrices of all clusters and divided each cluster's
  # cross-products by that sum, shrinking coefficients by about 1 / n_clusters
  # (and ran plain SRHT when "ihs" was requested).
  nvox <- 120L
  dset <- sketch_matrix_dataset(run_length = 60L, nvox = nvox, signal = 2,
                                seed = 51L)
  parcels <- rep(1:4, length.out = nvox)
  ctl <- fmri_lm_control(noise = noise_spec("ar1", pooling = "parcel",
                                            parcels = parcels))
  fit_by <- function(method) {
    set.seed(52)
    fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset,
            engine = "latent_sketch", control = ctl,
            lowrank = lowrank_control(parcels = parcels,
                                      time_sketch = list(method = method, m = 60L)))
  }
  fit_ihs <- fit_by("ihs")

  # Reference: per-cluster OLS on the cluster's whitened data, with the
  # cluster AR coefficients the fit used.
  X <- as.matrix(design_matrix(fit_ihs$model))
  Y <- fmridataset::collect_assay(dset)
  run_idx <- fmrireg:::.model_run_indices(fit_ihs$model, nrow(X))
  exact_first <- isTRUE(attr(fit_ihs, "config")$ar$exact_first)
  B_ref <- SE_ref <- matrix(NA_real_, ncol(X), nvox)
  for (g in as.character(1:4)) {
    J <- which(parcels == as.integer(g))
    w <- fmrireg:::ar_whiten_transform(X, Y[, J, drop = FALSE],
                                       fit_ihs$ar_coef[[g]], exact_first = exact_first,
                                       run_indices = run_idx)
    B_ref[, J] <- qr.solve(w$X, w$Y)
    s2 <- colSums((w$Y - w$X %*% B_ref[, J])^2) / (nrow(X) - ncol(X))
    SE_ref[, J] <- sqrt(diag(solve(crossprod(w$X)))) %o% sqrt(s2)
  }
  b <- .beta_payload(fit_ihs)
  expect_lt(max(abs(t(b$est) - B_ref) / SE_ref), 5e-3)
  expect_equal(unname(t(b$se)), SE_ref, tolerance = 1e-4)

  # Sketch-and-solve by cluster is unbiased for the per-cluster solution.
  fit_srht <- fit_by("srht")
  task <- 1:2
  slope <- coef(lm(as.numeric(fit_srht$betas_fixed[task, ]) ~
                     as.numeric(B_ref[task, ])))[[2]]
  expect_gt(slope, 0.85)
  expect_lt(slope, 1.15)
  expect_lt(abs(mean(fit_srht$sigma2) - 1), 0.15)

  # Cluster-specific covariance: post-hoc contrasts must be declared up front.
  expect_error(
    fit_contrasts(fit_srht, list(task = structure(1, colind = 1L))),
    "cluster-specific"
  )
})

test_that("by_cluster statistics stay aligned for unsorted labels and reject NA", {
  nvox <- 60L
  dset <- sketch_matrix_dataset(run_length = 60L, nvox = nvox, signal = 2,
                                seed = 61L)
  set.seed(62)
  parcels <- sample(c(7L, 2L, 5L), nvox, replace = TRUE)
  con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B",
                                    name = "A_vs_B"))
  ctl <- fmri_lm_control(noise = noise_spec("ar1", pooling = "parcel",
                                            parcels = parcels))
  set.seed(63)
  fit <- fmri_lm(onsets ~ hrf(condition, contrasts = con), block = ~run,
                 dataset = dset, engine = "latent_sketch", control = ctl,
                 lowrank = lowrank_control(parcels = parcels,
                                           time_sketch = list(method = "srht", m = 60L)))
  b <- .beta_payload(fit)
  expect_equal(unname(t(b$est)), unname(fit$betas_fixed))
  # Contrast rows are in voxel order: estimate = l'b and se from each voxel's
  # own cluster covariance.
  ct <- fit$result$contrasts
  hit <- grep("A_vs_B", ct$name, fixed = TRUE)
  expect_length(hit, 1L)
  cd <- ct$data[[hit]]
  l <- c(1, -1, rep(0, nrow(fit$betas_fixed) - 2L))
  expect_equal(cd$estimate, drop(l %*% fit$betas_fixed), tolerance = 1e-10)
  clusters <- fit$result$cluster_voxels
  se_expect <- numeric(nvox)
  for (i in seq_along(clusters)) {
    J <- clusters[[i]]
    V <- fit$result$covariance_by_cluster[[i]]
    se_expect[J] <- sqrt(drop(crossprod(l, V %*% l)) * fit$sigma2[J])
  }
  expect_equal(cd$se, se_expect, tolerance = 1e-10)

  parcels_na <- parcels; parcels_na[3] <- NA
  expect_error(
    fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset,
            engine = "latent_sketch",
            control = fmri_lm_control(noise = noise_spec(
              "ar1", pooling = "parcel", parcels = parcels_na)),
            lowrank = lowrank_control(parcels = parcels_na,
                                      time_sketch = list(method = "srht", m = 60L))),
    "must not contain NA"
  )
})

test_that("sketch fits report Satterthwaite df and cluster covariance scope", {
  nvox <- 60L
  dset <- sketch_matrix_dataset(run_length = 60L, nvox = nvox, seed = 71L)
  parcels <- rep(1:3, length.out = nvox)
  meta <- function(fit) fmrireg:::.fmri_lm_bids_inference_metadata(fit)

  set.seed(72)
  fit_g <- .sketch_fit(dset, "srht", m = 40L)
  expect_identical(fit_g$result$df$method, "satterthwaite")
  expect_identical(meta(fit_g)$DegreesOfFreedomMethod, "satterthwaite")
  expect_identical(meta(fit_g)$CovarianceScope, "shared")

  set.seed(72)
  fit_ihs <- .sketch_fit(dset, "ihs", m = 40L)
  expect_identical(meta(fit_ihs)$DegreesOfFreedomMethod, "residual")

  exact <- fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset)
  expect_identical(meta(exact)$DegreesOfFreedomMethod, "residual")

  ctl <- fmri_lm_control(noise = noise_spec("ar1", pooling = "parcel",
                                            parcels = parcels))
  set.seed(73)
  fit_c <- fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset,
                   engine = "latent_sketch", control = ctl,
                   lowrank = lowrank_control(parcels = parcels,
                                             time_sketch = list(method = "gaussian", m = 60L)))
  vm <- variance_model(fit_c)
  expect_identical(vm$covariance_scope, "cluster")
  expect_identical(vm$covariance, fit_c$result$covariance_by_cluster)
  expect_identical(meta(fit_c)$CovarianceScope, "cluster")
  expect_identical(meta(fit_c)$DegreesOfFreedomMethod, "satterthwaite")
  # One df for the map: the smallest cluster Satterthwaite df.
  expect_equal(meta(fit_c)$InferenceDegreesOfFreedom$Min, fit_c$result$rdf)
  expect_equal(meta(fit_c)$InferenceDegreesOfFreedom$Max, fit_c$result$rdf)
})

test_that("by_cluster HRF band uses the voxel's cluster covariance", {
  # by_cluster fits have no shared cov.unscaled; the band used to vanish
  # silently (and the curve was always drawn as significant).
  nvox <- 40L
  dset <- sketch_matrix_dataset(run_length = 60L, nvox = nvox, signal = 1,
                                seed = 81L)
  parcels <- rep(1:2, length.out = nvox)
  ctl <- fmri_lm_control(noise = noise_spec("ar1", pooling = "parcel",
                                            parcels = parcels))
  set.seed(82)
  fit <- fmri_lm(onsets ~ hrf(condition, basis = "spmg3"), block = ~run,
                 dataset = dset, engine = "latent_sketch", control = ctl,
                 lowrank = lowrank_control(parcels = parcels,
                                           time_sketch = list(method = "srht", m = 60L)))
  expect_null(fit$result$cov.unscaled)
  sa <- seq(0, 20, by = 0.5)
  G <- as.matrix(fmrihrf::HRF_SPMG3(sa))
  for (v in c(3L, 4L)) {  # one voxel in each cluster
    p <- ggplot2::autoplot(fit, type = "hrf", voxel = v, sample_at = sa)
    d <- p$data[p$data$condition == levels(p$data$condition)[1], ]
    expect_true(all(is.finite(d$lo)) && all(is.finite(d$hi)))
    k <- which(vapply(fit$result$cluster_voxels, function(g) v %in% g, logical(1)))
    V <- fit$result$covariance_by_cluster[[k]][1:3, 1:3] * fit$sigma2[v]
    half <- stats::qt(0.975, fit$result$rdf) * sqrt(rowSums((G %*% V) * G))
    expect_equal(drop(G %*% fit$betas_fixed[1:3, v]), d$mean, tolerance = 1e-8)
    expect_equal((d$hi - d$lo) / 2, half, tolerance = 1e-8)
  }
})

test_that("sketch-and-solve requires more sketch rows than design columns", {
  dset <- sketch_matrix_dataset()
  expect_error(.sketch_fit(dset, "srht", m = 8L), "must exceed the number of design columns")
  expect_error(lowrank_control(time_sketch = list(method = "ihs", iters = 0L)),
               "must be >= 1")
  expect_error(lowrank_control(time_sketch = list(method = "qr")),
               "must be one of")
})
