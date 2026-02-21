#' Internal: run low-rank/sketched engine under fmri_lm
#' @keywords internal
#' @noRd
.run_lowrank_engine <- function(fm, dataset, lowrank, cfg = NULL, ar_options = NULL) {
  if (is.null(cfg)) {
    cfg <- fmri_lm_control(ar_options = ar_options)
  }

  ar_opts <- cfg$ar %||% list()
  robust_opts <- cfg$robust %||% list()
  # Design (T x p)
  X <- as.matrix(design_matrix(fm))
  Tlen <- nrow(X); p <- ncol(X)
  # Residual df should reflect original time samples, not sketch rows.
  dfres_full <- max(1L, Tlen - p)
  varnames <- colnames(X)

  # Latent basis and loadings or full data path
  if (inherits(dataset, "latent_dataset")) {
    Z <- as.matrix(fmridataset::get_latent_scores(dataset))   # T x r
    # Extract LatentNeuroVec from backend
    lvec <- if (!is.null(dataset$lvec)) {
      dataset$lvec
    } else if (!is.null(dataset$backend) && !is.null(dataset$backend$data)) {
      dataset$backend$data[[1]]
    } else {
      stop("Cannot find LatentNeuroVec in latent_dataset")
    }
    lds <- lvec@loadings                             # V x r (dgCMatrix or matrix)
    if (inherits(lds, "Matrix")) {
      A <- Matrix::t(lds)                            # r x V (sparse)
    } else {
      A <- t(lds)
    }
    A_is_I <- FALSE
  } else {
    # Fallback: treat Z as full voxel data (T x V) and A as identity (V x V)
    Zfull <- as.matrix(fmridataset::get_data_matrix(dataset))
    Z <- Zfull
    A <- NULL
    A_is_I <- TRUE
  }

  # --- AR prewhitening options ---
  by_cluster <- isTRUE(ar_opts$by_cluster)
  ar_struct <- ar_opts$struct %||% "iid"
  ar_order <- switch(as.character(ar_struct),
                     "iid" = 0L,
                     "ar1" = 1L,
                     "ar2" = 2L,
                     "ar3" = 3L,
                     "ar4" = 4L,
                     "arp" = as.integer(ar_opts$p %||% 0L),
                     0L)
  ar_order <- if (is.finite(ar_order)) ar_order else 0L
  iter_gls <- as.integer(ar_opts$iter_gls %||% 1L)
  exact_first <- isTRUE(ar_opts$exact_first)
  shrink_c0 <- as.integer(ar_opts$shrink_c0 %||% 100L)
  no_whiten <- ar_order <= 0L
  ar_coef_store <- NULL

  # Build sketch (shared across all branches)
  sk <- lowrank$time_sketch %||% list(method = "gaussian", m = min(8L * p, Tlen))
  if (is.null(sk$m)) sk$m <- min(8L * p, Tlen)
  S <- make_time_sketch(Tlen, sk)

  # --- No-whitening short-circuit (diagnostics) ---
  if (no_whiten || ar_order <= 0L) {
    # Sketch and solve without temporal whitening
    if (identical(sk$method, "ihs")) {
      sol <- ihs_latent_solve(X, Z, m = sk$m, iters = as.integer(sk$iters %||% 3L))
      M <- sol$M; Ginv <- sol$Ginv
      B <- if (A_is_I) M else M %*% as.matrix(A)
      # Residuals via one SRHT draw
      plan <- make_srht_plan(Tlen, sk$m)
      Xs <- srht_apply(X, plan); Zs <- srht_apply(Z, plan)
      Rres <- Zs - Xs %*% M
      RA <- if (A_is_I) Rres else Rres %*% as.matrix(A)
      dfres <- dfres_full
      sigma2 <- colSums(RA * RA) / dfres
    } else if (identical(sk$method, "srht")) {
      plan <- make_srht_plan(Tlen, sk$m)
      Xs <- srht_apply(X, plan); Zs <- srht_apply(Z, plan)
      G <- crossprod(Xs)
      Ginv <- tryCatch(chol2inv(chol(G)), error = function(e) {
        ridge <- 1e-6 * sum(diag(G)) / max(1L, ncol(G))
        chol2inv(chol(G + diag(ridge, ncol(G))))
      })
      R <- crossprod(Xs, Zs)
      M <- Ginv %*% R
      B <- if (A_is_I) M else M %*% as.matrix(A)
      Rres <- Zs - Xs %*% M
      RA <- if (A_is_I) Rres else Rres %*% as.matrix(A)
      dfres <- dfres_full
      sigma2 <- colSums(RA * RA) / dfres
    } else if (identical(sk$method, "gaussian")) {
      Sg <- S
      Xs <- Sg %*% X; Zs <- Sg %*% Z
      G <- crossprod(Xs)
      Ginv <- tryCatch(chol2inv(chol(G)), error = function(e) {
        ridge <- 1e-6 * sum(diag(G)) / max(1L, ncol(G))
        chol2inv(chol(G + diag(ridge, ncol(G))))
      })
      R <- crossprod(Xs, Zs)
      M <- Ginv %*% R
      B <- if (A_is_I) M else M %*% as.matrix(A)
      Rres <- Zs - Xs %*% M
      RA <- if (A_is_I) Rres else Rres %*% as.matrix(A)
      dfres <- dfres_full
      sigma2 <- colSums(RA * RA) / dfres
    } else { # countsketch
      Xs <- as.matrix(S %*% X); Zs <- as.matrix(S %*% Z)
      G <- crossprod(Xs)
      Ginv <- tryCatch(chol2inv(chol(G)), error = function(e) {
        ridge <- 1e-6 * sum(diag(G)) / max(1L, ncol(G))
        chol2inv(chol(G + diag(ridge, ncol(G))))
      })
      R <- crossprod(Xs, Zs)
      M <- Ginv %*% R
      B <- if (A_is_I) M else M %*% as.matrix(A)
      Rres <- Zs - Xs %*% M
      RA <- if (A_is_I) Rres else Rres %*% as.matrix(A)
      dfres <- dfres_full
      sigma2 <- colSums(RA * RA) / dfres
    }
  } else if (by_cluster && !inherits(dataset, "latent_dataset") && !is.null(lowrank$parcels)) {
    # --- Grouped (by parcel) whitening path, full-voxel dataset only ---
    # Extract group ids per voxel
    gids <- if (inherits(lowrank$parcels, "ClusteredNeuroVol")) {
      as.integer(neuroim2::values(lowrank$parcels))
    } else {
      as.integer(lowrank$parcels)
    }
    if (length(gids) != ncol(Z)) stop("parcels/group ids length must equal number of voxels")
    ug <- sort(unique(gids))

    # Accumulators
    Gsum <- matrix(0, p, p)
    Rall <- matrix(0, p, ncol(Z))
    rss <- numeric(ncol(Z))

    # Estimate AR per group from parcel-mean residuals
    # Global precompute for OLS pinch
    XtX <- crossprod(X)
    Pinv <- tryCatch(chol2inv(chol(XtX)) %*% t(X), error = function(e) MASS::ginv(XtX) %*% t(X))
    # First pass: collect residuals and sizes for shrinkage
    res_per_group <- vector("list", length(ug))
    sizes <- integer(length(ug))
    names(res_per_group) <- as.character(ug)
    names(sizes) <- as.character(ug)
    for (g in ug) {
      Jg <- which(gids == g)
      if (length(Jg) == 0) next
      ybar_g <- rowMeans(Z[, Jg, drop = FALSE])
      beta_g <- Pinv %*% ybar_g
      resid_g <- ybar_g - drop(X %*% beta_g)
      res_per_group[[as.character(g)]] <- resid_g
      sizes[as.character(g)] <- length(Jg)
    }
    # Global phi for shrinkage
    phi_global <- .estimate_ar_parameters_routed(unlist(res_per_group, use.names = FALSE), ar_order)
    # Second pass: whiten, sketch, accumulate
    plan <- if (identical(sk$method, "srht") || identical(sk$method, "ihs")) make_srht_plan(Tlen, sk$m) else NULL
    phi_groups <- vector("list", length(ug))
    names(phi_groups) <- as.character(ug)
    for (g in ug) {
      Jg <- which(gids == g)
      if (length(Jg) == 0) next
      resid_g <- res_per_group[[as.character(g)]]
      alpha_g <- sizes[as.character(g)]/(sizes[as.character(g)] + shrink_c0)
      phi_g_raw <- .estimate_ar_parameters_routed(resid_g, ar_order)
      phi_g <- alpha_g * phi_g_raw + (1 - alpha_g) * phi_global
      phi_groups[[as.character(g)]] <- phi_g
      # Whiten group
      tmp <- ar_whiten_transform(X, Z[, Jg, drop = FALSE], phi_g, exact_first = exact_first)
      Xw_g <- tmp$X; Zw_g <- tmp$Y
      # Sketch
      if (identical(sk$method, "srht") || identical(sk$method, "ihs")) {
        Xs_g <- srht_apply(Xw_g, plan); Zs_g <- srht_apply(Zw_g, plan)
      } else if (identical(sk$method, "gaussian")) {
        Xs_g <- S %*% Xw_g; Zs_g <- S %*% Zw_g
      } else { # countsketch
        Xs_g <- as.matrix(S %*% Xw_g); Zs_g <- as.matrix(S %*% Zw_g)
      }
      # Accumulate cross-products
      Gsum <- Gsum + crossprod(Xs_g)
      Rall[, Jg] <- crossprod(Xs_g, Zs_g)
    }

    # Solve and betas
    Ginv <- tryCatch(chol2inv(chol(Gsum)), error = function(e) {
      ridge <- 1e-6 * sum(diag(Gsum)) / max(1L, ncol(Gsum))
      chol2inv(chol(Gsum + diag(ridge, ncol(Gsum))))
    })
    M <- Ginv %*% Rall
    B <- if (A_is_I) M else M %*% as.matrix(A)

    # Residuals per group for sigma2
    dfres <- dfres_full
    rss <- rss * 0
    for (g in ug) {
      Jg <- which(gids == g)
      if (length(Jg) == 0) next
      # Recompute whiten+sketch for residuals using the pooled phi estimates
      phi_g <- phi_groups[[as.character(g)]] %||% phi_global
      if (is.null(phi_g) || !length(phi_g)) {
        ybar_g <- rowMeans(Z[, Jg, drop = FALSE])
        beta_g <- Pinv %*% ybar_g
        resid_g <- ybar_g - drop(X %*% beta_g)
        phi_g <- .estimate_ar_parameters_routed(resid_g, ar_order)
      }
      tmp <- ar_whiten_transform(X, Z[, Jg, drop = FALSE], phi_g, exact_first = exact_first)
      Xw_g <- tmp$X; Zw_g <- tmp$Y
      if (identical(sk$method, "srht") || identical(sk$method, "ihs")) {
        Xs_g <- srht_apply(Xw_g, plan)
        Zs_g <- srht_apply(Zw_g, plan)
      } else if (identical(sk$method, "gaussian")) {
        Xs_g <- S %*% Xw_g
        Zs_g <- S %*% Zw_g
      } else { # countsketch
        Xs_g <- as.matrix(S %*% Xw_g)
        Zs_g <- as.matrix(S %*% Zw_g)
      }
      Eg <- Zs_g - Xs_g %*% M[, Jg, drop = FALSE]
      RA_g <- if (A_is_I) Eg else Eg %*% as.matrix(A)
      rss[Jg] <- colSums(RA_g * RA_g)
    }
    sigma2 <- rss / dfres
    attr(phi_groups, "global_phi") <- phi_global
    ar_coef_store <- phi_groups

  } else {
    # --- Global AR path ---
    # Estimate global AR from mean residuals over columns
    ybar <- rowMeans(Z)
    XtX <- crossprod(X)
    Pinv <- tryCatch(chol2inv(chol(XtX)) %*% t(X), error = function(e) MASS::ginv(XtX) %*% t(X))
    beta <- Pinv %*% ybar
    resid <- ybar - drop(X %*% beta)
    phi <- .estimate_ar_parameters_routed(resid, ar_order)

    if (iter_gls > 1L) {
      phi_current <- phi
      for (it in seq_len(iter_gls)) {
        tmp_iter <- ar_whiten_transform(X, Z, phi_current, exact_first = exact_first)
        Xw <- tmp_iter$X
        Zw <- tmp_iter$Y
        if (it < iter_gls) {
          XtX_iter <- crossprod(Xw)
          XtXinv_iter <- tryCatch(chol2inv(chol(XtX_iter)),
                                  error = function(e) MASS::ginv(XtX_iter))
          beta_iter <- XtXinv_iter %*% crossprod(Xw, Zw)
          resid_iter <- rowMeans(Z - X %*% beta_iter)
          phi_current <- .estimate_ar_parameters_routed(resid_iter, ar_order)
        }
      }
      phi <- phi_current
    }

    if (!exists("Xw", inherits = FALSE) || !exists("Zw", inherits = FALSE)) {
      tmp <- ar_whiten_transform(X, Z, phi, exact_first = exact_first)
      Xw <- tmp$X
      Zw <- tmp$Y
    }
    # Sketch and solve
    if (identical(sk$method, "ihs")) {
      if (A_is_I && !is.null(lowrank$landmarks)) {
        # Landmark solve + Nyström extension
        L <- as.integer(lowrank$landmarks)
        mask <- fmridataset::get_mask(dataset)
        coords <- neuroim2::index_to_coord(mask, which(as.vector(mask)))
        km_iter <- as.integer(lowrank$kmeans_iter_max %||% 1000L)
        km_nstart <- as.integer(lowrank$kmeans_nstart %||% 10L)
        km <- stats::kmeans(coords, centers = L, iter.max = km_iter, nstart = km_nstart)
        idx_lm <- as.integer(RANN::nn2(coords, km$centers, k = 1)$nn.idx[, 1])
        lcoords <- coords[idx_lm, , drop = FALSE]
        # Solve only on landmarks
        Zw_L <- Zw[, idx_lm, drop = FALSE]
        sol <- ihs_latent_solve(Xw, Zw_L, m = sk$m, iters = as.integer(sk$iters %||% 3L))
        M_L <- sol$M; Ginv <- sol$Ginv
        BL <- M_L  # p x L
        # Weights and extension
        W <- build_landmark_weights(coords, lcoords, k = as.integer(lowrank$k_neighbors %||% 16L))
        B <- extend_betas_landmarks(BL, W)
        # Variance propagate from landmark residuals via one SRHT draw
        plan <- make_srht_plan(Tlen, sk$m)
        Xs <- srht_apply(Xw, plan); Zs_L <- srht_apply(Zw_L, plan)
        Rres_L <- Zs_L - Xs %*% M_L
        dfres <- dfres_full
        sigma2_L <- colSums(Rres_L * Rres_L) / dfres
        W2 <- W; W2@x <- W2@x * W2@x
        sigma2 <- as.numeric(W2 %*% sigma2_L)
      } else {
        sol <- ihs_latent_solve(Xw, Zw, m = sk$m, iters = as.integer(sk$iters %||% 3L))
        M <- sol$M; Ginv <- sol$Ginv
        B <- if (A_is_I) M else M %*% as.matrix(A)
        # Residuals for sigma2 via one SRHT draw
        plan <- make_srht_plan(Tlen, sk$m)
        Xs <- srht_apply(Xw, plan); Zs <- srht_apply(Zw, plan)
        Rres <- Zs - Xs %*% M
        RA <- if (A_is_I) Rres else Rres %*% as.matrix(A)
        dfres <- dfres_full
        sigma2 <- colSums(RA * RA) / dfres
      }
    } else if (identical(sk$method, "srht")) {
      plan <- make_srht_plan(Tlen, sk$m)
      Xs <- srht_apply(Xw, plan)
      if (A_is_I && !is.null(lowrank$landmarks)) {
        L <- as.integer(lowrank$landmarks)
        mask <- fmridataset::get_mask(dataset)
        coords <- neuroim2::index_to_coord(mask, which(as.vector(mask)))
        km_iter <- as.integer(lowrank$kmeans_iter_max %||% 1000L)
        km_nstart <- as.integer(lowrank$kmeans_nstart %||% 10L)
        km <- stats::kmeans(coords, centers = L, iter.max = km_iter, nstart = km_nstart)
        idx_lm <- as.integer(RANN::nn2(coords, km$centers, k = 1)$nn.idx[, 1])
        lcoords <- coords[idx_lm, , drop = FALSE]
        Zs_L <- srht_apply(Zw[, idx_lm, drop = FALSE], plan)
        G <- crossprod(Xs)
        Ginv <- tryCatch(chol2inv(chol(G)), error = function(e) {
          ridge <- 1e-6 * sum(diag(G)) / max(1L, ncol(G))
          chol2inv(chol(G + diag(ridge, ncol(G))))
        })
        R_L <- crossprod(Xs, Zs_L)
        M_L <- Ginv %*% R_L
        BL <- M_L
        W <- build_landmark_weights(coords, lcoords, k = as.integer(lowrank$k_neighbors %||% 16L))
        B <- extend_betas_landmarks(BL, W)
        Rres_L <- Zs_L - Xs %*% M_L
        dfres <- dfres_full
        sigma2_L <- colSums(Rres_L * Rres_L) / dfres
        W2 <- W; W2@x <- W2@x * W2@x
        sigma2 <- as.numeric(W2 %*% sigma2_L)
      } else {
        Zs <- srht_apply(Zw, plan)
        G <- crossprod(Xs)
        Ginv <- tryCatch(chol2inv(chol(G)), error = function(e) {
          ridge <- 1e-6 * sum(diag(G)) / max(1L, ncol(G))
          chol2inv(chol(G + diag(ridge, ncol(G))))
        })
        R <- crossprod(Xs, Zs)
        M <- Ginv %*% R
        B <- if (A_is_I) M else M %*% as.matrix(A)
        Rres <- Zs - Xs %*% M
        RA <- if (A_is_I) Rres else Rres %*% as.matrix(A)
        dfres <- dfres_full
        sigma2 <- colSums(RA * RA) / dfres
      }
    } else if (identical(sk$method, "gaussian")) {
      Xs <- S %*% Xw
      if (A_is_I && !is.null(lowrank$landmarks)) {
        L <- as.integer(lowrank$landmarks)
        mask <- fmridataset::get_mask(dataset)
        coords <- neuroim2::index_to_coord(mask, which(as.vector(mask)))
        km_iter <- as.integer(lowrank$kmeans_iter_max %||% 1000L)
        km_nstart <- as.integer(lowrank$kmeans_nstart %||% 10L)
        km <- stats::kmeans(coords, centers = L, iter.max = km_iter, nstart = km_nstart)
        idx_lm <- as.integer(RANN::nn2(coords, km$centers, k = 1)$nn.idx[, 1])
        lcoords <- coords[idx_lm, , drop = FALSE]
        Zs_L <- S %*% Zw[, idx_lm, drop = FALSE]
        G <- crossprod(Xs)
        Ginv <- tryCatch(chol2inv(chol(G)), error = function(e) {
          ridge <- 1e-6 * sum(diag(G)) / max(1L, ncol(G))
          chol2inv(chol(G + diag(ridge, ncol(G))))
        })
        R_L <- crossprod(Xs, Zs_L)
        M_L <- Ginv %*% R_L
        BL <- M_L
        W <- build_landmark_weights(coords, lcoords, k = as.integer(lowrank$k_neighbors %||% 16L))
        B <- extend_betas_landmarks(BL, W)
        Rres_L <- Zs_L - Xs %*% M_L
        dfres <- dfres_full
        sigma2_L <- colSums(Rres_L * Rres_L) / dfres
        W2 <- W; W2@x <- W2@x * W2@x
        sigma2 <- as.numeric(W2 %*% sigma2_L)
      } else {
        Zs <- S %*% Zw
        G <- crossprod(Xs)
        Ginv <- tryCatch(chol2inv(chol(G)), error = function(e) {
          ridge <- 1e-6 * sum(diag(G)) / max(1L, ncol(G))
          chol2inv(chol(G + diag(ridge, ncol(G))))
        })
        R <- crossprod(Xs, Zs)
        M <- Ginv %*% R
        B <- if (A_is_I) M else M %*% as.matrix(A)
        Rres <- Zs - Xs %*% M
        RA <- if (A_is_I) Rres else Rres %*% as.matrix(A)
        dfres <- dfres_full
        sigma2 <- colSums(RA * RA) / dfres
      }
    } else { # countsketch
      Xs <- as.matrix(S %*% Xw)
      if (A_is_I && !is.null(lowrank$landmarks)) {
        L <- as.integer(lowrank$landmarks)
        mask <- fmridataset::get_mask(dataset)
        coords <- neuroim2::index_to_coord(mask, which(as.vector(mask)))
        km_iter <- as.integer(lowrank$kmeans_iter_max %||% 1000L)
        km_nstart <- as.integer(lowrank$kmeans_nstart %||% 10L)
        km <- stats::kmeans(coords, centers = L, iter.max = km_iter, nstart = km_nstart)
        idx_lm <- as.integer(RANN::nn2(coords, km$centers, k = 1)$nn.idx[, 1])
        lcoords <- coords[idx_lm, , drop = FALSE]
        Zs_L <- as.matrix(S %*% Zw[, idx_lm, drop = FALSE])
        G <- crossprod(Xs)
        Ginv <- tryCatch(chol2inv(chol(G)), error = function(e) {
          ridge <- 1e-6 * sum(diag(G)) / max(1L, ncol(G))
          chol2inv(chol(G + diag(ridge, ncol(G))))
        })
        R_L <- crossprod(Xs, Zs_L)
        M_L <- Ginv %*% R_L
        BL <- M_L
        W <- build_landmark_weights(coords, lcoords, k = as.integer(lowrank$k_neighbors %||% 16L))
        B <- extend_betas_landmarks(BL, W)
        Rres_L <- Zs_L - Xs %*% M_L
        dfres <- dfres_full
        sigma2_L <- colSums(Rres_L * Rres_L) / dfres
        W2 <- W; W2@x <- W2@x * W2@x
        sigma2 <- as.numeric(W2 %*% sigma2_L)
      } else {
        Zs <- as.matrix(S %*% Zw)
        G <- crossprod(Xs)
        Ginv <- tryCatch(chol2inv(chol(G)), error = function(e) {
          ridge <- 1e-6 * sum(diag(G)) / max(1L, ncol(G))
          chol2inv(chol(G + diag(ridge, ncol(G))))
        })
        R <- crossprod(Xs, Zs)
        M <- Ginv %*% R
        B <- if (A_is_I) M else M %*% as.matrix(A)
        Rres <- Zs - Xs %*% M
        RA <- if (A_is_I) Rres else Rres %*% as.matrix(A)
        dfres <- dfres_full
        sigma2 <- colSums(RA * RA) / dfres
      }
    }
    ar_coef_store <- list(phi)
  }

  # Build fmri_lm-like result structure compatible with downstream code
  # Use matrix-based beta_stats packager for consistent tibble output
  sigma <- sqrt(pmax(sigma2, 0))
  bstats <- beta_stats_matrix(Betas = B, XtXinv = Ginv, sigma = sigma,
                              dfres = dfres, varnames = varnames)

  # Event/baseline indices for coef() methods
  tmats <- term_matrices(fm)
  event_indices <- attr(tmats, "event_term_indices")
  baseline_indices <- attr(tmats, "baseline_term_indices")

  result <- list(
    betas = bstats,
    contrasts = dplyr::tibble(),
    event_indices = event_indices,
    baseline_indices = baseline_indices,
    sigma = sigma,
    rdf = dfres,
    ar_coef = ar_coef_store
  )

  ret <- list(
    result = result,
    model = fm,
    strategy = "sketch",
    bcons = list(),
    dataset = dataset,
    betas_fixed = B,
    sigma2 = sigma2,
    vcov_inv = Ginv,
    ar_coef = ar_coef_store
  )
  class(ret) <- "fmri_lm"
  attr(ret, "strategy") <- "sketch"
  attr(ret, "config") <- cfg
  ret
}

#' Internal: dispatch fmri_lm low-rank/sketch engine
#' @keywords internal
#' @noRd
fmri_lm_lowrank_dispatch <- function(formula_or_model, dataset, engine = NULL, lowrank = NULL,
                                     block = NULL, baseline_model = NULL,
                                     durations = 0, drop_empty = TRUE,
                                     cfg = NULL, ar_options = NULL) {
  if (is.null(engine)) return(NULL)
  engine <- match.arg(engine, c("latent_sketch", "sketch"))
  if (engine == "latent_sketch") {
    # Build fmri_model reusing the standard path
    fm <- if (inherits(formula_or_model, "fmri_model")) {
      formula_or_model
    } else {
      create_fmri_model(formula_or_model,
                        block = block,
                        baseline_model = baseline_model,
                        dataset = dataset,
                        drop_empty = drop_empty,
                        durations = durations)
    }
    return(.run_lowrank_engine(fm, dataset, lowrank, cfg = cfg, ar_options = ar_options))
  }
  NULL
}
