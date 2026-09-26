#' Internal: run low-rank/sketched engine under fmri_lm
#' @keywords internal
#' @noRd
.preflight_lowrank_engine <- function(fm, dataset, lowrank, cfg) {
  if (!inherits(fm, "fmri_model")) {
    stop("latent_sketch engine requires an 'fmri_model' object", call. = FALSE)
  }
  if (is.null(dataset)) {
    stop("latent_sketch engine requires a dataset", call. = FALSE)
  }
  if (!inherits(cfg, "fmri_lm_control")) {
    stop("latent_sketch engine requires an 'fmri_lm_control' object", call. = FALSE)
  }
  
  invisible(TRUE)
}

#' Internal: plugin-compatible preflight for low-rank engine
#' @keywords internal
#' @noRd
.preflight_lowrank_engine_plugin <- function(model, dataset, args, cfg) {
  .preflight_lowrank_engine(
    fm = model,
    dataset = dataset,
    lowrank = args$lowrank %||% list(),
    cfg = cfg
  )
}

#' Internal: map latent-space results back into voxel space
#' @keywords internal
#' @noRd
.lowrank_project_voxels <- function(x, A, A_is_I) {
  if (A_is_I) {
    return(x)
  }

  as.matrix(x %*% A)
}

#' Internal: blend fixed-order AR estimates, including run-specific lists
#' @keywords internal
#' @noRd
.lowrank_blend_ar <- function(local_phi, global_phi, alpha) {
  if (!is.list(local_phi) && !is.list(global_phi)) {
    return(alpha * local_phi + (1 - alpha) * global_phi)
  }

  if (!is.list(local_phi)) local_phi <- rep(list(local_phi), length(global_phi))
  if (!is.list(global_phi)) global_phi <- rep(list(global_phi), length(local_phi))
  if (length(local_phi) != length(global_phi)) {
    stop("run-specific local and global AR estimates have different lengths",
         call. = FALSE)
  }
  out <- Map(
    function(local, global) alpha * local + (1 - alpha) * global,
    local_phi, global_phi
  )
  names(out) <- names(local_phi) %||% names(global_phi)
  out
}

#' Internal: estimate low-rank AR once from initial OLS and whiten
#' @keywords internal
#' @noRd
.lowrank_whiten_initial_ols <- function(X, Z, ar_order, ar_opts,
                                       run_indices = NULL, censor = NULL) {
  proj <- .fast_preproject(X)
  residuals <- Z - X %*% (proj$Pinv %*% Z)
  phi <- .estimate_shared_ar_parameters(
    residuals = residuals,
    ar_order = ar_order,
    ar_opts = ar_opts,
    run_indices = run_indices,
    censor = censor,
    design = X
  )
  whitened <- ar_whiten_transform(
    X = X,
    Y = Z,
    phi = phi,
    exact_first = isTRUE(ar_opts$exact_first),
    censor = censor,
    run_indices = run_indices
  )
  list(phi = phi, X = whitened$X, Y = whitened$Y,
       residuals = residuals)
}

#' Internal: estimate and shrink parcel AR from OLS residual columns
#' @keywords internal
#' @noRd
.lowrank_group_ar_estimates <- function(residuals, sizes, ar_order, ar_opts,
                                        shrink_c0, design,
                                        run_indices = NULL, censor = NULL) {
  if (!is.list(residuals) || !length(residuals)) {
    stop("'residuals' must be a non-empty list of parcel residuals",
         call. = FALSE)
  }
  group_names <- names(residuals)
  if (is.null(group_names)) group_names <- as.character(seq_along(residuals))
  residual_matrix <- do.call(cbind, lapply(residuals, as.numeric))
  colnames(residual_matrix) <- group_names

  global_phi <- .estimate_shared_ar_parameters(
    residual_matrix, ar_order, ar_opts,
    run_indices = run_indices, censor = censor, design = design
  )
  phi_groups <- setNames(vector("list", length(residuals)), group_names)
  for (i in seq_along(residuals)) {
    local_phi <- .estimate_shared_ar_parameters(
      matrix(residuals[[i]], ncol = 1L), ar_order, ar_opts,
      run_indices = run_indices, censor = censor, design = design
    )
    group_size <- as.numeric(sizes[[group_names[[i]]]] %||% sizes[[i]])
    alpha <- group_size / (group_size + shrink_c0)
    phi_groups[[i]] <- .lowrank_blend_ar(local_phi, global_phi, alpha)
  }

  list(phi_groups = phi_groups, global_phi = global_phi)
}

#' Internal: validate and complete the time-sketch specification
#' @keywords internal
#' @noRd
.lowrank_resolve_time_sketch <- function(sk, p, Tlen) {
  sk <- sk %||% list()
  if (!is.list(sk)) stop("`time_sketch` must be a list", call. = FALSE)
  # Exact-name extraction: `sk$m` would partially match `sk$method`.
  sk <- list(method = sk[["method"]] %||% "gaussian", m = sk[["m"]],
             iters = sk[["iters"]], tol = sk[["tol"]])
  method <- sk$method
  if (!is.character(method) || length(method) != 1L ||
      !method %in% c("gaussian", "countsketch", "srht", "ihs")) {
    stop("`time_sketch$method` must be one of \"gaussian\", \"countsketch\", ",
         "\"srht\" or \"ihs\"", call. = FALSE)
  }
  sk$m <- as.integer(sk$m %||% min(8L * p, Tlen))
  if (length(sk$m) != 1L || is.na(sk$m) || sk$m < 1L || sk$m > Tlen) {
    stop(sprintf("`time_sketch$m` must be an integer in [1, %d]", Tlen),
         call. = FALSE)
  }
  if (sk$m <= p && !identical(method, "ihs")) {
    stop(sprintf(paste0(
      "`time_sketch$m` (%d) must exceed the number of design columns (%d): ",
      "sketch-and-solve needs residual degrees of freedom"), sk$m, p),
      call. = FALSE)
  }
  if (identical(method, "ihs")) {
    if (!is.null(sk$iters) && (length(sk$iters) != 1L || is.na(sk$iters) ||
                               sk$iters < 1L)) {
      stop("`time_sketch$iters` must be >= 1 for method = \"ihs\"", call. = FALSE)
    }
    if (!is.null(sk$tol) && (length(sk$tol) != 1L || !is.finite(sk$tol) ||
                             sk$tol < 0)) {
      stop("`time_sketch$tol` must be a non-negative number", call. = FALSE)
    }
  }
  sk
}

#' Internal: select landmark voxels and interpolation weights
#' @keywords internal
#' @noRd
.lowrank_landmarks <- function(dataset, lowrank) {
  L <- as.integer(lowrank$landmarks)
  mask <- .fmri_dataset_mask_space(dataset, "landmark selection")$mask
  coords <- neuroim2::index_to_coord(mask, which(as.vector(mask)))
  km_iter <- as.integer(lowrank$kmeans_iter_max %||% 1000L)
  km_nstart <- as.integer(lowrank$kmeans_nstart %||% 10L)
  km <- stats::kmeans(
    coords, centers = L, iter.max = km_iter, nstart = km_nstart,
    algorithm = "Lloyd"
  )
  idx <- as.integer(RANN::nn2(coords, km$centers, k = 1)$nn.idx[, 1])
  W <- build_landmark_weights(coords, coords[idx, , drop = FALSE],
                              k = as.integer(lowrank$k_neighbors %||% 16L))
  list(idx = idx, W = W)
}

#' Internal: residual variance from a time-sketched solve
#'
#' `sol$residuals` are sketched (sketch-and-solve) or full-data (IHS)
#' residuals of the latent/voxel columns; `sol$kappa` is their expected sum
#' of squares per unit noise variance, so the result is on the data scale.
#' @keywords internal
#' @noRd
.lowrank_sigma2 <- function(sol, A, A_is_I) {
  R <- sol$residuals
  if (A_is_I) return(colSums(R * R) / sol$kappa)
  # Voxel residual sums of squares diag(A' R'R A) via an r x r intermediate
  # rather than the dense (rows x V) projected residual matrix.
  RtR <- crossprod(R)
  as.numeric(Matrix::colSums(as.matrix(RtR %*% A) * A)) / sol$kappa
}

#' Internal: beta and contrast statistics for cluster-specific covariances
#'
#' Each cluster has its own whitened design and therefore its own coefficient
#' covariance; statistics are computed per cluster with the shared packagers
#' and stitched back into voxel order.
#' @keywords internal
#' @noRd
.lowrank_grouped_stats <- function(B, sigma2, groups, cov_list, dfres,
                                   varnames, ar_order, contrast_prep) {
  idx_all <- unlist(groups, use.names = FALSE)
  if (!identical(sort(idx_all), seq_len(ncol(B)))) {
    stop("internal error: cluster voxel sets must partition the voxels",
         call. = FALSE)
  }
  ord <- order(idx_all)
  sigma <- sqrt(pmax(sigma2, 0))

  pieces <- lapply(seq_along(groups), function(i) {
    J <- groups[[i]]
    beta_stats_matrix(
      Betas = B[, J, drop = FALSE], XtXinv = cov_list[[i]], sigma = sigma[J],
      dfres = dfres, varnames = varnames, ar_order = ar_order
    )
  })
  bstats <- pieces[[1L]]
  bdata <- bstats$data[[1L]]
  for (col in intersect(c("estimate", "se", "stat", "prob"), names(bdata))) {
    mats <- lapply(pieces, function(pc) pc$data[[1L]][[col]][[1L]])
    bdata[[col]] <- list(do.call(rbind, mats)[ord, , drop = FALSE])
  }
  if ("sigma" %in% names(bdata)) bdata$sigma <- list(sigma)
  bstats$data <- list(bdata)

  contrasts <- if (length(contrast_prep$standard) > 0L) {
    cpieces <- lapply(seq_along(groups), function(i) {
      J <- groups[[i]]
      fit_lm_contrasts_fast(
        B = B[, J, drop = FALSE], sigma2 = sigma2[J], XtXinv = cov_list[[i]],
        conlist = lapply(contrast_prep$simple, `[[`, "weights"),
        fconlist = lapply(contrast_prep$f, `[[`, "weights"),
        df = dfres, ar_order = ar_order
      )
    })
    stitched <- lapply(names(cpieces[[1L]]), function(nm) {
      out <- cpieces[[1L]][[nm]]
      data <- dplyr::bind_rows(lapply(cpieces, function(cp) cp[[nm]]$data[[1L]]))
      out$data <- list(data[ord, , drop = FALSE])
      out
    })
    dplyr::bind_rows(stitched)
  } else {
    empty_contrast_table()
  }
  list(bstats = bstats, contrasts = contrasts)
}

#' Internal: run low-rank/sketched engine under fmri_lm
#'
#' Inference follows the estimator that is returned. Sketch-and-solve
#' methods ("gaussian", "countsketch", "srht") report the conditional-on-sketch
#' covariance and an unbiased sketched residual variance with Satterthwaite
#' degrees of freedom (see `.lowrank_sketch_solve()`); "ihs" converges to the
#' full-data least-squares solution and reports exact OLS quantities.
#' @keywords internal
#' @noRd
.run_lowrank_engine <- function(fm, dataset, lowrank, cfg = NULL, ar_options = NULL) {
  if (is.null(cfg)) {
    cfg <- .fmri_lm_control_legacy(ar_options = ar_options)
  }
  lowrank <- lowrank %||% list()
  .preflight_lowrank_engine(fm, dataset, lowrank, cfg)

  ar_opts <- cfg$ar %||% list()
  # Design (T x p)
  X <- as.matrix(design_matrix(fm))
  Tlen <- nrow(X); p <- ncol(X)
  run_indices <- .model_run_indices(fm, Tlen)
  censor <- resolve_censor(cfg, dataset = dataset, n_time = Tlen)
  varnames <- colnames(X)

  # Latent basis and loadings or full data path
  if (.dset_is_latent(dataset)) {
    Z <- as.matrix(.dset_data_matrix(dataset))       # T x r component scores
    lds <- .dset_loadings(dataset)                   # V x r (dgCMatrix or matrix)
    if (inherits(lds, "Matrix")) {
      A <- Matrix::t(lds)                            # r x V (sparse)
    } else {
      A <- t(lds)
    }
    A_is_I <- FALSE
  } else {
    # Treat Z as full voxel data (T x V) and A as identity (V x V)
    Z <- as.matrix(.dset_data_matrix(dataset))
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
  exact_first <- isTRUE(ar_opts$exact_first)
  shrink_c0 <- as.integer(ar_opts$shrink_c0 %||% 100L)
  no_whiten <- ar_order <= 0L
  ar_coef_store <- NULL

  # One sketch (shared across all branches and clusters)
  sk <- .lowrank_resolve_time_sketch(lowrank$time_sketch, p, Tlen)
  op <- .lowrank_sketch_operator(Tlen, sk)

  grouped <- NULL
  if (no_whiten) {
    # The AR paths below warn about a rank-deficient X through
    # .fast_preproject(); this path has no preliminary OLS, so warn here.
    sol <- .lowrank_time_solve(X, Z, sk, op, warn = TRUE)
    B <- .lowrank_project_voxels(sol$M, A, A_is_I)
    sigma2 <- .lowrank_sigma2(sol, A, A_is_I)
    cov_unscaled <- sol$cov_unscaled
    dfres <- sol$df
    solve_info <- list(sol)
  } else if (by_cluster && !.dset_is_latent(dataset) && !is.null(lowrank$parcels)) {
    # --- Grouped (by parcel) whitening path, full-voxel dataset only ---
    gids <- if (inherits(lowrank$parcels, "ClusteredNeuroVol")) {
      as.integer(neuroim2::values(lowrank$parcels))
    } else {
      as.integer(lowrank$parcels)
    }
    if (length(gids) != ncol(Z)) stop("parcels/group ids length must equal number of voxels")
    if (anyNA(gids)) {
      stop("parcels/group ids must not contain NA: every voxel needs a cluster ",
           "for by_cluster AR whitening", call. = FALSE)
    }
    ug <- sort(unique(gids))
    groups <- lapply(ug, function(g) which(gids == g))
    names(groups) <- as.character(ug)

    # Estimate AR per group from parcel-mean OLS residuals, shrunk to global
    Pinv <- .fast_preproject(X)$Pinv
    res_per_group <- lapply(groups, function(Jg) {
      ybar_g <- rowMeans(Z[, Jg, drop = FALSE])
      ybar_g - drop(X %*% (Pinv %*% ybar_g))
    })
    sizes <- vapply(groups, length, integer(1))
    group_ar <- .lowrank_group_ar_estimates(
      residuals = res_per_group,
      sizes = sizes,
      ar_order = ar_order,
      ar_opts = ar_opts,
      shrink_c0 = shrink_c0,
      design = X,
      run_indices = run_indices,
      censor = censor
    )
    phi_global <- group_ar$global_phi
    phi_groups <- group_ar$phi_groups

    # Whiten and solve each cluster with its own design. Every cluster has a
    # different whitened design, hence its own normal equations and its own
    # coefficient covariance; the sketch itself is shared.
    M <- matrix(0, p, ncol(Z))
    sigma2 <- numeric(ncol(Z))
    rss_cluster <- numeric(ncol(Z))
    # Estimability is judged per cluster on that cluster's whitened design;
    # each covariance carries its own rank attributes.
    cov_list <- vector("list", length(groups))
    solve_info <- vector("list", length(groups))
    for (i in seq_along(groups)) {
      Jg <- groups[[i]]
      tmp <- ar_whiten_transform(
        X, Z[, Jg, drop = FALSE], phi_groups[[names(groups)[i]]],
        exact_first = exact_first, censor = censor,
        run_indices = run_indices
      )
      sol_g <- .lowrank_time_solve(tmp$X, tmp$Y, sk, op)
      M[, Jg] <- sol_g$M
      sigma2[Jg] <- .lowrank_sigma2(sol_g, NULL, TRUE)
      rss_cluster[Jg] <- sigma2[Jg] * sol_g$kappa
      cov_list[[i]] <- sol_g$cov_unscaled
      solve_info[[i]] <- sol_g
    }
    B <- M
    # One residual df for the whole map: the smallest across clusters.
    dfres <- min(vapply(solve_info, `[[`, numeric(1), "df"))
    cov_unscaled <- NULL
    grouped <- list(groups = groups, cov_list = cov_list)
    attr(phi_groups, "global_phi") <- phi_global
    ar_coef_store <- phi_groups
  } else {
    # --- Global AR path ---
    # Estimate the configured shared covariance from the full OLS residual
    # matrix. pooled_acvf targets a typical column; mean_series explicitly
    # targets the coherent spatial component.
    initial_ar <- .lowrank_whiten_initial_ols(
      X = X,
      Z = Z,
      ar_order = ar_order,
      ar_opts = ar_opts,
      run_indices = run_indices,
      censor = censor
    )
    phi <- initial_ar$phi
    Xw <- initial_ar$X
    Zw <- initial_ar$Y
    if (A_is_I && !is.null(lowrank$landmarks)) {
      # Landmark solve + Nystrom extension
      lm <- .lowrank_landmarks(dataset, lowrank)
      sol <- .lowrank_time_solve(Xw, Zw[, lm$idx, drop = FALSE], sk, op)
      B <- extend_betas_landmarks(sol$M, lm$W)
      # Propagate landmark residual variances through squared weights
      sigma2_L <- .lowrank_sigma2(sol, NULL, TRUE)
      W2 <- lm$W; W2@x <- W2@x * W2@x
      sigma2 <- as.numeric(W2 %*% sigma2_L)
    } else {
      sol <- .lowrank_time_solve(Xw, Zw, sk, op)
      B <- .lowrank_project_voxels(sol$M, A, A_is_I)
      sigma2 <- .lowrank_sigma2(sol, A, A_is_I)
    }
    cov_unscaled <- sol$cov_unscaled
    dfres <- sol$df
    solve_info <- list(sol)
    ar_coef_store <- if (is.list(phi)) phi else list(phi)
  }

  # Build fmri_lm-like result structure compatible with downstream code
  sigma <- sqrt(pmax(sigma2, 0))
  contrast_prep <- prepare_fmri_lm_contrasts(fm)
  if (is.null(grouped)) {
    bstats <- beta_stats_matrix(
      Betas = B,
      XtXinv = cov_unscaled,
      sigma = sigma,
      dfres = dfres,
      varnames = varnames,
      ar_order = ar_order
    )
    contrast_results <- if (length(contrast_prep$standard) > 0L) {
      dplyr::bind_rows(
        fit_lm_contrasts_fast(
          B = B,
          sigma2 = sigma2,
          XtXinv = cov_unscaled,
          conlist = lapply(contrast_prep$simple, `[[`, "weights"),
          fconlist = lapply(contrast_prep$f, `[[`, "weights"),
          df = dfres,
          ar_order = ar_order
        )
      )
    } else {
      empty_contrast_table()
    }
  } else {
    gs <- .lowrank_grouped_stats(
      B = B, sigma2 = sigma2, groups = grouped$groups,
      cov_list = grouped$cov_list, dfres = dfres, varnames = varnames,
      ar_order = ar_order, contrast_prep = contrast_prep
    )
    bstats <- gs$bstats
    contrast_results <- gs$contrasts
  }

  # Event/baseline indices for coef() methods
  tmats <- term_matrices(fm)
  event_indices <- attr(tmats, "event_term_indices")
  baseline_indices <- attr(tmats, "baseline_term_indices")
  # `rss` is the residual sum of squares of the rows actually fitted: the
  # sketched residuals ||r_s||^2 for sketch-and-solve, the full-data RSS for
  # IHS. Its expectation is sigma^2 * kappa with kappa = tr(PK) (the residual
  # count T - p for IHS), so resvar = rss / kappa; kappa is not the
  # Satterthwaite rdf. sigma2 = ||r_s||^2 / kappa exactly, so this recovers
  # the fitted RSS without re-forming the residuals (for landmark fits it is
  # the landmark RSS interpolated like sigma2).
  kappa <- vapply(solve_info, function(s) as.numeric(s$kappa), numeric(1))
  rss <- if (is.null(grouped)) sigma2 * kappa else rss_cluster

  # Keep-but-aliased: the reported coefficients of aliased columns are NA
  # (the stats above already are); B itself keeps zeros there so that
  # estimable contrasts are not poisoned by 0 * NA.
  betas_report <- B
  if (is.null(grouped)) {
    al <- attr(cov_unscaled, "aliased", exact = TRUE)
    if (length(al)) betas_report[al, ] <- NA_real_
  } else {
    for (i in seq_along(grouped$groups)) {
      al <- attr(grouped$cov_list[[i]], "aliased", exact = TRUE)
      if (length(al)) betas_report[al, grouped$groups[[i]]] <- NA_real_
    }
  }

  sketch_info <- list(
    method = sk$method,
    m = sk$m,
    df = dfres,
    kappa = kappa,
    inference = if (identical(sk$method, "ihs")) "ols" else "sketch_conditional",
    iters = vapply(solve_info, function(s) as.integer(s$iters), integer(1)),
    converged = vapply(solve_info, function(s) as.logical(s$converged), logical(1))
  )

  if (identical(sk$method, "ihs") && any(!sketch_info$converged, na.rm = TRUE)) {
    warning(sprintf(paste0(
      "IHS did not reach tol = %g OLS standard errors within %d iterations ",
      "in %d of %d solve(s); coefficients are not at the least-squares ",
      "solution. Increase `time_sketch$iters` or `time_sketch$m`."),
      sk$tol %||% .lowrank_ihs_default_tol,
      as.integer(sk$iters %||% .lowrank_ihs_default_iters),
      sum(!sketch_info$converged, na.rm = TRUE), length(sketch_info$converged)),
      call. = FALSE)
  }

  result <- list(
    betas = bstats,
    contrasts = contrast_results,
    event_indices = event_indices,
    baseline_indices = baseline_indices,
    cov.unscaled = cov_unscaled,
    sigma = sigma,
    rdf = dfres,
    rss = rss,
    resvar = sigma2,
    ar_coef = ar_coef_store,
    sketch = sketch_info,
    # Sketch-and-solve rdf is a Satterthwaite df (see .lowrank_sketch_solve());
    # IHS reports exact OLS df.
    df_method = if (identical(sk$method, "ihs")) "residual" else "satterthwaite"
  )
  if (!is.null(grouped)) {
    result$covariance_by_cluster <- grouped$cov_list
    result$cluster_voxels <- grouped$groups
    result$contrast_scope <- list(
      allowed_colind = integer(0),
      mode = "error",
      reason = paste0(
        "Parcel-pooled (by_cluster) sketch fits have a cluster-specific ",
        "coefficient covariance; specify contrasts in the model so they are ",
        "computed at fit time."
      )
    )
  }

  ret <- list(
    result = result,
    model = fm,
    strategy = "sketch",
    bcons = contrast_prep$processed,
    dataset = dataset,
    betas_fixed = betas_report,
    sigma2 = sigma2,
    vcov_inv = cov_unscaled,
    ar_coef = ar_coef_store,
    sketch = sketch_info
  )
  class(ret) <- "fmri_lm"
  attr(ret, "strategy") <- "sketch"
  attr(ret, "config") <- cfg
  ret
}


#' Internal: plugin-compatible fit for low-rank engine
#' @keywords internal
#' @noRd
.fit_lowrank_engine_plugin <- function(model, dataset, args, cfg) {
  .run_lowrank_engine(
    fm = model,
    dataset = dataset,
    lowrank = args$lowrank %||% list(),
    cfg = cfg
  )
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
  if (engine == "sketch") {
    engine <- "latent_sketch"
  }

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
  .run_lowrank_engine(fm, dataset, lowrank, cfg = cfg, ar_options = ar_options)
}
