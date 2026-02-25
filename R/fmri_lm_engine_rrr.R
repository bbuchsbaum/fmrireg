#' Internal: preflight checks for reduced-rank GLS engine
#' @keywords internal
#' @noRd
.preflight_rrr_gls_engine <- function(model, dataset, args, cfg) {
  if (!inherits(model, "fmri_model")) {
    stop("rrr_gls engine requires an 'fmri_model' object", call. = FALSE)
  }
  if (is.null(dataset)) {
    stop("rrr_gls engine requires a dataset", call. = FALSE)
  }
  if (!inherits(cfg, "fmri_lm_config")) {
    stop("rrr_gls engine requires an 'fmri_lm_config' object", call. = FALSE)
  }

  args <- .rrr_normalize_args(args)

  tmats <- term_matrices(model)
  event_indices <- attr(tmats, "event_term_indices")
  if (is.null(event_indices) || length(event_indices) == 0L) {
    stop("rrr_gls engine requires at least one event/task regressor", call. = FALSE)
  }

  robust_type <- cfg$robust$type %||% FALSE
  if (!identical(robust_type, FALSE)) {
    stop("rrr_gls does not support robust fitting; set robust = FALSE", call. = FALSE)
  }

  if (isTRUE(cfg$ar$voxelwise)) {
    stop("rrr_gls supports only shared (non-voxelwise) temporal covariance", call. = FALSE)
  }

  if (isTRUE(cfg$ar$by_cluster)) {
    stop("rrr_gls does not support parcel-specific AR whitening", call. = FALSE)
  }

  invisible(TRUE)
}


#' Internal: reduced-rank GLS engine fit
#' @keywords internal
#' @noRd
.fit_rrr_gls_engine <- function(model, dataset, args, cfg) {
  .preflight_rrr_gls_engine(model, dataset, args, cfg)
  args <- .rrr_normalize_args(args)

  dataset <- dataset %||% model$dataset
  Y <- .rrr_extract_response_matrix(dataset)
  mode(Y) <- "numeric"

  X_full <- as.matrix(design_matrix(model))
  mode(X_full) <- "numeric"
  if (nrow(X_full) != nrow(Y)) {
    stop("Row mismatch between design matrix and data matrix", call. = FALSE)
  }

  tmats <- term_matrices(model)
  event_indices <- as.integer(attr(tmats, "event_term_indices"))
  baseline_indices <- as.integer(attr(tmats, "baseline_term_indices") %||% integer(0))
  varnames <- colnames(X_full) %||% attr(tmats, "varnames")
  if (is.null(varnames) || length(varnames) != ncol(X_full)) {
    varnames <- paste0("beta_", seq_len(ncol(X_full)))
  }

  X_task <- X_full[, event_indices, drop = FALSE]
  Z_nuisance <- if (length(baseline_indices) > 0L) {
    X_full[, baseline_indices, drop = FALSE]
  } else {
    NULL
  }

  run_indices <- .rrr_derive_run_indices(model, nrow(X_full))
  white <- .rrr_whiten_data(
    X_task = X_task,
    Z_nuisance = Z_nuisance,
    Y = Y,
    cfg = cfg,
    run_indices = run_indices
  )

  resid <- .rrr_residualize_against_nuisance(
    Xw_task = white$Xw_task,
    Yw = white$Yw,
    Zw = white$Zw
  )
  X0 <- resid$X0
  Y0 <- resid$Y0
  rank_z <- resid$rank_z

  rank_control <- list(
    rank_mode = args$rank_mode,
    rank = args$rank,
    energy_keep = args$energy_keep,
    rss_budget = args$rss_budget
  )
  task_fit <- .rrr_fit_task_subspace(X0, Y0, rank_control)
  B_task <- task_fit$B_task
  C_task <- task_fit$C_task
  V_r <- task_fit$V_r
  D_task <- task_fit$XtXinv_task

  rank_x <- qr(X0, LAPACK = TRUE)$rank
  dfres <- max(1L, nrow(X0) - rank_z - rank_x)

  if (task_fit$rank_used > 0L && ncol(V_r) > 0L) {
    Y_r <- Y0 %*% V_r
    E_r <- Y_r - X0 %*% C_task
    Sigma_r <- crossprod(E_r) / dfres
    s2_cond <- rowSums((V_r %*% Sigma_r) * V_r)
    s2_cond <- pmax(as.numeric(s2_cond), .Machine$double.eps)
  } else {
    Sigma_r <- matrix(0, nrow = 0L, ncol = 0L)
    s2_cond <- rep(.Machine$double.eps, ncol(Y0))
  }

  se_event_cond <- outer(
    sqrt(s2_cond),
    sqrt(pmax(diag(D_task), 0))
  )

  contrast_prep <- prepare_fmri_lm_contrasts(model)
  simple_weights <- lapply(contrast_prep$simple, `[[`, "weights")
  if (length(simple_weights)) {
    names(simple_weights) <- names(contrast_prep$simple)
  }
  f_weights <- lapply(contrast_prep$f, `[[`, "weights")
  if (length(f_weights)) {
    names(f_weights) <- names(contrast_prep$f)
  }

  scoped_contrasts <- .rrr_filter_task_contrasts(
    simple_weights = simple_weights,
    f_weights = f_weights,
    event_indices = event_indices,
    policy = args$contrast_policy
  )

  boot <- NULL
  if (identical(args$se_mode, "bootstrap")) {
    boot <- .rrr_bootstrap_task_se(
      X0 = X0,
      Y0 = Y0,
      B_task_hat = B_task,
      rank_control = rank_control,
      t_contrasts_remapped = scoped_contrasts$simple,
      nboot = args$bootstrap_n,
      block_size = args$bootstrap_block_size %||% 1L,
      seed = args$bootstrap_seed
    )
  }

  se_event <- if (!is.null(boot)) t(boot$se_beta) else se_event_cond

  B_full <- matrix(0, nrow = ncol(X_full), ncol = ncol(Y0))
  B_full[event_indices, ] <- B_task

  if (!is.null(white$Zw) && ncol(white$Zw) > 0L && length(baseline_indices) > 0L) {
    proj_z <- .fast_preproject(white$Zw)
    gamma_hat <- proj_z$Pinv %*% (white$Yw - white$Xw_task %*% B_task)
    B_full[baseline_indices, ] <- gamma_hat
  }

  fitted_w <- white$Xw_task %*% B_task
  if (!is.null(white$Zw) && ncol(white$Zw) > 0L && length(baseline_indices) > 0L) {
    fitted_w <- fitted_w + white$Zw %*% B_full[baseline_indices, , drop = FALSE]
  }
  resid_w <- white$Yw - fitted_w
  rss <- colSums(resid_w^2)
  resvar <- pmax(rss / dfres, .Machine$double.eps)

  sigma <- sqrt(s2_cond)

  bstats <- .rrr_build_beta_stats(
    B_full = B_full,
    event_indices = event_indices,
    se_event = se_event,
    dfres = dfres,
    varnames = varnames,
    sigma = sigma
  )

  contrast_results <- if (length(scoped_contrasts$simple) > 0L || length(scoped_contrasts$f) > 0L) {
    fit_lm_contrasts_fast(
      B = B_task,
      sigma2 = s2_cond,
      XtXinv = D_task,
      conlist = scoped_contrasts$simple,
      fconlist = scoped_contrasts$f,
      df = dfres
    )
  } else {
    list()
  }

  contrast_results <- .rrr_restore_contrast_metadata(
    contrast_results = contrast_results,
    original_simple = scoped_contrasts$original_simple,
    original_f = scoped_contrasts$original_f
  )

  if (!is.null(boot) && length(boot$se_t) > 0L) {
    for (i in seq_along(contrast_results)) {
      cres <- contrast_results[[i]]
      nm <- cres$name[1]
      if (!identical(cres$type[1], "contrast") || is.null(boot$se_t[[nm]])) {
        next
      }
      data_tbl <- cres$data[[1]]
      est <- as.numeric(data_tbl$estimate)
      se_boot <- as.numeric(boot$se_t[[nm]])
      stat <- ifelse(abs(se_boot) < .Machine$double.eps^0.5, 0, est / se_boot)
      prob <- 2 * pt(-abs(stat), df = dfres)

      data_tbl$se <- se_boot
      data_tbl$stat <- stat
      data_tbl$prob <- prob
      data_tbl$sigma <- sigma
      contrast_results[[i]]$data[[1]] <- data_tbl
    }
  }

  combined_contrasts <- if (length(contrast_results) > 0L) {
    dplyr::bind_rows(contrast_results)
  } else {
    tibble::tibble()
  }

  D_full <- matrix(0, nrow = ncol(X_full), ncol = ncol(X_full))
  D_full[event_indices, event_indices] <- D_task
  dimnames(D_full) <- list(varnames, varnames)

  result <- list(
    contrasts = combined_contrasts,
    betas = bstats,
    event_indices = event_indices,
    baseline_indices = baseline_indices,
    cov.unscaled = D_full,
    sigma = sigma,
    rss = rss,
    rdf = dfres,
    resvar = resvar,
    ar_coef = if (length(white$phi) > 0L) list(white$phi) else NULL,
    contrast_scope = list(
      allowed_colind = event_indices,
      mode = "error",
      reason = "rrr_gls supports inference only for event/task parameters."
    )
  )

  ret <- list(
    result = result,
    model = model,
    strategy = "engine",
    bcons = contrast_prep$processed,
    dataset = dataset,
    ar_coef = result$ar_coef,
    betas_fixed = B_full,
    sigma2 = s2_cond,
    vcov_inv = D_full,
    rrr = list(
      rank_mode = args$rank_mode,
      rank_used = task_fit$rank_used,
      singular_values = task_fit$singular_values,
      energy_retained = task_fit$energy_retained,
      delta_rss = task_fit$delta_rss,
      se_mode = args$se_mode,
      dropped_contrasts = scoped_contrasts$dropped,
      bootstrap = if (!is.null(boot)) {
        list(nboot = args$bootstrap_n, block_size = args$bootstrap_block_size %||% 1L)
      } else {
        NULL
      }
    )
  )

  class(ret) <- "fmri_lm"
  attr(ret, "config") <- cfg
  attr(ret, "strategy") <- "engine"
  attr(ret, "engine") <- "rrr_gls"
  ret
}


#' Internal: normalize/validate reduced-rank engine arguments
#' @keywords internal
#' @noRd
.rrr_normalize_args <- function(args) {
  if (is.null(args)) {
    args <- list()
  } else if (!is.list(args)) {
    args <- list(args)
  }

  out <- utils::modifyList(
    list(
      rank = NULL,
      rank_mode = "fixed",
      energy_keep = 0.99,
      rss_budget = NULL,
      se_mode = "conditional",
      bootstrap_n = 200L,
      bootstrap_block_size = NULL,
      bootstrap_seed = NULL,
      contrast_policy = "warn_drop"
    ),
    args
  )

  if (!is.null(out$energy) && is.null(args$energy_keep)) {
    out$energy_keep <- out$energy
  }
  if (!is.null(out$nboot) && is.null(args$bootstrap_n)) {
    out$bootstrap_n <- out$nboot
  }

  out$rank_mode <- match.arg(as.character(out$rank_mode), c("fixed", "energy", "rss_budget"))
  out$se_mode <- match.arg(as.character(out$se_mode), c("conditional", "bootstrap"))
  out$contrast_policy <- match.arg(as.character(out$contrast_policy), c("warn_drop", "drop", "error"))

  if (!is.null(out$rank)) {
    out$rank <- as.integer(out$rank)
    if (length(out$rank) != 1L || !is.finite(out$rank) || out$rank < 1L) {
      stop("rrr_gls `rank` must be a positive integer", call. = FALSE)
    }
  }

  out$energy_keep <- as.numeric(out$energy_keep)[1]
  if (out$rank_mode == "energy" &&
      (!is.finite(out$energy_keep) || out$energy_keep <= 0 || out$energy_keep > 1)) {
    stop("rrr_gls `energy_keep` must be in (0, 1] when rank_mode = 'energy'", call. = FALSE)
  }

  if (!is.null(out$rss_budget)) {
    out$rss_budget <- as.numeric(out$rss_budget)[1]
  }
  if (out$rank_mode == "rss_budget" &&
      (is.null(out$rss_budget) || !is.finite(out$rss_budget) || out$rss_budget < 0)) {
    stop("rrr_gls `rss_budget` must be a non-negative scalar when rank_mode = 'rss_budget'", call. = FALSE)
  }

  out$bootstrap_n <- as.integer(out$bootstrap_n %||% 200L)
  if (out$se_mode == "bootstrap" && (!is.finite(out$bootstrap_n) || out$bootstrap_n < 2L)) {
    stop("rrr_gls `bootstrap_n` must be >= 2 when se_mode = 'bootstrap'", call. = FALSE)
  }

  if (!is.null(out$bootstrap_block_size)) {
    out$bootstrap_block_size <- as.integer(out$bootstrap_block_size)
    if (!is.finite(out$bootstrap_block_size) || out$bootstrap_block_size < 1L) {
      stop("rrr_gls `bootstrap_block_size` must be a positive integer", call. = FALSE)
    }
  }

  if (!is.null(out$bootstrap_seed)) {
    out$bootstrap_seed <- as.integer(out$bootstrap_seed)[1]
  }

  out
}


#' Internal: extract response matrix from dataset
#' @keywords internal
#' @noRd
.rrr_extract_response_matrix <- function(dataset) {
  if (inherits(dataset, "latent_dataset")) {
    Z <- as.matrix(fmridataset::get_latent_scores(dataset))
    lvec <- if (!is.null(dataset$lvec)) {
      dataset$lvec
    } else if (!is.null(dataset$backend) && !is.null(dataset$backend$data)) {
      dataset$backend$data[[1]]
    } else {
      stop("rrr_gls: cannot find latent loadings in latent_dataset", call. = FALSE)
    }
    L <- as.matrix(lvec@loadings)
    return(Z %*% t(L))
  }

  as.matrix(fmridataset::get_data_matrix(dataset))
}


#' Internal: derive run indices from model sampling frame
#' @keywords internal
#' @noRd
.rrr_derive_run_indices <- function(model, n) {
  sframe <- tryCatch(model$event_model$sampling_frame, error = function(e) NULL)
  if (is.null(sframe)) {
    return(NULL)
  }

  runs <- tryCatch(fmrihrf::blockids(sframe), error = function(e) NULL)
  if (!is.numeric(runs) || length(runs) != n) {
    return(NULL)
  }

  split_idx <- split(seq_len(n), runs)
  lapply(split_idx, as.integer)
}


#' Internal: map cfg$ar into AR order
#' @keywords internal
#' @noRd
.rrr_ar_order_from_cfg <- function(cfg) {
  ar_struct <- as.character(cfg$ar$struct %||% "iid")
  order <- switch(
    ar_struct,
    iid = 0L,
    ar1 = 1L,
    ar2 = 2L,
    ar3 = 3L,
    ar4 = 4L,
    arp = as.integer(cfg$ar$p %||% 0L),
    0L
  )
  if (!is.finite(order)) {
    0L
  } else {
    as.integer(order)
  }
}


#' Internal: apply global AR whitening to design and data
#' @keywords internal
#' @noRd
.rrr_whiten_data <- function(X_task, Z_nuisance, Y, cfg, run_indices = NULL) {
  X_ref <- if (!is.null(Z_nuisance) && ncol(Z_nuisance) > 0L) {
    cbind(X_task, Z_nuisance)
  } else {
    X_task
  }

  ar_order <- .rrr_ar_order_from_cfg(cfg)
  if (ar_order <= 0L) {
    return(list(
      Xw_task = X_task,
      Zw = Z_nuisance,
      Yw = Y,
      phi = numeric(0),
      ar_order = 0L
    ))
  }

  iter_gls <- as.integer(cfg$ar$iter_gls %||% 1L)
  if (!is.finite(iter_gls) || iter_gls < 1L) {
    iter_gls <- 1L
  }
  exact_first <- isTRUE(cfg$ar$exact_first)
  censor <- cfg$ar$censor %||% NULL

  estimate_phi <- function(Xm, Ym) {
    proj <- .fast_preproject(Xm)
    mean_y <- matrix(rowMeans(Ym), ncol = 1)
    beta <- proj$Pinv %*% mean_y
    resid <- rowMeans(Ym) - drop(Xm %*% beta)
    .estimate_ar_parameters_routed(
      residuals_vec = resid,
      ar_order = ar_order,
      run_indices = run_indices,
      censor = censor
    )
  }

  phi <- estimate_phi(X_ref, Y)
  Xw_ref <- NULL
  Yw <- NULL

  for (it in seq_len(iter_gls)) {
    tmp <- ar_whiten_transform(
      X = X_ref,
      Y = Y,
      phi = phi,
      exact_first = exact_first,
      censor = censor
    )
    Xw_ref <- as.matrix(tmp$X)
    Yw <- as.matrix(tmp$Y)

    if (it < iter_gls) {
      phi <- estimate_phi(Xw_ref, Yw)
    }
  }

  p_task <- ncol(X_task)
  Xw_task <- Xw_ref[, seq_len(p_task), drop = FALSE]
  Zw <- if (ncol(Xw_ref) > p_task) {
    Xw_ref[, seq.int(p_task + 1L, ncol(Xw_ref)), drop = FALSE]
  } else {
    NULL
  }

  list(
    Xw_task = Xw_task,
    Zw = Zw,
    Yw = Yw,
    phi = phi,
    ar_order = ar_order
  )
}


#' Internal: residualize task design and response against nuisance design
#' @keywords internal
#' @noRd
.rrr_residualize_against_nuisance <- function(Xw_task, Yw, Zw) {
  if (is.null(Zw) || ncol(Zw) == 0L) {
    return(list(X0 = Xw_task, Y0 = Yw, rank_z = 0L))
  }

  qr_z <- qr(Zw, LAPACK = TRUE)
  rank_z <- qr_z$rank
  if (rank_z <= 0L) {
    return(list(X0 = Xw_task, Y0 = Yw, rank_z = 0L))
  }

  Qz <- qr.Q(qr_z, complete = FALSE)[, seq_len(rank_z), drop = FALSE]
  Y0 <- Yw - Qz %*% crossprod(Qz, Yw)
  X0 <- Xw_task - Qz %*% crossprod(Qz, Xw_task)

  list(X0 = X0, Y0 = Y0, rank_z = rank_z)
}


#' Internal: choose truncation rank from singular values
#' @keywords internal
#' @noRd
.rrr_choose_rank <- function(d, rank_control) {
  d <- as.numeric(d)
  d <- d[is.finite(d) & d > 0]
  if (length(d) == 0L) {
    return(0L)
  }

  mode <- rank_control$rank_mode %||% "fixed"
  rank_user <- rank_control$rank
  energy_keep <- rank_control$energy_keep %||% 0.99
  rss_budget <- rank_control$rss_budget

  r <- switch(
    mode,
    fixed = {
      if (is.null(rank_user)) length(d) else as.integer(rank_user)
    },
    energy = {
      ss <- d^2
      total <- sum(ss)
      if (total <= 0) {
        length(d)
      } else {
        target <- as.numeric(energy_keep)[1]
        which(cumsum(ss) / total >= target)[1]
      }
    },
    rss_budget = {
      ss <- d^2
      deltas <- c(if (length(ss) > 1L) rev(cumsum(rev(ss)))[2:length(ss)] else numeric(0), 0)
      budget <- as.numeric(rss_budget %||% 0)[1]
      idx <- which(deltas <= budget)[1]
      if (is.na(idx)) length(ss) else idx
    },
    length(d)
  )

  r <- as.integer(r %||% length(d))
  r <- max(1L, min(r, length(d)))
  r
}


#' Internal: fit reduced-rank task model in residualized space
#' @keywords internal
#' @noRd
.rrr_fit_task_subspace <- function(X0, Y0, rank_control) {
  X0 <- as.matrix(X0)
  Y0 <- as.matrix(Y0)
  k <- ncol(X0)
  V <- ncol(Y0)

  qr_x <- qr(X0, LAPACK = TRUE)
  rank_x <- qr_x$rank
  if (rank_x <= 0L) {
    return(list(
      B_task = matrix(0, nrow = k, ncol = V),
      C_task = matrix(0, nrow = k, ncol = 0L),
      V_r = matrix(0, nrow = V, ncol = 0L),
      XtXinv_task = .fast_preproject(X0)$XtXinv,
      singular_values = numeric(0),
      rank_used = 0L,
      energy_retained = 0,
      delta_rss = 0
    ))
  }

  Qx_full <- qr.Q(qr_x, complete = FALSE)
  Qx <- Qx_full[, seq_len(rank_x), drop = FALSE]
  R_full <- qr.R(qr_x, complete = FALSE)
  R <- R_full[seq_len(rank_x), seq_len(rank_x), drop = FALSE]
  piv <- qr_x$pivot[seq_len(rank_x)]

  G <- crossprod(Qx, Y0)
  S <- tcrossprod(G)
  eig <- eigen(S, symmetric = TRUE)
  sing <- sqrt(pmax(eig$values, 0))

  tol <- max(dim(G)) * .Machine$double.eps * max(c(sing, 1))
  keep <- sing > tol
  if (!any(keep)) {
    return(list(
      B_task = matrix(0, nrow = k, ncol = V),
      C_task = matrix(0, nrow = k, ncol = 0L),
      V_r = matrix(0, nrow = V, ncol = 0L),
      XtXinv_task = .fast_preproject(X0)$XtXinv,
      singular_values = sing,
      rank_used = 0L,
      energy_retained = 0,
      delta_rss = 0
    ))
  }

  U <- eig$vectors[, keep, drop = FALSE]
  d <- sing[keep]
  rank_used <- .rrr_choose_rank(d, rank_control)

  U_r <- U[, seq_len(rank_used), drop = FALSE]
  d_r <- d[seq_len(rank_used)]

  V_r <- t(G) %*% U_r
  V_r <- sweep(V_r, 2L, d_r, `/`)

  UDr <- sweep(U_r, 2L, d_r, `*`)
  C_pivot <- backsolve(R, UDr)
  C_task <- matrix(0, nrow = k, ncol = rank_used)
  C_task[piv, ] <- C_pivot
  B_task <- C_task %*% t(V_r)

  ss <- d^2
  total_ss <- sum(ss)
  kept_ss <- sum(ss[seq_len(rank_used)])
  delta_rss <- if (rank_used < length(ss)) sum(ss[(rank_used + 1L):length(ss)]) else 0
  energy <- if (total_ss > 0) kept_ss / total_ss else 1

  list(
    B_task = B_task,
    C_task = C_task,
    V_r = V_r,
    XtXinv_task = .fast_preproject(X0)$XtXinv,
    singular_values = sing,
    rank_used = rank_used,
    energy_retained = energy,
    delta_rss = delta_rss
  )
}


#' Internal: scope contrasts to event/task regressors
#' @keywords internal
#' @noRd
.rrr_filter_task_contrasts <- function(simple_weights, f_weights, event_indices, policy = "warn_drop") {
  policy <- match.arg(policy, c("warn_drop", "drop", "error"))
  event_indices <- as.integer(event_indices)
  event_map <- setNames(seq_along(event_indices), as.character(event_indices))
  dropped <- character(0)

  remap <- function(weight_list) {
    if (length(weight_list) == 0L) {
      return(list(remapped = list(), original = list()))
    }

    remapped <- list()
    original <- list()

    for (nm in names(weight_list)) {
      w <- weight_list[[nm]]
      colind <- attr(w, "colind")
      if (is.null(colind)) {
        colind <- if (is.matrix(w)) seq_len(ncol(w)) else seq_len(length(w))
      }
      colind <- as.integer(colind)

      if (!all(colind %in% event_indices)) {
        dropped <<- c(dropped, nm)
        next
      }

      mapped <- as.integer(event_map[as.character(colind)])
      w2 <- w
      attr(w2, "colind") <- mapped

      remapped[[nm]] <- w2
      original[[nm]] <- w
    }

    list(remapped = remapped, original = original)
  }

  t_out <- remap(simple_weights)
  f_out <- remap(f_weights)

  if (length(dropped) > 0L) {
    dropped <- unique(dropped)
    msg <- paste0(
      "rrr_gls supports contrasts on event/task parameters only; dropped contrasts: ",
      paste(dropped, collapse = ", ")
    )
    if (policy == "error") {
      stop(msg, call. = FALSE)
    }
    if (policy == "warn_drop") {
      warning(msg, call. = FALSE)
    }
  }

  list(
    simple = t_out$remapped,
    f = f_out$remapped,
    original_simple = t_out$original,
    original_f = f_out$original,
    dropped = unique(dropped)
  )
}


#' Internal: restore original contrast metadata after remapping
#' @keywords internal
#' @noRd
.rrr_restore_contrast_metadata <- function(contrast_results, original_simple, original_f) {
  if (length(contrast_results) == 0L) {
    return(contrast_results)
  }

  for (i in seq_along(contrast_results)) {
    cres <- contrast_results[[i]]
    nm <- cres$name[1]
    orig <- original_simple[[nm]] %||% original_f[[nm]]
    if (is.null(orig)) {
      next
    }
    cres$conmat[[1]] <- orig
    orig_colind <- attr(orig, "colind")
    if (!is.null(orig_colind)) {
      cres$colind[[1]] <- as.integer(orig_colind)
    }
    contrast_results[[i]] <- cres
  }

  contrast_results
}


#' Internal: build beta statistics table with task-only inference
#' @keywords internal
#' @noRd
.rrr_build_beta_stats <- function(B_full, event_indices, se_event, dfres, varnames, sigma) {
  B_full <- as.matrix(B_full)
  V <- ncol(B_full)
  p <- nrow(B_full)

  estimate <- t(B_full)
  se <- matrix(NA_real_, nrow = V, ncol = p)
  se[, event_indices] <- se_event

  stat <- estimate / se
  prob <- 2 * pt(-abs(stat), df = dfres)

  colnames(estimate) <- varnames
  colnames(se) <- varnames
  colnames(stat) <- varnames
  colnames(prob) <- varnames

  dplyr::tibble(
    type = "beta",
    name = "parameter_estimates",
    stat_type = "tstat",
    df.residual = dfres,
    conmat = list(NULL),
    colind = list(NULL),
    data = list(dplyr::tibble(
      estimate = list(estimate),
      se = list(se),
      stat = list(stat),
      prob = list(prob),
      sigma = list(sigma)
    ))
  )
}


#' Internal: bootstrap task-parameter SEs with row/block resampling
#' @keywords internal
#' @noRd
.rrr_bootstrap_task_se <- function(X0, Y0, B_task_hat, rank_control,
                                   t_contrasts_remapped, nboot,
                                   block_size = 1L, seed = NULL) {
  X0 <- as.matrix(X0)
  Y0 <- as.matrix(Y0)
  B_task_hat <- as.matrix(B_task_hat)

  n <- nrow(X0)
  k <- ncol(X0)
  V <- ncol(Y0)
  nboot <- as.integer(nboot)
  block_size <- as.integer(block_size %||% 1L)
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  fitted <- X0 %*% B_task_hat
  resid <- Y0 - fitted

  mean_beta <- matrix(0, nrow = k, ncol = V)
  m2_beta <- matrix(0, nrow = k, ncol = V)

  t_vectors <- list()
  if (length(t_contrasts_remapped) > 0L) {
    for (nm in names(t_contrasts_remapped)) {
      w <- t_contrasts_remapped[[nm]]
      colind <- as.integer(attr(w, "colind"))
      l <- numeric(k)
      l[colind] <- as.numeric(w)
      t_vectors[[nm]] <- l
    }
  }
  mean_t <- lapply(t_vectors, function(x) numeric(V))
  m2_t <- lapply(t_vectors, function(x) numeric(V))

  blocks <- if (block_size > 1L) {
    starts <- seq.int(1L, n, by = block_size)
    lapply(starts, function(s) seq.int(s, min(n, s + block_size - 1L)))
  } else {
    NULL
  }

  sample_idx <- function() {
    if (is.null(blocks)) {
      return(sample.int(n, size = n, replace = TRUE))
    }
    idx <- integer(0)
    while (length(idx) < n) {
      idx <- c(idx, blocks[[sample.int(length(blocks), 1L)]])
    }
    idx[seq_len(n)]
  }

  for (b in seq_len(nboot)) {
    idx <- sample_idx()
    Yb <- fitted + resid[idx, , drop = FALSE]
    fit_b <- .rrr_fit_task_subspace(X0, Yb, rank_control)
    B_b <- fit_b$B_task

    delta_b <- B_b - mean_beta
    mean_beta <- mean_beta + delta_b / b
    m2_beta <- m2_beta + delta_b * (B_b - mean_beta)

    if (length(t_vectors) > 0L) {
      for (nm in names(t_vectors)) {
        est_b <- as.numeric(t_vectors[[nm]] %*% B_b)
        delta_t <- est_b - mean_t[[nm]]
        mean_t[[nm]] <- mean_t[[nm]] + delta_t / b
        m2_t[[nm]] <- m2_t[[nm]] + delta_t * (est_b - mean_t[[nm]])
      }
    }
  }

  se_beta <- sqrt(pmax(m2_beta / (nboot - 1L), 0))
  se_t <- lapply(m2_t, function(x) sqrt(pmax(x / (nboot - 1L), 0)))

  list(
    se_beta = se_beta,
    se_t = se_t
  )
}


#' Internal: register built-in engines
#' @keywords internal
#' @noRd
.register_builtin_engines <- function() {
  if (is.null(get_engine("rrr_gls"))) {
    register_engine(
      name = "rrr_gls",
      preflight = .preflight_rrr_gls_engine,
      fit = .fit_rrr_gls_engine,
      capabilities = list(
        rank_constrained = TRUE,
        whitening = "shared_ar",
        se_modes = c("conditional", "bootstrap"),
        contrast_scope = "event_only"
      )
    )
  }

  invisible(TRUE)
}
