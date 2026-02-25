.fmrigds_map_robust <- function(r) {
  if (identical(r, "t")) "huber" else r %||% "none"
}

.fmrigds_coerce_var <- function(var, se) {
  if (!is.null(var)) return(var)
  if (!is.null(se)) return(se^2)
  stop("Reducer requires 'var' or 'se' assay", call. = FALSE)
}

.fmrigds_meta_method <- function(method_id, opts_method = NULL) {
  if (!is.null(opts_method)) return(opts_method)
  switch(method_id,
         fe = "fe",
         dl = "dl",
         pm = "pm",
         reml = "reml",
         method_id)
}

.fmrigds_make_meta_reducer <- function(method_id) {
  force(method_id)
  function(beta, var, X, z, p, df, df1, df2, opts) {
    opts <- opts %||% list()
    V <- .fmrigds_coerce_var(var, opts$se_override %||% NULL)
    robust <- .fmrigds_map_robust(opts$robust)
    method <- .fmrigds_meta_method(method_id, opts$method)
    if (!is.null(opts$contrasts)) {
      Cmat <- opts$contrasts
      res <- fmrireg::fmri_meta_fit_contrasts(
        Y = beta, V = V, X = X, Cmat = Cmat,
        method = method,
        robust = robust,
        n_threads = getOption("fmrireg.num_threads", 0)
      )
      out <- list(
        beta = t(res$beta),
        se   = t(res$se),
        z    = t(res$z),
        p    = 2 * stats::pnorm(abs(t(res$z)), lower.tail = FALSE),
        con_est = t(res$c_beta),
        con_se  = t(res$c_se),
        con_z   = t(res$c_z)
      )
    } else {
      return_cov <- isTRUE(identical(opts$return_cov, "tri")) && method %in% c("pm", "reml")
      if (return_cov) {
        res <- fmrireg::fmri_meta_fit_cov(
          Y = beta, V = V, X = X,
          method = method,
          robust = robust,
          n_threads = getOption("fmrireg.num_threads", 0)
        )
      } else {
        res <- fmrireg::fmri_meta_fit(
          Y = beta, V = V, X = X,
          method = method,
          robust = robust,
          n_threads = getOption("fmrireg.num_threads", 0)
        )
      }

      out <- list(
        beta = t(res$beta),
        se   = t(res$se),
        z    = t(res$z),
        p    = 2 * stats::pnorm(abs(t(res$z)), lower.tail = FALSE)
      )
      if (!is.null(res$tau2)) out$tau2 <- as.numeric(res$tau2)
      if (!is.null(res$I2_fe)) out$I2   <- as.numeric(res$I2_fe)
      if (!is.null(res$Q_fe)) out$Q     <- as.numeric(res$Q_fe)
      if (!is.null(res$df))   out$df    <- as.numeric(res$df)
      if (return_cov && !is.null(res$cov_tri)) out$cov_tri <- t(res$cov_tri)
    }
    out
  }
}

.fmrigds_make_meta_reg <- function(method_default) {
  force(method_default)
  function(beta, var, X, z, p, df, df1, df2, opts) {
    opts <- opts %||% list()
    method <- .fmrigds_meta_method(method_default, opts$method)
    robust <- .fmrigds_map_robust(opts$robust)
    return_cov <- isTRUE(identical(opts$return_cov, "tri")) && method %in% c("pm", "reml")
    if (return_cov) {
      res <- fmrireg::fmri_meta_fit_cov(
        Y = beta, V = var, X = X,
        method = method,
        robust = robust,
        n_threads = getOption("fmrireg.num_threads", 0)
      )
    } else {
      res <- fmrireg::fmri_meta_fit(
        Y = beta, V = var, X = X,
        method = method,
        robust = robust,
        n_threads = getOption("fmrireg.num_threads", 0)
      )
    }
    out <- list(
      coef   = res$beta,
      se_coef= res$se,
      t_coef = res$z,
      p_coef = 2 * stats::pnorm(abs(res$z), lower.tail = FALSE),
      tau2   = if (!is.null(res$tau2)) as.numeric(res$tau2) else NULL,
      I2     = if (!is.null(res$I2_fe)) as.numeric(res$I2_fe) else NULL,
      Q      = if (!is.null(res$Q_fe))  as.numeric(res$Q_fe)  else NULL,
      df     = if (!is.null(res$df))    as.numeric(res$df)    else NULL
    )
    if (return_cov && !is.null(res$cov_tri)) out$cov_tri <- res$cov_tri
    out
  }
}

.fmrigds_spatial_posthoc_compute <- function(arrays, group, alpha) {
  # Get p or z from either top-level arrays or arrays$assays
  p_arr <- arrays[["p"]]
  if (is.null(p_arr) && !is.null(arrays$assays)) p_arr <- arrays$assays[["p"]]
  z_arr <- arrays[["z"]]
  if (is.null(z_arr) && !is.null(arrays$assays)) z_arr <- arrays$assays[["z"]]

  # Normalize to [P x S x K]
  if (is.null(p_arr) && is.null(z_arr)) {
    stop("fdr:spatial requires p (or z to derive p)", call. = FALSE)
  }
  if (!is.null(z_arr)) {
    z_dims <- dim(z_arr)
    if (is.null(z_dims)) z_arr <- array(z_arr, dim = c(length(z_arr), 1L, 1L))
    if (length(dim(z_arr)) == 2L) z_arr <- array(z_arr, dim = c(dim(z_arr)[1L], 1L, dim(z_arr)[2L]))
    # keep z for empirical-null handling inside spatial_fdr(z=...)
  } else {
    p_dims <- dim(p_arr)
    if (is.null(p_dims)) p_arr <- array(p_arr, dim = c(length(p_arr), 1L, 1L))
    if (length(dim(p_arr)) == 2L) p_arr <- array(p_arr, dim = c(dim(p_arr)[1L], 1L, dim(p_arr)[2L]))
  }

  dims <- if (!is.null(z_arr)) dim(z_arr) else dim(p_arr)
  if (is.null(dims) || length(dims) != 3L) {
    stop("fdr:spatial expected assay dimensions [features x subjects x contrasts]", call. = FALSE)
  }
  Q <- array(NA_real_, dim = dims)
  for (j in seq_len(dims[2L])) {
    for (k in seq_len(dims[3L])) {
      if (!is.null(z_arr)) {
        zv <- z_arr[, j, k]
        sres <- fmrireg::spatial_fdr(z = zv, group = group, alpha = alpha)
      } else {
        pv <- p_arr[, j, k]
        sres <- fmrireg::spatial_fdr(p = pv, group = group, alpha = alpha)
      }
      Q[, j, k] <- sres$q
    }
  }
  Q
}

.fmrigds_spatial_posthoc <- function(arrays, opts) {
  group <- opts$group
  if (is.null(group)) stop("fdr:spatial requires options$group (length n_features)", call. = FALSE)
  alpha <- opts$alpha %||% 0.05
  q <- .fmrigds_spatial_posthoc_compute(arrays, group, alpha)
  list(q = q)
}

.fmrigds_required_reducers <- c(
  "meta:fe", "meta:re:dl", "meta:re:pm", "meta:re:reml",
  "frg:meta:fe", "frg:meta:re:dl", "frg:meta:re:pm", "frg:meta:re:reml",
  "meta:fe_reg", "meta:re_reg"
)

.fmrigds_check_reducer <- function(name, expected) {
  ns <- asNamespace("fmrigds")
  get_reducer <- get("get_reducer", envir = ns)
  r <- try(get_reducer(name), silent = TRUE)
  if (inherits(r, "try-error") || is.null(r) || is.null(r$provides)) return(FALSE)
  setequal(sort(r$provides), sort(expected))
}

.register_fmrireg_bindings <- function(force = FALSE) {
  if (!requireNamespace("fmrigds", quietly = TRUE)) return(invisible(FALSE))

  reducers <- try(fmrigds::list_reducers(), silent = TRUE)
  need_reducers <- force || inherits(reducers, "try-error") ||
    !all(.fmrigds_required_reducers %in% reducers) ||
    !all(c(
      .fmrigds_check_reducer("meta:fe_reg", c("coef","se_coef","t_coef","p_coef","I2","Q","df")),
      .fmrigds_check_reducer("meta:re_reg", c("coef","se_coef","t_coef","p_coef","tau2","I2","Q","df","cov_tri"))
    ))

  # Always ensure our posthoc is registered; reducers only when needed

  meta_reducers <- list(
    "meta:fe"        = .fmrigds_make_meta_reducer("fe"),
    "meta:re:dl"     = .fmrigds_make_meta_reducer("dl"),
    "meta:re:pm"     = .fmrigds_make_meta_reducer("pm"),
    "meta:re:reml"   = .fmrigds_make_meta_reducer("reml"),
    "frg:meta:fe"    = .fmrigds_make_meta_reducer("fe"),
    "frg:meta:re:dl" = .fmrigds_make_meta_reducer("dl"),
    "frg:meta:re:pm" = .fmrigds_make_meta_reducer("pm"),
    "frg:meta:re:reml" = .fmrigds_make_meta_reducer("reml")
  )

  meta_reg_reducers <- list(
    "meta:fe_reg" = .fmrigds_make_meta_reg("fe"),
    "meta:re_reg" = .fmrigds_make_meta_reg("pm")
  )

  if (need_reducers) {
    for (nm in names(meta_reducers)) {
      provides <- if (grepl("re:(pm|reml)", nm)) {
        c("beta","se","z","p","tau2","I2","Q","df","cov_tri")
      } else if (grepl("re:dl", nm, fixed = TRUE)) {
        c("beta","se","z","p","tau2","I2","Q","df")
      } else {
        c("beta","se","z","p","I2","Q","df")
      }
      try(
        fmrigds::register_reducer(
          name = nm,
          fun = meta_reducers[[nm]],
          requires = c("beta", "var"),
          provides = provides
        ),
        silent = TRUE
      )
    }
    for (nm in names(meta_reg_reducers)) {
      provides <- if (nm == "meta:re_reg") c("coef","se_coef","t_coef","p_coef","tau2","I2","Q","df","cov_tri") else c("coef","se_coef","t_coef","p_coef","I2","Q","df")
      try(
        fmrigds::register_reducer(
          name = nm,
          fun = meta_reg_reducers[[nm]],
          requires = c("beta", "var", "X"),
          provides = provides
        ),
        silent = TRUE
      )
    }
  }

  try(
    fmrigds::register_posthoc(
      name = "fdr:spatial",
      fun = .fmrigds_spatial_posthoc,
      requires = c("p"),
      provides = c("q"),
      overwrite = TRUE
    ),
    silent = TRUE
  )

  invisible(TRUE)
}
