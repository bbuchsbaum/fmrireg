# robust = "t" used to be mapped silently to "huber"; it is now refused
# (see .meta_match_robust()).
.fmrigds_map_robust <- function(r) {
  .meta_match_robust(r)
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

# The assay fields each fmrireg reducer declares (`provides`) and returns.
# fmrigds pre-allocates one array per declared non-regression field and fills
# it from the reducer output, and silently ignores undeclared fields, so the
# declaration and the output must agree exactly. The packed coefficient
# covariance of the regression reducers is not an assay: fmrigds carries a
# reducer's `cov_tri` element (K(K+1)/2 x features) as cov:<a>:<b> parameter
# arrays plus an attachment. Declaring it as an assay made fmrigds recycle the
# packed matrix into a meaningless "cov_tri" assay.
.fmrigds_reducer_provides <- function(name) {
  is_reg <- grepl("_reg$", name)
  is_fe <- grepl("(^|:)fe(_reg)?$", name)
  base <- if (is_reg) c("coef", "se_coef", "t_coef", "p_coef") else c("beta", "se", "z", "p")
  het <- if (is_fe) c("I2", "Q", "df") else c("tau2", "I2", "Q", "df")
  c(base, het)
}

.fmrigds_stop_reducer_contrasts <- function() {
  stop("fmrireg's fmrigds meta reducers do not compute fit-time contrasts ",
       "(fmrigds cannot carry them as reducer outputs); use ",
       "fmri_meta(..., contrasts = ) or contrast() on the fitted model.",
       call. = FALSE)
}

.fmrigds_make_meta_reducer <- function(method_id,
                                       name = if (identical(method_id, "fe")) "meta:fe" else paste0("meta:re:", method_id)) {
  force(method_id)
  provides <- .fmrigds_reducer_provides(name)
  function(beta, var, X, z, p, df, df1, df2, opts) {
    opts <- opts %||% list()
    if (!is.null(opts$contrasts)) .fmrigds_stop_reducer_contrasts()
    V <- .fmrigds_coerce_var(var, opts$se_override %||% NULL)
    robust <- .fmrigds_map_robust(opts$robust)
    method <- .fmrigds_meta_method(method_id, opts$method)
    res <- fmrireg::fmri_meta_fit(
      Y = beta, V = V, X = X,
      method = method,
      robust = robust,
      n_threads = getOption("fmrireg.num_threads", 0)
    )
    out <- list(
      beta = t(res$beta),
      se   = t(res$se),
      z    = t(res$z),
      p    = 2 * stats::pnorm(abs(t(res$z)), lower.tail = FALSE),
      tau2 = as.numeric(res$tau2),
      I2   = as.numeric(res$I2_fe),
      Q    = as.numeric(res$Q_fe),
      df   = as.numeric(res$df)
    )
    out[provides]
  }
}

.fmrigds_make_meta_reg <- function(method_default,
                                   name = if (identical(method_default, "fe")) "meta:fe_reg" else "meta:re_reg") {
  force(method_default)
  provides <- .fmrigds_reducer_provides(name)
  function(beta, var, X, z, p, df, df1, df2, opts) {
    opts <- opts %||% list()
    if (!is.null(opts$contrasts)) .fmrigds_stop_reducer_contrasts()
    method <- .fmrigds_meta_method(method_default, opts$method)
    robust <- .fmrigds_map_robust(opts$robust)
    # The covariance is always returned (whatever opts$return_cov says), so
    # contrast() on the result is exact for every method.
    res <- fmrireg::fmri_meta_fit_cov(
      Y = beta, V = var, X = X,
      method = method,
      robust = robust,
      n_threads = getOption("fmrireg.num_threads", 0)
    )
    out <- list(
      coef    = res$beta,
      se_coef = res$se,
      t_coef  = res$z,
      p_coef  = 2 * stats::pnorm(abs(res$z), lower.tail = FALSE),
      tau2    = as.numeric(res$tau2),
      I2      = as.numeric(res$I2_fe),
      Q       = as.numeric(res$Q_fe),
      df      = as.numeric(res$df)
    )
    out <- out[provides]
    out$cov_tri <- res$cov_tri
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


.fmrigds_meta_reducer_names <- c(
  "meta:fe", "meta:re:dl", "meta:re:pm", "meta:re:reml",
  "frg:meta:fe", "frg:meta:re:dl", "frg:meta:re:pm", "frg:meta:re:reml"
)
.fmrigds_meta_reg_names <- c("meta:fe_reg", "meta:re_reg")
.fmrigds_required_reducers <- c(.fmrigds_meta_reducer_names, .fmrigds_meta_reg_names)

# Register fmrireg's reducers and the fdr:spatial posthoc with fmrigds
# (overwriting any existing registrations of these names). Errors propagate;
# .ensure_fmrigds_registered() reports them.
.register_fmrireg_bindings <- function() {
  if (!requireNamespace("fmrigds", quietly = TRUE)) return(invisible(FALSE))

  for (nm in .fmrigds_meta_reducer_names) {
    method_id <- sub("^.*:", "", nm)
    fmrigds::register_reducer(
      name = nm,
      fun = .fmrigds_make_meta_reducer(method_id, name = nm),
      requires = c("beta", "var"),
      provides = .fmrigds_reducer_provides(nm)
    )
  }
  for (nm in .fmrigds_meta_reg_names) {
    method_default <- if (nm == "meta:re_reg") "pm" else "fe"
    fmrigds::register_reducer(
      name = nm,
      fun = .fmrigds_make_meta_reg(method_default, name = nm),
      requires = c("beta", "var", "X"),
      provides = .fmrigds_reducer_provides(nm)
    )
  }
  fmrigds::register_posthoc(
    name = "fdr:spatial",
    fun = .fmrigds_spatial_posthoc,
    requires = c("p"),
    provides = c("q"),
    overwrite = TRUE
  )

  invisible(TRUE)
}
