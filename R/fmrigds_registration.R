# fmrireg's contributions to the fmrigds registries. Everything is registered
# once, at load, under a namespaced name, and never replaces an fmrigds name.
#
# - reducer "fmrireg:meta_reg": the meta-regression estimator behind
#   fmri_meta() (FE / DerSimonian-Laird / Paule-Mandel, Huber IRLS, and the
#   relative "equal"/"custom" weights). It wraps .meta_fit_core() unchanged.
# - posthoc "fmrireg:fdr:spatial": an alias of fmrigds's "fdr:sabha", which is
#   the port of spatial_fdr().

.fmrireg_meta_reg_name <- "fmrireg:meta_reg"
.fmrireg_posthoc_name <- "fmrireg:fdr:spatial"

# Reducer kernel. fmrigds passes one contrast at a time: `beta` and `var` are
# subjects x features, `X` is the design built from the plan's formula.
# Options: method ("pm", "fe", "dl"; "reml" runs PM), robust ("none",
# "huber"; "t" is accepted and fits without robust weights, as fmri_meta()
# always has), weights ("1/var" = inverse variance with known sampling
# variances; "equal"/"custom" = relative WLS weights, see .meta_fit_core),
# custom_weights (fmri_meta's weights_custom; named so that it cannot
# partially match `weights`), huber_c, robust_iter, n_threads.
.fmrireg_meta_reg_kernel <- function(beta, var, X, z, p, df, df1, df2, opts) {
  Y <- as.matrix(beta)
  if (is.null(X)) {
    X <- matrix(1, nrow(Y), 1L, dimnames = list(NULL, "(Intercept)"))
  }
  weights <- opts[["weights"]] %||% "1/var"
  if (identical(weights, "1/var")) weights <- "ivw"
  fit <- .meta_fit_core(
    Y = Y, V = var, X = X,
    method = .meta_resolve_method(opts[["method"]] %||% "pm"),
    robust = opts[["robust"]] %||% "none",
    weights = weights, weights_custom = opts[["custom_weights"]],
    huber_c = opts[["huber_c"]] %||% 1.345, robust_iter = opts[["robust_iter"]] %||% 2,
    n_threads = opts[["n_threads"]] %||% .fmrireg_num_threads()
  )
  terms <- colnames(X) %||% paste0("X", seq_len(ncol(X)))
  coef_rows <- function(M) {
    M <- as.matrix(M)
    rownames(M) <- terms
    M
  }
  # I2 is undefined when Q = 0 (an exact fit); report 0 there.
  Q <- as.numeric(fit$Q_fe)
  df_res <- as.numeric(fit$df)
  I2 <- as.numeric(fit$I2_fe)
  fill <- is.na(I2) & !is.na(Q) & !is.na(df_res)
  I2[fill] <- pmax((Q[fill] - df_res[fill]) / pmax(Q[fill], .Machine$double.eps), 0)
  list(
    coef = coef_rows(fit$beta),
    se_coef = coef_rows(fit$se),
    z_coef = coef_rows(fit$z),
    p_coef = coef_rows(2 * stats::pnorm(abs(fit$z), lower.tail = FALSE)),
    cov_tri = fit$cov_tri,
    tau2 = as.numeric(fit$tau2),
    Q = Q,
    I2 = I2,
    df_res = df_res
  )
}

.register_fmrireg_meta_reg <- function() {
  if (!requireNamespace("fmrigds", quietly = TRUE) ||
      !is.null(fmrigds::get_reducer(.fmrireg_meta_reg_name))) {
    return(invisible(FALSE))
  }
  fmrigds::register_reducer(
    name = .fmrireg_meta_reg_name,
    fun = .fmrireg_meta_reg_kernel,
    requires = c("beta", "var"),
    provides = c("coef", "se_coef", "z_coef", "p_coef", "cov_tri",
                 "tau2", "Q", "I2", "df_res"),
    options_schema = list(
      weights = c("1/var", "equal", "custom"),
      method = c("pm", "fe", "dl", "reml"),
      robust = c("none", "huber", "t")
    ),
    model_contract = list(uses_X = TRUE, estimands = "linear",
                          weight_mode = "model_specific",
                          missingness = "samplewise")
  )
  invisible(TRUE)
}

.register_fmrireg_posthoc <- function() {
  if (!requireNamespace("fmrigds", quietly = TRUE) ||
      !is.null(fmrigds::get_posthoc(.fmrireg_posthoc_name))) {
    return(invisible(FALSE))
  }
  sabha <- fmrigds::get_posthoc("fdr:sabha")
  if (is.null(sabha)) return(invisible(FALSE))
  fmrigds::register_posthoc(
    name = .fmrireg_posthoc_name,
    fun = sabha$fun,
    requires = sabha$requires,
    provides = sabha$provides,
    overwrite = FALSE
  )
  invisible(TRUE)
}
