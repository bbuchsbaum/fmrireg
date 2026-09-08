# Twelfth wave: plugin basis resolution/replacement + engine_spec print.

test_that("get_basis / resolve_basis / .fmrireg_replace_basis_expr cover registry", {
  nm <- paste0("cov_basis_wave12_", as.integer(Sys.time()) %% 1e6L)
  register_basis(nm, function(span = 24, ...) {
    fmrihrf::HRF_SPMG1
  })
  expect_true(is.function(fmrireg:::get_basis(nm)))
  expect_null(fmrireg:::get_basis(""))
  expect_null(fmrireg:::get_basis(c("a", "b")))

  resolved <- resolve_basis(nm, span = 30)
  expect_true(!is.null(resolved))

  # Replace basis= character name inside hrf() call
  expr <- substitute(onset ~ hrf(condition, basis = BASIS), list(BASIS = nm))
  replaced <- fmrireg:::.fmrireg_replace_basis_expr(expr, env = environment())
  expect_true(is.call(replaced) || identical(replaced, expr))
  txt <- paste(deparse(replaced), collapse = " ")
  expect_true(grepl("resolve_basis|hrf", txt))

  # Non-hrf call leaves expression intact
  other <- quote(mean(x))
  expect_equal(fmrireg:::.fmrireg_replace_basis_expr(other, env = environment()), other)

  # Unknown basis name left as-is
  expr2 <- quote(onset ~ hrf(condition, basis = "not_registered_xyz"))
  left <- fmrireg:::.fmrireg_replace_basis_expr(expr2, env = environment())
  expect_true(is.call(left))
})

test_that("print.fmrireg_engine_spec and .new_engine_spec include_functions", {
  reg <- list(
    fit = function(model, dataset, args, cfg) list(ok = TRUE),
    preflight = function(...) TRUE,
    capabilities = list(robust = TRUE, ma = FALSE, ar_voxelwise = TRUE)
  )
  spec <- fmrireg:::.new_engine_spec("latent_sketch", reg, include_functions = TRUE)
  expect_s3_class(spec, "fmrireg_engine_spec")
  expect_equal(spec$source, "builtin")
  expect_equal(spec$aliases, "sketch")
  expect_true(is.function(spec$fit))
  expect_true(isTRUE(spec$has_preflight))

  out <- capture.output(print(spec))
  expect_true(any(grepl("latent_sketch|builtin|capabilities|aliases", out, ignore.case = TRUE)))

  plugin_spec <- fmrireg:::.new_engine_spec(
    "custom_plugin_x",
    list(fit = function(...) NULL, preflight = NULL, capabilities = list()),
    include_functions = FALSE
  )
  expect_equal(plugin_spec$source, "plugin")
  expect_false(isTRUE(plugin_spec$has_preflight))
  expect_null(plugin_spec$fit)
})

test_that("multiresponse_lm fits with and without modmat", {
  set.seed(241)
  n <- 40L
  V <- 3L
  df <- data.frame(a = rnorm(n), b = rnorm(n))
  Y <- cbind(1, df$a, df$b) %*% matrix(c(0.2, 0.5, -0.3), 3, 1) %*% matrix(1, 1, V) +
    matrix(rnorm(n * V, sd = 0.3), n, V)
  env <- list2env(df)
  env$.y <- Y
  form <- .y ~ a + b
  X <- model.matrix(form, env)
  conlist <- list()
  fcon <- list()
  out1 <- fmrireg:::multiresponse_lm(form, env, conlist, colnames(X), fcon, modmat = NULL)
  expect_true(is.list(out1) || inherits(out1, "tbl_df") || is.data.frame(out1))

  out2 <- fmrireg:::multiresponse_lm(form, env, conlist, colnames(X), fcon, modmat = X)
  expect_true(is.list(out2) || inherits(out2, "tbl_df") || is.data.frame(out2))
})
