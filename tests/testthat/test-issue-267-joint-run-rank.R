# Joint estimation checks rank on the joint design only (bbuchsbaum/fmrireg#267).
# The per-run fits that only estimate the noise model used to warn that a
# condition absent from one run made "the design matrix" rank deficient.

rank_binding <- function(seed = 1) {
  set.seed(seed)
  ev <- data.frame(onset = rep(c(5, 25, 45, 65), 2),
                   cond = factor(c("A", "B", "C", "A", "B", "A", "B", "A")),
                   run = rep(1:2, each = 4))
  list(id = "s1", scans = matrix(rnorm(160 * 5), 160, 5), TR = 2,
       run_length = c(80L, 80L), events = ev)
}

rank_warnings <- function(expr) {
  msgs <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    if (grepl("rank deficient", conditionMessage(w))) msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}

fit_267 <- function(control, formula = onset ~ hrf(cond), b = rank_binding()) {
  rank_warnings(run_job(instantiate(fmri_template(formula, ~ run, control = control), b)))
}

test_that("a level missing from one run raises no rank warning under joint estimation", {
  J <- estimation_spec("joint")
  controls <- list(
    ar1 = fmri_lm_control(estimation = J, noise = noise_spec("ar1")),
    ar2_global = fmri_lm_control(estimation = J, noise = noise_spec("ar2", pooling = "global")),
    robust = fmri_lm_control(estimation = J, robust = robust_spec("huber")),
    ar1_robust = fmri_lm_control(estimation = J, noise = noise_spec("ar1"),
                                 robust = robust_spec("huber"))
  )
  for (nm in names(controls)) {
    res <- fit_267(controls[[nm]])
    expect_length(res$warnings, 0L)
    expect_true(all(is.finite(coef(res$value))), info = nm)
    # the joint design is full rank
    X <- as.matrix(design_matrix(res$value$model))
    expect_equal(qr(X)$rank, ncol(X), info = nm)
  }
})

test_that("trialwise joint fits are not reported rank deficient per run", {
  res <- fit_267(fmri_lm_control(estimation = estimation_spec("joint"),
                                 noise = noise_spec("ar1")),
                 formula = onset ~ trialwise())
  expect_length(res$warnings, 0L)
  expect_true(all(is.finite(coef(res$value))))
})

test_that("a rank-deficient joint design still warns, once", {
  b <- rank_binding()
  b$events$dup <- b$events$cond
  res <- fit_267(fmri_lm_control(estimation = estimation_spec("joint"),
                                 noise = noise_spec("ar1")),
                 formula = onset ~ hrf(cond) + hrf(dup), b = b)
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "Aliased columns: dup")
})

test_that("runwise_meta still reports the run whose design is deficient", {
  res <- fit_267(fmri_lm_control(estimation = estimation_spec("runwise_meta"),
                                 noise = noise_spec("ar1")))
  expect_gt(length(res$warnings), 0L)
  expect_match(res$warnings[1], "Aliased columns: .*C")
})
