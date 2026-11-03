# fmri_template(contrasts = ) reach the fit, the contrast table and
# write_results() (bbuchsbaum/fmrireg#265). Before the fix they were computed
# post hoc into an attribute only reduce_contrasts() read, and the default
# runwise_meta + AR configuration refused them outright.

two_run_binding <- function(seed = 1) {
  set.seed(seed)
  ev <- data.frame(onset = rep(c(5, 25, 45, 65), 2),
                   cond = factor(rep(c("A", "B", "A", "B"), 2)),
                   run = rep(1:2, each = 4))
  list(id = "sub-01", scans = matrix(rnorm(160 * 4), 160, 4), TR = 2,
       run_length = c(80L, 80L), events = ev)
}

cs_265 <- contrast_set(
  pair_contrast(~ cond == "A", ~ cond == "B", name = "AvB"),
  contrast(~ A - B, name = "AminusB")
)

ctl_265 <- function(scope) {
  fmri_lm_control(estimation = estimation_spec(scope), noise = noise_spec("ar1"))
}

test_that("template contrasts equal the same contrasts written in the formula", {
  b <- two_run_binding()
  for (scope in c("runwise_meta", "joint")) {
    tmpl <- fmri_template(onset ~ hrf(cond), ~ run, contrasts = cs_265,
                          control = ctl_265(scope))
    in_formula <- fmri_template(onset ~ hrf(cond, contrasts = cs_265), ~ run,
                                control = ctl_265(scope))
    fit <- run_job(instantiate(tmpl, b))
    ref <- run_job(instantiate(in_formula, b))
    expect_setequal(fit$result$contrasts$name, c("AvB", "AminusB"))
    expect_equal(coef(fit, type = "contrasts"), coef(ref, type = "contrasts"),
                 info = scope)
    expect_equal(standard_error(fit, type = "contrasts"),
                 standard_error(ref, type = "contrasts"), info = scope)
    # the two spellings of A - B agree
    expect_equal(coef(fit, type = "contrasts")$AvB,
                 coef(fit, type = "contrasts")$AminusB, info = scope)
  }
})

test_that("reduce_contrasts reports template contrasts with se and stat", {
  b <- two_run_binding()
  tmpl <- fmri_template(onset ~ hrf(cond), ~ run, contrasts = cs_265,
                        control = ctl_265("runwise_meta"),
                        reducer = reduce_contrasts())
  out <- run_job(instantiate(tmpl, b))
  expect_setequal(unique(out$term), c("AvB", "AminusB"))
  expect_equal(nrow(out), 2L * 4L)
  expect_true(all(is.finite(out$se)) && all(out$se > 0))
  expect_equal(out$stat, out$estimate / out$se, tolerance = 1e-8)
})

test_that("a weight-vector template contrast is computed post hoc into the contrast table", {
  b <- two_run_binding()
  w <- c(cond_cond.A = 1, cond_cond.B = -1)
  fit <- run_job(instantiate(
    fmri_template(onset ~ hrf(cond), ~ run, contrasts = list(w = w),
                  control = fmri_lm_control(estimation = estimation_spec("joint"))),
    b))
  ref <- run_job(instantiate(
    fmri_template(onset ~ hrf(cond, contrasts = cs_265[1]), ~ run,
                  control = fmri_lm_control(estimation = estimation_spec("joint"))),
    b))
  expect_equal(fit$result$contrasts$name, "w")
  expect_equal(coef(fit, type = "contrasts")$w, coef(ref, type = "contrasts")$AvB)
  expect_equal(standard_error(fit, type = "contrasts")$w,
               standard_error(ref, type = "contrasts")$AvB)
})

test_that("a template contrast that names a formula contrast is an error", {
  b <- two_run_binding()
  tmpl <- fmri_template(onset ~ hrf(cond, contrasts = cs_265), ~ run,
                        contrasts = cs_265[1])
  expect_error(build_model(instantiate(tmpl, b)), "duplicates a contrast")
})

test_that("contrast(~ A - B) resolves post hoc against level labels", {
  b <- two_run_binding()
  fit <- run_job(instantiate(
    fmri_template(onset ~ hrf(cond), ~ run,
                  control = fmri_lm_control(estimation = estimation_spec("joint"))),
    b))
  ph <- fit_contrasts(fit, list(contrast(~ A - B, name = "d")))
  ph2 <- fit_contrasts(fit, list(pair_contrast(~ cond == "A", ~ cond == "B", name = "d")))
  expect_equal(ph$d$estimate, ph2$d$estimate)
  expect_equal(ph$d$se, ph2$d$se)
})

test_that("reduce_write_results writes template contrast maps (runwise_meta + AR)", {
  skip_if_not_installed("bidser")
  skip_if_not_installed("RNifti")
  skip_if_not_installed("jsonlite")
  root <- withr::local_tempdir()
  suppressMessages(make_bids_fixture(root, subjects = "01"))
  proj <- suppressMessages(bidser::bids_project(root, fmriprep = TRUE))
  mani <- from_bids(proj, task = "stroop", space = "MNI152NLin2009cAsym",
                    mask = file.path(root, "brain_mask.nii.gz"))
  out <- file.path(root, "derivatives", "glm")
  tmpl <- fmri_template(
    onset ~ hrf(trial_type), ~ run,
    contrasts = contrast_set(pair_contrast(~ trial_type == "incongruent",
                                           ~ trial_type == "congruent", name = "ivc")),
    control = ctl_265("runwise_meta"),
    reducer = reduce_write_results("nifti", stats = c("beta", "se", "tstat"), path = out)
  )
  res <- run_jobs(instantiate(tmpl, mani), on_error = "stop")
  expect_true(all(res$ok))
  for (stat in c("beta", "se", "tstat")) {
    expect_length(list.files(out, sprintf("^sub-01_.*_desc-%s_bold\\.nii\\.gz$", stat)), 1L)
  }
})
