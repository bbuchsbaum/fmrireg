# Review 2026-10, stream C, C4: post-hoc fit_contrasts.fmri_lm() resolved
# contrast specs against the first event term only, fell back to positional
# column indices when names did not match, and silently dropped (with a
# warning) every spec that was not a pair_contrast.

.posthoc_fixture <- function() {
  set.seed(1)
  n <- 160L
  ne <- 24L
  etab <- data.frame(
    onset = seq(4, by = 6, length.out = ne),
    kind = rep(c("c", "t"), ne / 2),
    cond = factor(rep(c("A", "B", "C"), length.out = ne)),
    task = factor(rep(c("x", "x", "y", "y"), length.out = ne)),
    run = 1L
  )
  Y <- matrix(rnorm(n * 4), n, 4)
  dset <- matrix_frame(Y, TR = 1, run_length = n, event_table = etab)
  specs <- list(
    AvB = pair_contrast(~ cond == "A", ~ cond == "B", name = "AvB"),
    unitc = unit_contrast(~ cond, name = "unitc"),
    maincond = oneway_contrast(~ cond, name = "maincond"),
    XvY = pair_contrast(~ task == "x", ~ task == "y", name = "XvY")
  )
  c1 <- do.call(contrast_set, unname(specs[c("AvB", "unitc", "maincond")]))
  c2 <- contrast_set(specs$XvY)
  fit <- fmri_lm(
    onset ~ hrf(cond, subset = kind == "c", contrasts = c1) +
      hrf(task, subset = kind == "t", contrasts = c2),
    block = ~ run, dataset = dset
  )
  list(fit = fit, specs = specs, Y = Y)
}

test_that("C4: post-hoc contrasts equal in-fit contrasts for pair, unit and oneway specs on several terms", {
  s <- .posthoc_fixture()
  fit <- s$fit
  ph <- fit_contrasts(fit, s$specs)
  expect_setequal(names(ph), names(s$specs))

  tst <- stats(fit, type = "contrasts")
  est <- coef(fit, type = "contrasts")
  se <- standard_error(fit, type = "contrasts")
  pv <- p_values(fit, type = "contrasts")
  for (nm in c("AvB", "unitc", "XvY")) {
    expect_equal(ph[[nm]]$type, "contrast")
    expect_equal(ph[[nm]]$estimate, est[[nm]], tolerance = 1e-8)
    expect_equal(ph[[nm]]$se, se[[nm]], tolerance = 1e-8)
    expect_equal(ph[[nm]]$stat, tst[[nm]], tolerance = 1e-8)
    expect_equal(ph[[nm]]$prob, pv[[nm]], tolerance = 1e-8)
  }
  expect_equal(ph$maincond$type, "Fcontrast")
  expect_equal(ph$maincond$stat, stats(fit, type = "F")[["maincond"]], tolerance = 1e-8)

  # Independent reference for the second-term contrast (term 2, not term 1).
  X <- as.matrix(design_matrix(fit$model))
  jx <- which(colnames(X) == "task_task.x")
  jy <- which(colnames(X) == "task_task.y")
  for (v in seq_len(ncol(s$Y))) {
    m <- lm(s$Y[, v] ~ 0 + X)
    l <- numeric(ncol(X)); l[jx] <- 1; l[jy] <- -1
    tv <- sum(l * coef(m)) / sqrt(drop(t(l) %*% vcov(m) %*% l))
    expect_equal(ph$XvY$stat[v], tv, tolerance = 1e-6)
  }
})

test_that("C4: named weight vectors map by design column name, not position", {
  s <- .posthoc_fixture()
  fit <- s$fit
  w <- c(task_task.x = 1, task_task.y = -1)
  ph <- fit_contrasts(fit, list(XvY = w))
  expect_equal(ph$XvY$stat, stats(fit, type = "contrasts")[["XvY"]], tolerance = 1e-8)

  # Order of the names does not matter.
  ph2 <- fit_contrasts(fit, list(XvY = w[2:1]))
  expect_equal(ph2$XvY$stat, ph$XvY$stat)
})

test_that("C4: unmatched names and unresolvable specs are errors, not silent drops", {
  s <- .posthoc_fixture()
  fit <- s$fit
  expect_error(fit_contrasts(fit, list(bad = c(nope_a = 1, nope_b = -1))),
               "bad")
  # Unnamed weights shorter than the coefficient vector have no defined
  # mapping; the old code applied them to the first columns.
  expect_error(fit_contrasts(fit, list(short = c(1, -1))), "short")
  # A spec naming a variable no event term has is an error ...
  expect_error(
    fit_contrasts(fit, list(novar = pair_contrast(~ colour == "red", ~ colour == "blue",
                                                  name = "novar"))),
    "novar"
  )
  # ... and an empty level mask is reported, as it is for in-fit contrasts.
  expect_warning(
    fit_contrasts(fit, list(nolevel = pair_contrast(~ cond == "Z", ~ cond == "A",
                                                    name = "nolevel"))),
    "empty"
  )
})
