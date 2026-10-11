# Chunked voxelwise contrasts + pull_stat invalid-index warning + preproject edges.

test_that("pull_stat warns on invalid event indices and recovers columns", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  fit2 <- fit
  b <- fit2$result$betas$data[[1]]$estimate[[1]]
  fit2$result$event_indices <- c(ncol(b) + 5L, ncol(b) + 6L)
  expect_warning(
    out <- fmrireg:::pull_stat(fit2, "betas", "estimate"),
    "No valid event indices|valid"
  )
  expect_true(inherits(out, "tbl_df") || is.data.frame(out) || is.matrix(out))
})

test_that(".fast_preproject rank-deficient aliased path", {
  X_def <- cbind(1, c(0, 1, 0, 1, 0, 1), c(0, 1, 0, 1, 0, 1))
  colnames(X_def) <- c("Intercept", "A", "A_dup")
  expect_warning(
    proj <- fmrireg:::.fast_preproject(X_def),
    "rank deficient|Aliased"
  )
  expect_true(is.list(proj))
  expect_true(!is.null(proj$XtXinv))
  expect_false(isTRUE(proj$is_full_rank))
})
