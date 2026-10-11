# Performance pass 2, P2-2: coef_images() resolves the coefficient names
# (design_matrix()) and the mask space once, not once per coefficient.

.p2_coef_fixture <- function(spatial = TRUE) {
  set.seed(42)
  dims <- c(4, 3, 2)
  n <- 60
  etab <- data.frame(
    onset = c(4, 14, 24, 34, 44, 52),
    condition = factor(rep(c("A", "B", "C"), 2)),
    run = 1L
  )
  if (spatial) {
    arr <- array(rnorm(prod(dims) * n), c(dims, n))
    vec <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, n)))
    mask_arr <- array(FALSE, dims)
    mask_arr[c(1, 2, 5, 7, 8, 11, 13, 17, 20, 24)] <- TRUE
    mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dim = dims))
    dset <- neurovec_frame(scans = list(vec), mask = mask, TR = 2, event_table = etab)
  } else {
    mask_arr <- NULL
    dset <- matrix_frame(matrix(rnorm(n * 5), n, 5), TR = 2, run_length = n,
                         event_table = etab)
  }
  fit <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset,
                 progress = FALSE)
  list(fit = fit, mask_arr = mask_arr)
}

# Independent reference: estimate / statistic columns from the stored result,
# placed into the mask by hand.
.p2_expected_image <- function(fit, mask_arr, j, element) {
  mat <- fit$result$betas$data[[1]][[element]][[1]]
  ei <- fit$result$event_indices
  vals <- mat[, ei[j]]
  if (is.null(mask_arr)) return(unname(vals))
  arr <- array(NA_real_, dim(mask_arr))
  arr[which(mask_arr)] <- vals
  arr
}

.p2_img_array <- function(img) {
  array(as.vector(as.array(img)), dim(img))
}

test_that("P2-2: coef_images() equals per-coefficient reconstruction", {
  fx <- .p2_coef_fixture(spatial = TRUE)
  nms <- coef_names(fx$fit, type = "estimates")
  expect_length(nms, 3L)
  for (st in c("estimate", "se", "tstat")) {
    element <- c(estimate = "estimate", se = "se", tstat = "stat")[[st]]
    imgs <- coef_images(fx$fit, statistic = st, type = "estimates")
    expect_identical(names(imgs), nms)
    for (j in seq_along(nms)) {
      expect_s4_class(imgs[[j]], "NeuroVol")
      expect_identical(.p2_img_array(imgs[[j]]),
                       .p2_expected_image(fx$fit, fx$mask_arr, j, element))
      expect_identical(dim(neuroim2::space(imgs[[j]])), dim(fx$mask_arr))
      # Same object as the single-coefficient accessor.
      expect_identical(imgs[[j]], coef_image(fx$fit, coef = nms[j], statistic = st))
    }
  }
  sub <- coef_images(fx$fit, coefs = nms[c(3, 1)])
  expect_identical(names(sub), nms[c(3, 1)])
  expect_identical(.p2_img_array(sub[[1]]), .p2_expected_image(fx$fit, fx$mask_arr, 3, "estimate"))
})

test_that("P2-2: coef_images() on a non-spatial dataset returns raw vectors", {
  fx <- .p2_coef_fixture(spatial = FALSE)
  imgs <- coef_images(fx$fit)
  nms <- coef_names(fx$fit, type = "estimates")
  expect_identical(names(imgs), nms)
  for (j in seq_along(nms)) {
    expect_identical(imgs[[j]], .p2_expected_image(fx$fit, NULL, j, "estimate"))
  }
})

test_that("P2-2 guard: names and mask space are resolved once per coef_images() call", {
  fx <- .p2_coef_fixture(spatial = TRUE)
  n_names <- 0L
  n_space <- 0L
  real_names <- get("coef_names", envir = asNamespace("fmrireg"))
  real_space <- get(".fmri_dataset_mask_space", envir = asNamespace("fmrireg"))
  testthat::local_mocked_bindings(
    coef_names = function(...) {
      n_names <<- n_names + 1L
      real_names(...)
    },
    .fmri_dataset_mask_space = function(...) {
      n_space <<- n_space + 1L
      real_space(...)
    }
  )
  imgs <- coef_images(fx$fit)
  expect_length(imgs, 3L)
  expect_identical(n_names, 1L)
  expect_identical(n_space, 1L)
})
