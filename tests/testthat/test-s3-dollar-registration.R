# Issue C: the `$` methods for lazy chunks and pre-projected designs must be
# registered S3 methods (NAMESPACE `S3method("$", ...)`), not merely functions
# that happen to be visible from inside the namespace. Code outside fmrireg
# (future workers, user code, other packages) only reaches them through the
# registration, so these tests dispatch from an environment that cannot see
# the fmrireg namespace.

.outside_dollar <- function(name) {
  # function(obj) obj$<name>, defined in (and dispatching from) baseenv()
  eval(call("function", as.pairlist(alist(obj = )),
            call("$", quote(obj), as.name(name))),
       baseenv())
}

test_that("`$` methods for data_chunk and fmri_preproject are registered", {
  for (cls in c("data_chunk", "fmri_preproject")) {
    m <- utils::getS3method("$", cls, optional = TRUE, envir = baseenv())
    expect_true(is.function(m), info = cls)
    expect_identical(m, get(paste0("$.", cls), envir = asNamespace("fmrireg")),
                     info = cls)
  }
})

test_that("`chunk$data` reads a lazy chunk when called from outside fmrireg", {
  set.seed(1)
  Y <- matrix(rnorm(40 * 6), 40, 6)
  ev <- data.frame(onset = c(4, 24, 44, 64), condition = factor("a"), run = 1)
  dset <- matrix_frame(Y, TR = 2, run_length = 40, event_table = ev)
  chunks <- fmrireg:::.dset_chunks(dset, nchunks = 2)
  ch <- chunks[[2]]
  expect_null(.subset2(ch, "data"))
  got <- .outside_dollar("data")(ch)
  expect_false(is.null(got))
  expect_equal(unname(as.matrix(got)), unname(Y[, ch$voxel_ind]))
})

test_that("`proj$Pinv` is computed when called from outside fmrireg", {
  set.seed(2)
  n <- 30
  X <- cbind(1, matrix(rnorm(n * 3), n, 3))
  proj <- fmrireg:::.fast_preproject(X)
  pinv <- .outside_dollar("Pinv")(proj)
  expect_false(is.null(pinv))
  expect_equal(c(pinv), c(MASS::ginv(X)), tolerance = 1e-10)
})
