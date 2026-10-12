# write_results(format = "h5"): hdf5r's close_all() (which always runs a full
# gc()) vs .h5_close() (closes the file directly when nothing else is open).
# With fmristore before its handle-closing fix, .h5_close() falls back to
# close_all() and the two rows are the same.
devtools::load_all("/home/user/fmrireg", quiet = TRUE)
set.seed(1)
dims <- c(40, 40, 20); nt <- 120
scans <- lapply(1:2, function(r) neuroim2::NeuroVec(array(rnorm(prod(dims) * nt), c(dims, nt)),
                                                    neuroim2::NeuroSpace(c(dims, nt))))
mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dims))
ev <- data.frame(onset = rep(seq(5, 105, by = 10), 2),
                 condition = factor(rep(c("A", "B", "C"), length.out = 22)),
                 run = rep(1:2, each = 11))
dset <- neurovec_frame(scans = scans, mask = mask, TR = 1, event_table = ev)
con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B", name = "AvB"),
                    pair_contrast(~ condition == "A", ~ condition == "C", name = "AvC"),
                    pair_contrast(~ condition == "B", ~ condition == "C", name = "BvC"))
fit <- fmri_lm(onset ~ hrf(condition, contrasts = con), block = ~ run, dataset = dset,
               progress = FALSE)
junk <- replicate(200, rnorm(5e4), simplify = FALSE)   # a realistic R heap for gc() to walk
time_write <- function(close_fun) {
  ns <- asNamespace("fmrireg")
  old <- get(".h5_close", ns); unlockBinding(".h5_close", ns)
  assign(".h5_close", close_fun, ns); on.exit({assign(".h5_close", old, ns); lockBinding(".h5_close", ns)})
  out <- tempfile(); dir.create(out)
  t <- system.time(write_results(fit, path = out, subject = "01", task = "t",
                                 space = "MNI152NLin2009cAsym", format = "h5",
                                 strategy = "by_contrast"))[["elapsed"]]
  list(t = t, files = sort(list.files(out, recursive = TRUE)), dir = out)
}
cat("fmristore", format(packageVersion("fmristore")),
    "- open ids after a write are closed by fmristore:",
    { h <- fmristore::write_labeled_vec(scans[[1]][[1]] * 0 + 1 + neuroim2::NeuroVec(array(0, c(dims, 1)), neuroim2::NeuroSpace(c(dims, 1))), mask, "a", file = tempfile(fileext = ".h5"))
      n <- length(h$get_obj_ids()); h$close_all(); n == 1L }, "\n")
for (rep in 1:3) {
  old <- time_write(function(h5) h5$close_all())
  new <- time_write(fmrireg:::.h5_close)
  stopifnot(identical(old$files, new$files))
  same <- all(vapply(old$files[grepl("h5$", old$files)], function(f) {
    a <- fmristore::read_labeled_vec(file.path(old$dir, f))
    b <- fmristore::read_labeled_vec(file.path(new$dir, f))
    on.exit({try(close(a), silent = TRUE); try(close(b), silent = TRUE)})
    identical(a@labels, b@labels) &&
      all(vapply(a@labels, function(l) identical(as.numeric(a[[l]][mask]),
                                                 as.numeric(b[[l]][mask])), logical(1)))
  }, logical(1)))
  cat(sprintf("rep %d: close_all %.2f s, .h5_close %.2f s, same files/contents: %s\n",
              rep, old$t, new$t, same))
}
