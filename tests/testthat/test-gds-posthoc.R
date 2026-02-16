test_that("posthoc fdr:spatial matches spatial_fdr() per coefficient", {
  skip("GDS spatial FDR parity under review")
  skip_if_not_installed("fmrigds")
  set.seed(42)
  # Build a tiny gds with z assay (features x coef)
  P <- 50; K <- 2
  z_mat <- matrix(rnorm(P*K), nrow = P, ncol = K)
  # as_gds expects 3D arrays: [sample x subject x contrast]
  z_arr <- array(z_mat, dim = c(P, 1, K))
  # Two groups of equal size
  grp <- rep(1:5, each = 10)
  g <- fmrigds::as_gds(list(z = z_arr), space = NULL)
  # Apply spatial FDR via posthoc
  pl <- fmrigds::as_plan(g)
  pl <- fmrigds::posthoc(pl, method = "fdr:spatial", options = list(group = grp))
  gq <- fmrigds::compute(pl)
  q_assay <- fmrigds::assay(gq, "q")
  if (is.null(dim(q_assay))) {
    q_assay <- matrix(q_assay, ncol = 1)
  } else if (length(dim(q_assay)) == 3 && dim(q_assay)[2] == 1) {
    q_assay <- q_assay[, 1, , drop = TRUE]
    if (is.null(dim(q_assay))) q_assay <- matrix(q_assay, ncol = 1)
  }
  # Compare column-wise
  for (k in 1:K) {
    sres <- spatial_fdr(z = z_mat[,k], group = grp)
    expect_equal(as.numeric(q_assay[,k]), as.numeric(sres$q))
  }
})
