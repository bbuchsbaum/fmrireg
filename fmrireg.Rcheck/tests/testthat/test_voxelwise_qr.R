context("Voxelwise QR contrast engine")

set.seed(123)

n <- 20
p <- 3
V <- 5

X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
Y <- matrix(rnorm(n * V), n, V)

proj <- .fast_preproject(X)
fit <- .fast_lm_matrix(X, Y, proj)

Betas <- fit$betas
sigma_vec <- sqrt(fit$sigma2)
qr_obj <- proj$qr

XtXinv_list <- replicate(V, proj$XtXinv, simplify = FALSE)
qr_list <- replicate(V, qr_obj, simplify = FALSE)

conlist <- list(A = c(0, 1, 0))
attr(conlist$A, "colind") <- 2L
fconlist <- list()

res_xtx <- fit_lm_contrasts_voxelwise(Betas, sigma_vec^2, XtXinv_list,
                                       conlist, fconlist, proj$dfres)
res_qr <- fit_lm_contrasts_voxelwise_qr(Betas, qr_list, sigma_vec,
                                        conlist, fconlist, proj$dfres)

expect_equal(res_qr$A$data[[1]], res_xtx$A$data[[1]], tolerance = 1e-8)
