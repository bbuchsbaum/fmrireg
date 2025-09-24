# Legacy fMRI beta estimation helpers preserved for reference.
# These functions were removed from the main package surface when
# estimate_betas() was simplified to support only mixed, lss, and ols.

ridge_betas <- function(X, Y, penalty_factor = rep(1, ncol(X)), lambda = .01) {
  with_package("glmnet")
  fit <- glmnet::glmnet(X, Y, penalty.factor = penalty_factor, alpha = 0, lambda = lambda)
  coef(fit)[, 1, drop = FALSE]
}

pls_betas <- function(X, Y, ncomp = 3) {
  with_package("pls")
  dx <- list(X = as.matrix(X), Y = Y)
  fit <- pls::plsr(Y ~ X, data = dx, ncomp = ncomp, method = "simpls", scale = TRUE, maxit = 1500)
  coef(fit, ncomp = ncomp)[, , 1]
}

pls_global_betas <- function(X, Y, ncomp = 3) {
  with_package("pls")
  dx <- list(X = as.matrix(X), Y = Y)
  fit <- pls::plsr(Y ~ X, data = dx, ncomp = ncomp, method = "widekernelpls", scale = TRUE, maxit = 1500)
  coef(fit, ncomp = ncomp)[, , 1]
}

slm_betas <- function(X, Y) {
  with_package("care")
  slm.1 <- care::slm(X, Y, verbose = FALSE)
  b2 <- coef(slm.1)[, -(1:2)]
  b1 <- coef(slm.1)[, 1]
  b1 + b2
}

mixed_betas_cpp <- function(X, Y, ran_ind, fixed_ind) {
  X_fixed <- if (is.null(fixed_ind) || length(fixed_ind) == 0) {
    matrix(1, nrow = nrow(X), ncol = 1)
  } else {
    X[, fixed_ind, drop = FALSE]
  }

  fit <- fmrilss::mixed_solve(Y = as.matrix(Y), Z = X[, ran_ind, drop = FALSE],
                              X = X_fixed, bounds = c(1e-05, .2))
  c(fit$u, fit$b)
}

run_estimate_betas_legacy <- function(bdes, dset, method,
                                      block, maxit = 100,
                                      ncomp = 4, fracs = .5,
                                      progress = TRUE,
                                      ...) {
  method <- match.arg(method, c("lss", "lss_naive", "lss_cpp", "mixed", "mixed_cpp", "pls", "pls_global", "ols"))

  dotargs <- list(...)
  xdat <- build_design_data(bdes)

  if (method == "mixed") {
    vecs <- masked_vectors(dset)
    res <- map_voxels(vecs, function(v) {
      v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
      mixed_betas(xdat$X, v0, ran_ind = seq_len(ncol(bdes$dmat_ran)),
                  fixed_ind = if (!is.null(bdes$dmat_fixed)) {
                    (ncol(bdes$dmat_ran) + 1):(ncol(bdes$dmat_ran) + ncol(bdes$dmat_fixed))
                  } else {
                    NULL
                  })
    }, .progress = progress)
    list(beta_matrix = as.matrix(res), estimated_hrf = NULL)
  } else if (method == "mixed_cpp") {
    vecs <- masked_vectors(dset)
    res <- map_voxels(vecs, function(v) {
      v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
      mixed_betas_cpp(as.matrix(xdat$X), v0,
                      ran_ind = seq_len(ncol(bdes$dmat_ran)),
                      fixed_ind = if (!is.null(bdes$dmat_fixed)) {
                        (ncol(bdes$dmat_ran) + 1):(ncol(bdes$dmat_ran) + ncol(bdes$dmat_fixed))
                      } else {
                        NULL
                      })
    }, .progress = progress)

    list(beta_matrix = as.matrix(res), estimated_hrf = NULL)

  } else if (method == "lss_naive") {
    fmrilss::lss_naive(dset, bdes)
  } else if (method == "lss") {
    data_matrix <- get_data_matrix(dset)
    dmat_base <- as.matrix(bdes$dmat_base)
    dmat_fixed <- if (!is.null(bdes$fixed_ind)) as.matrix(bdes$dmat_fixed) else NULL
    dmat_ran <- as.matrix(bdes$dmat_ran)

    if (!is.null(dmat_fixed)) {
      nuisance_matrix <- cbind(dmat_base, dmat_fixed)
    } else {
      nuisance_matrix <- dmat_base
    }

    beta_matrix_ran <- fmrilss::lss(Y = data_matrix, X = dmat_ran, Z = NULL,
                                    Nuisance = nuisance_matrix, method = "r_optimized")

    if (!is.null(bdes$fixed_ind) && length(bdes$fixed_ind) > 0) {
      data_matrix <- get_data_matrix(dset)
      mask_idx <- which(fmridataset::get_mask(dset) > 0)
      vecs <- neuroim2::vectors(data_matrix, subset = mask_idx)
      X_base_fixed <- cbind(as.matrix(bdes$dmat_base), as.matrix(bdes$dmat_fixed))

      beta_matrix_fixed <- map_voxels(vecs, function(v) {
        fit <- lm.fit(X_base_fixed, v)
        coef(fit)[(ncol(bdes$dmat_base) + 1):length(coef(fit))]
      }, .progress = progress)

      beta_matrix <- rbind(beta_matrix_ran, beta_matrix_fixed)
    } else {
      beta_matrix <- beta_matrix_ran
    }

    list(beta_matrix = beta_matrix, estimated_hrf = NULL)
  } else if (method == "lss_cpp") {
    data_matrix <- get_data_matrix(dset)
    dmat_base <- as.matrix(bdes$dmat_base)
    dmat_fixed <- if (!is.null(bdes$fixed_ind)) as.matrix(bdes$dmat_fixed) else NULL
    dmat_ran <- as.matrix(bdes$dmat_ran)

    if (!is.null(dmat_fixed)) {
      nuisance_matrix <- cbind(dmat_base, dmat_fixed)
    } else {
      nuisance_matrix <- dmat_base
    }

    beta_matrix <- fmrilss::lss(Y = data_matrix, X = dmat_ran, Z = NULL,
                                Nuisance = nuisance_matrix, method = "cpp_optimized")
    list(beta_matrix = beta_matrix, estimated_hrf = NULL)
  } else if (method == "pls") {
    vecs <- masked_vectors(dset)
    res <- map_voxels(vecs, function(v) {
      v0 <- resid(lsfit(xdat$Base, v, intercept = FALSE))
      pls_betas(xdat$X, v0, ncomp = ncomp)
    }, .progress = progress)

    list(beta_matrix = as.matrix(res), estimated_hrf = NULL)
  } else if (method == "pls_global") {
    vecs <- masked_vectors(dset)
    Y <- map_voxels(vecs, function(v) v, .progress = progress)

    if (ncomp < log(ncol(Y))) {
      warning("'ncomp' for pls_global method is less than log(nvoxels), consider increasing.")
    }

    Y0 <- resid(lsfit(xdat$Base, Y, intercept = FALSE))
    list(beta_matrix = pls_global_betas(xdat$X, Y0, ncomp = ncomp), estimated_hrf = NULL)
  } else if (method == "ols") {
    vecs <- masked_vectors(dset)
    Y <- map_voxels(vecs, function(v) v, .progress = progress)
    Y0 <- resid(lsfit(xdat$Base, Y, intercept = FALSE))
    beta_matrix <- ols_betas(xdat$X, Y0)
    list(beta_matrix = as.matrix(beta_matrix), estimated_hrf = NULL)
  } else {
    stop("Invalid method. Supported methods are 'lss', 'lss_naive', 'mixed', 'mixed_cpp', 'pls', 'pls_global', and 'ols'")
  }
}
