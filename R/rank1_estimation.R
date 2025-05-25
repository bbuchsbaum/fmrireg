##############################################################################
## fast‑rank1.R      – vectorised, QR‑cached rank‑1 GLM solvers             ##
##############################################################################
#' @include utils-internal.R
NULL

############################################################
## helpers shared by both public functions                ##
############################################################
.r1_postproc <- function(h, beta, hrf_basis, hrf_ref, flip_sign) {
  hrf_est <- drop(hrf_basis %*% h)
  sc      <- max(abs(hrf_est)); if (sc < 1e-12) sc <- 1
  h       <- h   / sc
  beta    <- beta* sc
  hrf_est <- hrf_est / sc

  if (flip_sign && sum(hrf_est * hrf_ref) < 0) {
    hrf_est <- -hrf_est; beta <- -beta; h <- -h
  }
  list(h = h, beta = beta, hrf_est = hrf_est)
}

############################################################
##  r1_glm_betas()  – stacked design                      ##
############################################################
#' Rank‑1 GLM (stacked design)
#' @inheritParams r1_glm_betas
#' @param method  `"als"` (default) or `"svd"`
#' @export
r1_glm_betas <- function(X, y, Z = NULL,
                         hrf_basis,
                         hrf_ref,
                         maxit               = 100,
                         flip_sign           = TRUE,
                         use_box_constraints = FALSE,
                         method              = c("als", "svd"))
{
  method <- match.arg(method)
  n  <- length(y)
  m  <- ncol(hrf_basis)
  k  <- ncol(X) / m;  stopifnot(k == round(k))
  q  <- if (is.null(Z)) 0L else ncol(Z)
  if (is.null(Z)) Z <- matrix(0, n, 0)

  if (method == "svd" && use_box_constraints)
    stop("`method=\"svd\"` cannot enforce box‑constraints; set `use_box_constraints = FALSE` or use `method=\"als\"`.")

  ##########################################################################
  ## FAST   engine 1: SVD  (+ two OLS)                                    ##
  ##########################################################################
  if (method == "svd") {
    # Step‑1: dominant singular vector of X'y  (reshaped m×k)
    Vh     <- matrix(crossprod(X, y), m, k)
    sv     <- svd(Vh, nu = 1, nv = 1)
    h_raw  <- sv$u
    beta   <- sv$d[1] * sv$v[, 1]

    # Step‑2: optional nuisance regressors
    if (q) {
      res <- y - X %*% as.vector(kronecker(beta, h_raw))
      omega <- qr.coef(qr(Z), res)
    } else omega <- numeric(0)

    out  <- .r1_postproc(h_raw, beta, hrf_basis, hrf_ref, flip_sign)
    return(list(beta      = out$beta,
                h         = out$h,
                omega     = omega,
                converged = TRUE,
                value     = NA_real_))
  }

  ##########################################################################
  ## FAST   engine 2: ALS  (+ L‑BFGS‑B polish)                            ##
  ##########################################################################
  ## ---- split X once so we can re‑use BLAS                               ##
  Xb   <- lapply(seq_len(k), \(j) X[, ((j-1)*m + 1):(j*m), drop = FALSE])
  Xt   <- lapply(Xb, t)
  Xall <- X %*% kronecker(matrix(1, k, 1), diag(m))
  Xtall<- t(Xall)

  ## ---- ALS warm start  --------------------------------------------------
  h    <- colMeans(hrf_basis)          # non‑zero start
  beta <- rep(0, k);  z <- rep(0, k)

  for (iter in 1:3) {
    # 1) beta,z  given h
    D <- cbind(do.call(cbind, lapply(Xb, \(B) B %*% h)),
               Xall %*% h,
               Z)
    coef <- qr.coef(qr(D), y)
    beta <- coef[        1:k]
    z    <- rep(coef[k+1], k)
    omega<- if (q) coef[-(1:(k+1))] else numeric(0)

    # 2)    h    given beta,z
    y2   <- y - Z %*% omega
    pred <- rowSums(mapply(\(B,b) (B %*% h) * b, Xb, beta, SIMPLIFY = FALSE))
    h    <- qr.coef(qr(Xall * sum(z)), y2 - pred)
  }

  ## ---- L‑BFGS‑B  polish  ------------------------------------------------
  par0   <- c(h, beta, z, omega)
  lower  <- upper <- rep(-Inf, length(par0))
  if (use_box_constraints && m == 3) {
    lower[m+1] <- upper[m+1] <- 1
    lower[m+2:3] <- -1; upper[m+2:3] <- 1
  }

  htemp  <- matrix(0, n, k)          # reused workspace

  obj_grad <- function(w) {
    h    <- w[1:m]
    beta <- w[m + 1L + 0:(k-1)]
    z    <- w[m + k + 1L + 0:(k-1)]
    omega<- if (q) w[-(1:(m+2*k))] else numeric(0)

    for (j in seq_len(k))
      htemp[, j] <- Xb[[j]] %*% h

    pred <- htemp %*% beta + (Xall %*% h) * sum(z) +
            if (q) Z %*% omega else 0
    res  <- y - pred

    f    <- 0.5 * sum(res^2) - 1e-6 * h[1]^2          # tiny ridge

    # gradients (all BLAS)
    grad_beta <- -colSums(htemp * res)
    g_z       <- -sum((Xall %*% h) * res)
    grad_z    <- rep(g_z, k)
    grad_h    <- -Reduce(`+`, Map(\(B,b) b * (t(B) %*% res), Xt, beta)) -
                 sum(z) * (Xtall %*% res) + h
    grad_omega<- if (q) -crossprod(Z, res) else numeric(0)

    list(value    = f,
         gradient = c(grad_h, grad_beta, grad_z, grad_omega))
  }

  opt <- optim(par0,
               fn = \(p) obj_grad(p)$value,
               gr = \(p) obj_grad(p)$gradient,
               method  = "L-BFGS-B",
               lower   = lower,
               upper   = upper,
               control = list(maxit = maxit, factr = 1e7))

  h      <- opt$par[1:m]
  beta   <- opt$par[m + 1L + 0:(k-1)]
  omega  <- if (q) opt$par[-(1:(m+2*k))] else numeric(0)

  out <- .r1_postproc(h, beta, hrf_basis, hrf_ref, flip_sign)

  list(beta      = out$beta,
       h         = out$h,
       omega     = omega,
       converged = opt$convergence == 0,
       value     = opt$value)
}

############################################################
##  r1_glms_betas() – Mumford/split designs               ##
############################################################
#' Rank‑1 GLM (Mumford / split designs)
#' @inheritParams r1_glms_betas
#' @param method `"als"` (default) or `"svd"`
#' @export
r1_glms_betas <- function(X_list,
                          X_all = NULL,
                          y, Z = NULL,
                          hrf_basis,
                          hrf_ref,
                          maxit  = 100,
                          method = c("als", "svd"))
{
  method <- match.arg(method)
  n  <- length(y)
  k  <- length(X_list)
  m  <- ncol(hrf_basis)
  q  <- if (is.null(Z)) 0L else ncol(Z)
  if (is.null(Z)) Z <- matrix(0, n, 0)

  if (is.null(X_all))
    X_all <- do.call(cbind, X_list) %*%
             kronecker(matrix(1, k, 1), diag(m))

  if (method == "svd") {
    # build the stacked design once, reuse r1_glm_betas()
    X_stack <- do.call(cbind, X_list)
    fit <- r1_glm_betas(X_stack, y, Z,
                        hrf_basis, hrf_ref,
                        method = "svd")
    attr(fit$beta, "estimated_hrf")  <- fit$h %*% t(hrf_basis)
    attr(fit$beta, "basis_weights")  <- fit$h
    return(fit$beta)
  }

  ## -------------- ALS + L‑BFGS‑B  ----------------------------------------
  Xt      <- lapply(X_list, t)
  Xt_all  <- t(X_all)
  htemp   <- matrix(0, n, k)

  init_par<- c(rep(1, m), rep(1, k), rep(1, k))  # h, beta, z

  obj_grad <- function(w) {
    h    <- w[1:m]
    beta <- w[m + 1L + 0:(k-1)]
    z    <- w[m + k + 1L + 0:(k-1)]

    for (j in seq_len(k))
      htemp[, j] <- X_list[[j]] %*% h

    pred <- htemp %*% beta + (X_all %*% h) * sum(z)
    res  <- y - pred
    f    <- 0.5 * sum(res^2) - 1e-6*h[1]^2

    grad_beta <- -colSums(htemp * res)
    g_z       <- -sum((X_all %*% h) * res)
    grad_z    <- rep(g_z, k)
    grad_h    <- -Reduce(`+`, Map(\(B,b) b * (t(B) %*% res), Xt, beta)) -
                 sum(z) * (Xt_all %*% res) + h

    list(value    = f,
         gradient = c(grad_h, grad_beta, grad_z))
  }

  opt <- optim(init_par,
               fn  = \(p) obj_grad(p)$value,
               gr  = \(p) obj_grad(p)$gradient,
               method = "L-BFGS-B",
               control = list(maxit = maxit, factr = 1e7))

  h     <- opt$par[1:m]
  beta  <- opt$par[m + 1L + 0:(k-1)]
  out   <- .r1_postproc(h, beta, hrf_basis, hrf_ref, TRUE)

  attr(out$beta, "estimated_hrf") <- out$hrf_est
  attr(out$beta, "basis_weights") <- out$h
  out$beta
}