
ridge_betas <- function(X, Y, penalty_factor=rep(1:ncol(X)), lambda=.01) {
  fit <- glmnet(X, Y, penalty.factor=penalty_factor, alpha=0,lambda=lambda)
  coef(fit)[,1,drop=FALSE]
}

pls_betas <- function(X, Y, ncomp=3) {
  dx <- data.frame(X=X, Y=Y)
  fit <- pls::plsr(Y ~ X, data=dx, ncomp=ncomp, method="simpls", scale=TRUE)
  coef(fit, ncomp=ncomp)[,,1]
}

pls_global_betas <- function(X, Y, ncomp=3) {
  dx <- data.frame(X=X, Y=Y)
  fit <- pls::plsr(Y ~ X, data=dx, ncomp=ncomp, method="widekernelpls", scale=TRUE)
  coef(fit, ncomp=ncomp)[,,1]
}

slm_betas <- function(X, Y) {
  slm.1 <- care::slm(X, Y, verbose=FALSE)
  b2 <- coef(slm.1)[,-(1:2)]
  b1 <- coef(slm.1)[,1]
  b1 + b2
}

#' @importFrom rrBLUP mixed.solve
mixed_betas <- function(X, Y, ran_ind, fixed_ind) {
  fit <- rrBLUP::mixed.solve(Y, Z=X[,ran_ind], X=X[,c(fixed_ind)], bounds=c(c(1e-05, .2)))
  c(fit$u, fit$b)
}

gen_beta_design <- function(fixed, ran, block, bmod, dset) {
  
  emod_fixed <- event_model(fixed, data=dset$event_table, block=block, sampling_frame=dset$sampling_frame)
  emod_ran <- event_model(ran, data=dset$event_table, block=block, sampling_frame=dset$sampling_frame)
  
  dmat_fixed <- design_matrix(emod_fixed)
  dmat_ran <- design_matrix(emod_ran)
  
  dmat_base <- design_matrix(bmod)
  #dmat_all <- cbind(scale(dmat_ran), scale(dmat_fixed), scale(dmat_base))
  dmat_all <- cbind(dmat_ran, dmat_fixed, dmat_base)
  
  
  start_ran <- 1
  start_fixed <- ncol(dmat_ran)+1
  start_base <- start_fixed + ncol(dmat_fixed)
  ran_ind <- 1:ncol(dmat_ran)
  
  fixed_ind <- start_fixed:(start_base-1)
  base_ind <- start_base:(ncol(dmat_all))
  
  
  list(bmod=bmod,
       emod_fixed=emod_fixed,
       emod_ran=emod_ran,
       dmat_fixed=dmat_fixed,
       dmat_ran=dmat_ran,
       dmat_base=dmat_base,
       ran_ind=ran_ind,
       fixed_ind=fixed_ind,
       base_ind=base_ind)
}
  

run_estimate_betas <- function(bdes, dset, method, ncomp=3) {

  betas <- with(bdes, {
    vecs <- vectors(get_data(dset), subset = which(get_mask(dset) >0))
    if (method == "slm") {
      ## TODO FIXME
      X  <- cbind(dmat_ran, dmat_fixed)
      Base <- as.matrix(dmat_base)
      
      res <- do.call(cbind, furrr::future_map(vecs, function(v) { 
        v0 <- resid(lsfit(Base, v, intercept = FALSE))
        slm_betas(X, v0) 
      }))
      as.matrix(res)
      ##neuroim2::NeuroVec(as.matrix(res), neuroim2::add_dim(neuroim2::space(mask),nrow(res)))
      
    } else if (method == "mixed") {
      X  <- cbind(dmat_ran, dmat_fixed)
      X[is.na(X)] <- 0
      Base <- as.matrix(dmat_base)
      res <-
        do.call(cbind, furrr::future_map(vecs, function(v) {
          v0 <- resid(lsfit(Base, v, intercept = FALSE))
          mixed_betas(X, v0, ran_ind = ran_ind, fixed_ind = fixed_ind)
        }))
      as.matrix(res)
    } else if (method == "pls") {
      X <- cbind(scale(dmat_ran), scale(dmat_fixed))
      X[is.na(X)] <- 0
      
      Base <- as.matrix(dmat_base)
      
      res <-
        do.call(cbind, furrr::future_map(vecs, function(v) {
          v0 <- resid(lsfit(Base, v, intercept = FALSE))
          pls_betas(X, v0, ncomp = ncomp)
        }))
      as.matrix(res)
    } else if (method == "pls_global") {
  
      X <- cbind(scale(dmat_ran), scale(dmat_fixed))
      X[is.na(X)] <- 0
      
      Base <- as.matrix(dmat_base)
      Y <- do.call(cbind, lapply(vecs, function(v) v))
      
      if (ncomp < log(ncol(Y))) {
        warning("'ncomp' for pls_global method is less than log(nvoxels), consider increasing.")
      }
      
      Y0 <- resid(lsfit(Base, Y, intercept = FALSE))
      pls_global_betas(X, Y0, ncomp=ncomp)
    } else {
      X <- cbind(scale(dmat_ran), scale(dmat_fixed))
      Base <- as.matrix(dmat_base)
      
      res <- Reduce("+", furrr::future_map(1:niter, function(iter) {
        slight <- neuroim2::random_searchlight(mask, radius = radius)
        mset <- Reduce("+", furrr::future_map(slight, function(s) {
          cds <- neuroim2::coords(s)
          Y <- neuroim2::series(bvec, cds)
          Y0 <- resid(lsfit(Base, Y, intercept = FALSE))
          
          dx <- list(Y = Y0, X = X)
          res <-
            pls::plsr(
              Y ~ X,
              data = dx,
              ncomp = ncomp,
              validation = "none",
              method = "simpls"
            )
          idx <- neuroim2::grid_to_index(mask, cds)
          B <- coef(res)[, , 1, drop = FALSE]
          m <- Matrix::sparseMatrix(
            i = rep(1:nrow(B), length(idx)),
            j = rep(idx, each = nrow(B)),
            x = as.numeric(B),
            dims = c(nrow(B), prod(dim(mask)))
          )
          #com <- colMeans(cds)
          #D <- sqrt(rowSums(sweep(cds, 2, com, "-")^2))
          m
          
        }))
      })) / niter
      
      as.matrix(res)[, mask != 0]
    }
  })
  
  
}


#' estimate trialwise betas for an fMRI dataset.
#' 
#' @param dataset the fmri_dataset instance
#' @param fixed the fixed factors
#' @param ran the random (trialwise) factors
#' @param block the block factor
#' @param basemod \code{baseline_model} instance to regress out of data before bet estimation
#' @param niter number of searchlight iterations for method "pls_searchlight"
#' @param ncomp number of pls components for method "pls" and "pls_searchlight" and "pls_global"
#' @param lambda lambda parameter (not currently used)
#' @import pls
#' @importFrom care slm
#' @examples 
#' 
#' facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
#' facesdes$frun <- factor(facedes$run)
#' scans <- paste0("rscan0", 1:6, ".nii")
#' 
#' dset <- fmri_dataset(scans=scans,mask="mask.nii",TR=1.5,run_length=rep(436,6),event_table=facedes)
#' fixed = onset ~ hrf(run)
#' ran = onset ~ trialwise()
#' block = ~ run
#' 
#' 
## TODO trialwise(durations=4) failed
estimate_betas.fmri_dataset <- function(dataset,fixed, ran, block,  
                           method=c("mixed", "pls", "pls_searchlight", "pls_global"), 
                           basemod=NULL, 
                           radius=8, niter=8, ncomp=4, lambda=.01) {
  
  dset <- dataset
  bvec <- get_data(dset)
  mask <- get_mask(dset)
  
  bmod <- if (is.null(basemod)) {
    baseline_model("constant", sframe=dset$sampling_frame)
  } else {
    basemod
  }
  
  
  bdes <- gen_beta_design(fixed, ran, block, bmod, dset)
  betas <- run_estimate_betas(bdes, dset, method, ncomp=ncomp)
  
  nbetas <- nrow(betas)
  ospace_ran <- neuroim2::add_dim(neuroim2::space(mask), length(bdes$ran_ind))
  ospace_fixed <- neuroim2::add_dim(neuroim2::space(mask), length(bdes$fixed_ind))
  
  ran <- neuroim2::NeuroVec(as.matrix(betas[bdes$ran_ind,,drop=FALSE]), ospace_ran, mask=mask)
  fixed <- neuroim2::NeuroVec(as.matrix(betas[bdes$fixed_ind,,drop=FALSE]), ospace_fixed, mask=mask)
  #baseline <- neuroim2::NeuroVec(as.matrix(betas[base_ind,,drop=FALSE]), ospace, mask=mask)
  
  ret <- list(betas_fixed=fixed,
       betas_ran=ran,
       #design=bdes$dmat_all,
       design_ran=bdes$dmat_ran,
       design_fixed=bdes$dmat_fixed,
       design_base=bdes$dmat_base)
  class(ret) <- "fmri_betas"
  ret
  
}

#' @export 
estimate_betas.matrix_dataset <- function(dataset,fixed, ran, block,  
                                        method=c("mixed", "pls", "pls_searchlight"), 
                                        basemod=NULL,
                                        radius=8, niter=8, ncomp=4, lambda=.01) {
  
  
  dset <- dataset
  mask <- get_mask(dset)
  
  bmod <- if (is.null(basemod)) {
    baseline_model("constant", sframe=dset$sampling_frame)
  } else {
    basemod
  }

  bdes <- gen_beta_design(fixed, ran, block, bmod, dset)
  betas <- run_estimate_betas(bdes, dset, method, ncomp=ncomp)
  
  ran <- as.matrix(betas[bdes$ran_ind,,drop=FALSE])
  fixed <- as.matrix(betas[bdes$fixed_ind,,drop=FALSE])
  
  ret <- list(betas_fixed=fixed,
              betas_ran=ran,
              #design=bdes$dmat_all,
              design_ran=bdes$dmat_ran,
              design_fixed=bdes$dmat_fixed,
              design_base=bdes$dmat_base)
  
  class(ret) <-  c("fmri_betas")
  ret
  
}


estimate_betas.fmri_latent_dataset <- function(dataset,fixed, ran, block,  
                                        method=c("mixed", "pls"), 
                                        ncomp=4, lambda=.01) {
  

}


#' @importFrom mgcv gam s
estimate_hrf <- function(form, fixed=NULL, block, dataset, 
                           bs=c("tp", "ts", "cr", "ps"), 
                           rsam=seq(0,20,by=1),
                           basedeg=5, nuisance_list=NULL) {
  dset <- dataset
  bvec <- get_data(dset)
  mask <- get_mask(dset)
  
  onset_var <- lazyeval::f_lhs(form)
  dvars <- lazyeval::f_rhs(form)
  
  
  bmod <- baseline_model("bs", degree=basedeg, sframe=dset$sampling_frame, nuisance_list=nuisance_list)
  
  if (!is.null(fixed)) {
    emod_fixed <- event_model(fixed, data=dset$event_table, block=block, sampling_frame=dset$sampling_frame)
    X_fixed <- as.matrix(design_matrix(emat_fixed))
    has_fixed=TRUE
  } else {
    has_fixed=FALSE
  }
  
  emat_cond <- event_model(form, data=dset$event_table, block=block, 
                           sampling_frame=dset$sampling_frame)
  

  X_base <- as.matrix(design_matrix(bmod))
  X_cond <- as.matrix(design_matrix(emat_cond))
  #browser()
  res <- do.call(cbind, furrr::future_map(neuroim2::vectors(bvec, subset=which(mask>0)), function(v) {
    gam.1 <- if (has_fixed) {
      gam(v ~ s(X_cond, bs=bs) + X_fixed + X_base)
    } else {
      gam(v ~ s(X_cond, bs=bs) + X_base)
    }
    
    time <- rsam
    xb <- matrix(colMeans(X_base), length(time),ncol(X_base), byrow=TRUE)
    predict(gam.1, list(X_cond=time, X_base=xb))
  }))
  
  res
  
}
  
  

  
  
  
                             

