
ridge_betas <- function(X, Y, penalty_factor=rep(1:ncol(X)), lambda=.01) {
  fit <- glmnet(X, Y, penalty.factor=penalty_factor, alpha=0,lambda=lambda)
  coef(fit)[,1]
}


#' @examples 
#' facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
#' facesdes$frun <- factor(facedes$run)
#' scans <- paste0("rscan0", 1:6, ".nii")
#' dset <- fmri_dataset(scans=scans,mask="mask.nii",TR=1.5,run_length=rep(436,6),event_table=facedes)
#' fixed = onset ~ hrf(run)
#' ran = onset ~ trialwise()
#' block = ~ run
estimate_betas <- function(dset, fixed, ran, block, method=c("ridge", "pls"), basedeg=5, nuisance_list=NULL, radius=8, niter=20, ncomp=4, lambda=.01) {
  bvec <- get_data(dset)
  mask <- get_mask(dset)
  
  
  bmod <- baseline_model("bs", degree=basedeg, sframe=dset$sampling_frame, nuisance_list=nuisance_list)
  
  emod_fixed <- event_model(fixed, data=dset$event_table, block=block, sampling_frame=dset$sampling_frame)
  emod_ran <- event_model(ran, data=dset$event_table, block=block, sampling_frame=dset$sampling_frame)
  
  dmat_fixed <- design_matrix(emod_fixed)
  dmat_ran <- design_matrix(emod_ran)
  
  dmat_base <- design_matrix(bmod)
  
  dmat_all <- cbind(dmat_ran, dmat_fixed, dmat_base)
  ran_ind <- 1:ncol(dmat_ran)
  
  betas <- if (method == "ridge") {
    X <- as.matrix(dmat_all)
    pfac <-  c(rep(1, ncol(dmat_ran)), rep(0, ncol(dmat_fixed)), rep(0, ncol(dmat_base)))
    res <- do.call(cbind, furrr::future_map(neuroim2::vectors(bvec, subset=which(mask>0)), function(v) {
      fit <- ridge_betas(X, v, penalty_factor = pfac, lambda=lambda)
      coef(fit)[,1]
    }))
    
    neuroim2::NeuroVec(as.matrix(res), neuroim2::add_dim(neuroim2::space(mask),nrow(res)))
    
  } else {
  
    res <- Reduce("+", furrr::future_map(1:niter, function(iter) {
      slight <- neuroim2::random_searchlight(mask, radius=radius)
      mset <- Reduce("+", purrr::map(slight, function(s) {
        Y <- neuroim2::series(bvec, neuroim2::coords(s))
        dx <- list(Y=Y, X=as.matrix(dmat_all))
        res <- pls::plsr(Y ~ X, data=dx, ncomp=ncomp, validation="none", method="simpls")
      #browser()
        idx = neuroim2::grid_to_index(mask, s@coords)
        B <- coef(res)[,,1,drop=FALSE]
        m=Matrix::sparseMatrix(i=rep(1:nrow(B), length(idx)), 
                           j=rep(idx, each=nrow(B)), x=as.numeric(B), 
                           dims=c(nrow(B), prod(dim(mask))))
        com <- colMeans(s@coords)
        D <- sqrt(rowSums(sweep(s@coords, 2, com, "-")^2))
        m
              
      }))
    }))/niter
  
    neuroim2::NeuroVec(as.matrix(res), neuroim2::add_dim(neuroim2::space(mask),nrow(res)))
  }
  
  list(betas=betas,
       design=dmat_all,
       design_ran=dmat_ran,
       design_fixed=dmat_fixed,
       design_base=dmat_base)
  
}