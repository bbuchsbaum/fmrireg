
#' @keywords internal
#' @noRd
setup_hrf_library <- function(hrflib = NULL, rsam = seq(0, 24, by = 1), ncomp = 5) {
  if (is.null(hrflib)) {
    params <- expand.grid(
      h1 = seq(1, 3, by = 0.33),
      h2 = seq(3, 7, by = 0.33),
      h3 = seq(6, 8, by = 1),
      h4 = seq(5, 9, by = 1),
      f1 = seq(0, 0.2, by = 0.1),
      f2 = seq(0, 0.2, by = 0.1)
    )
    hrflib <- gen_hrf_library(hrf_half_cosine, params)
  }
  
  L <- hrflib(rsam)
  pca_L <- multivarious::pca(L, preproc = multivarious::pass())
  Lk <- pca_L$u[, 1:ncomp, drop = FALSE]
  
  basis_set <- lapply(1:ncol(Lk), function(i) {
    gen_empirical_hrf(rsam, Lk[, i])
  })
  
  list(L = L, Lk = Lk, basis_set = basis_set)
}



#' @export
find_best_hrf <- function(dataset, onset_var, hrflib = NULL, rsam = seq(0, 24, by = 1),
                          ncomp = 5, nsplits = 3, block = NULL, basemod = NULL, 
                          cluster_series = TRUE, 
                          cluster_quantile = .5,
                          dist_method = "euclidean",
                          spatial_smoothing = FALSE,
                          lambda = 0.1,
                          n_iterations = 5,
                          k_neighbors = 6,
                          nboot = 1) {
  UseMethod("find_best_hrf")
}


#' @keywords internal
#' @noRd
compute_best_hrfs <- function(X_proj, L_proj, L, cluster_series = TRUE, cluster_quantile = .5) {
  basis_euc <- proxy::dist(L_proj, X_proj, method = "euclidean")
  best_match_indices <- apply(as.matrix(basis_euc), 2, which.min)
  
  if (!cluster_series) {
    return(list(
      best_hrf_indices = best_match_indices,
      L_best_hrfs = L[, best_match_indices, drop = FALSE],
      L_proj = L_proj,
      X_proj = X_proj
    ))
  }
  
  unique_best_hrfs <- unique(best_match_indices)
  L_best_hrfs <- L[, unique_best_hrfs, drop = FALSE]
  hrf_dist <- dist(t(L_best_hrfs))
  hclust_result <- hclust(hrf_dist, method = "average")
  clusters <- cutree(hclust_result, h = quantile(hrf_dist, cluster_quantile))
  
  cluster_medoid <- function(cluster_indices) {
    sub_dist <- as.matrix(hrf_dist)[cluster_indices, cluster_indices]
    medoid_index <- cluster_indices[which.min(rowSums(sub_dist))]
    return(unique_best_hrfs[medoid_index])
  }
  
  unique_clusters <- unique(clusters)
  cluster_representatives <- sapply(unique_clusters, function(cluster) {
    cluster_indices <- which(clusters == cluster)
    cluster_medoid(cluster_indices)
  })
  
  reassigned_hrfs <- sapply(best_match_indices, function(idx) {
    cluster <- clusters[which(unique_best_hrfs == idx)]
    cluster_representatives[as.character(cluster)]
  })
  
  list(
    best_hrf_indices = unique_best_hrfs,
    clusters = clusters,
    cluster_representatives = cluster_representatives,
    reassigned_hrfs = reassigned_hrfs,
    L_best_hrfs = L_best_hrfs,
    L_proj = L_proj,
    X_proj = X_proj
  )
}

#' @keywords internal
#' @noRd
smooth_hrf_assignments <- function(reassigned_hrfs, L_best_hrfs, voxel_coords, 
                                   lambda = 0.1, n_iterations = 5, k = 6) {
  n_voxels <- length(reassigned_hrfs)
  if (n_voxels == 0) return(reassigned_hrfs)
  
  # Build neighbor graph
  dists <- as.matrix(dist(voxel_coords))
  neighbor_indices <- apply(dists, 1, function(row) {
    neigh <- order(row)
    neigh <- neigh[neigh != which.min(row)] # remove self
    head(neigh, k)
  })
  
  get_hrf_shape <- function(idx) L_best_hrfs[, idx, drop = FALSE]
  
  for (iter in seq_len(n_iterations)) {
    new_assignments <- reassigned_hrfs
    for (v in seq_len(n_voxels)) {
      v_hrf_idx <- reassigned_hrfs[v]
      neigh_idx <- neighbor_indices[, v]
      neigh_hrfs <- reassigned_hrfs[neigh_idx]
      
      # Most common neighbor HRF
      hrf_freq <- sort(table(neigh_hrfs), decreasing = TRUE)
      candidate_hrf <- as.numeric(names(hrf_freq)[1])
      
      if (candidate_hrf != v_hrf_idx) {
        v_shape <- get_hrf_shape(v_hrf_idx)
        neigh_shapes <- L_best_hrfs[, neigh_hrfs, drop = FALSE]
        current_cost <- mean(apply(neigh_shapes, 2, function(x) sqrt(sum((v_shape - x)^2))))
        
        cand_shape <- get_hrf_shape(candidate_hrf)
        candidate_cost <- mean(apply(neigh_shapes, 2, function(x) sqrt(sum((cand_shape - x)^2))))
        
        if (candidate_cost < current_cost) {
          new_assignments[v] <- candidate_hrf
        }
      }
    }
    reassigned_hrfs <- new_assignments
  }
  
  reassigned_hrfs
}



#' @export
find_best_hrf.matrix_dataset <- function(dataset, onset_var, hrflib = NULL, rsam = seq(0, 24, by = 1),
                                         ncomp = 5, nsplits = 3, block = NULL, basemod = NULL, 
                                         cluster_series = TRUE, 
                                         cluster_quantile = .5,
                                         dist_method = "euclidean",
                                         spatial_smoothing = FALSE,
                                         lambda = 0.1,
                                         n_iterations = 5,
                                         k_neighbors = 6,
                                         nboot = 1) {
  
  # A helper function to run one iteration of the pipeline (no bootstrap)
  run_pipeline <- function(ds) {
    # Setup HRF library and basis functions
    hrf_setup <- setup_hrf_library(hrflib, rsam, ncomp)
    hrf_list <- do.call(gen_hrf_set, hrf_setup$basis_set)
    
    # Setup models
    if (is.null(basemod)) {
      bmod <- baseline_model("constant", sframe = ds$sampling_frame)
    } else {
      bmod <- basemod
    }
    
    # Create split variable
    split_ids <- rep(1:nsplits, length.out = nrow(ds$event_table))
    fac_var <- factor(paste0("split_", split_ids))
    ds$event_table <- ds$event_table %>% mutate(.splitvar = fac_var)
    
    # Setup event model
    form <- as.formula(paste0(onset_var, " ~ ", "hrf(.splitvar, basis = hrf_list)"))
    emod <- event_model(
      x = form,
      formula = form,
      data = ds$event_table,
      block = block,
      sampling_frame = ds$sampling_frame
    )
    
    # Get data and compute residuals
    X <- get_data_matrix(ds)
    X_base <- as.matrix(design_matrix(bmod))
    X_resid <- resid(lsfit(X_base, X, intercept = FALSE))
    
    # Dimensionality reduction
    pca_X <- multivarious::pca(X_resid, preproc = multivarious::center())
    ncomp_X <- min(ncomp, ncol(pca_X$u))
    Mk <- pca_X$u[, 1:ncomp_X, drop = FALSE]
    
    # Compute projections
    Rk <- design_matrix(emod)
    C <- t(Rk) %*% Mk
    svd_C <- svd(C)
    
    X_proj <- t(X_resid) %*% Mk %*% svd_C$v
    L_proj <- t(hrf_setup$L) %*% hrf_setup$Lk %*% svd_C$u
    
    # Compute best HRFs using chosen distance measure
    res <- compute_best_hrfs(X_proj, L_proj, hrf_setup$L, cluster_series, cluster_quantile, dist_method = dist_method)
    
    # Return result plus setup objects for reference
    list(
      res = res,
      L_best_hrfs = res$L_best_hrfs,
      dataset = ds
    )
  }
  
  # If no bootstrap, just run once
  if (nboot == 1) {
    out <- run_pipeline(dataset)
    final_assignments <- if (!is.null(out$res$reassigned_hrfs)) out$res$reassigned_hrfs else out$res$best_hrf_indices
    
    # Apply smoothing if requested
    if (spatial_smoothing && !is.null(dataset$voxel_coords)) {
      final_assignments <- smooth_hrf_assignments(
        final_assignments,
        out$L_best_hrfs,
        dataset$voxel_coords,
        lambda = lambda,
        n_iterations = n_iterations,
        k = k_neighbors
      )
      out$res$reassigned_hrfs_smoothed <- final_assignments
    }
    
    return(out$res)
  }
  
  # If nboot > 1, bootstrap
  n_voxels <- nrow(get_data_matrix(dataset))
  bootstrap_assignments <- matrix(NA, nrow = n_voxels, ncol = nboot)
  last_run <- NULL
  
  for (b in seq_len(nboot)) {
    # Resample event_table
    boot_events <- dataset$event_table[sample(seq_len(nrow(dataset$event_table)), replace = TRUE), ]
    ds_b <- dataset
    ds_b$event_table <- boot_events
    run_b <- run_pipeline(ds_b)
    
    final_assignments_b <- if (!is.null(run_b$res$reassigned_hrfs)) run_b$res$reassigned_hrfs else run_b$res$best_hrf_indices
    bootstrap_assignments[, b] <- final_assignments_b
    last_run <- run_b
  }
  
  # Majority vote for each voxel
  # For each voxel, choose the HRF that appears most frequently across bootstraps
  final_assignments <- apply(bootstrap_assignments, 1, function(x) {
    tt <- sort(table(x), decreasing = TRUE)
    as.numeric(names(tt)[1])
  })
  
  # If smoothing is requested
  if (spatial_smoothing && !is.null(dataset$voxel_coords)) {
    final_assignments <- smooth_hrf_assignments(
      final_assignments,
      last_run$L_best_hrfs,
      dataset$voxel_coords,
      lambda = lambda,
      n_iterations = n_iterations,
      k = k_neighbors
    )
  }
  
  # Construct a final result object similar to res
  # We will use the last_run$res as a template and replace HRF assignments
  final_res <- last_run$res
  # We assume final_assignments are analogous to reassigned_hrfs
  final_res$reassigned_hrfs <- final_assignments
  
  if (spatial_smoothing) {
    final_res$reassigned_hrfs_smoothed <- final_assignments
  }
  
  final_res
}

#' @export
find_best_hrf.fmri_mem_dataset <- function(dataset, onset_var, hrflib = NULL, rsam = seq(0, 24, by = 1),
                                           ncomp = 5, nsplits = 3, block = NULL, basemod = NULL, 
                                           cluster_series = TRUE, 
                                           cluster_quantile = .5,
                                           dist_method = "euclidean",
                                           spatial_smoothing = FALSE,
                                           lambda = 0.1,
                                           n_iterations = 5,
                                           k_neighbors = 6,
                                           nboot = 1) {
  find_best_hrf.matrix_dataset(
    dataset = dataset,
    onset_var = onset_var,
    hrflib = hrflib,
    rsam = rsam,
    ncomp = ncomp,
    nsplits = nsplits,
    block = block,
    basemod = basemod,
    cluster_series = cluster_series,
    cluster_quantile = cluster_quantile,
    dist_method = dist_method,
    spatial_smoothing = spatial_smoothing,
    lambda = lambda,
    n_iterations = n_iterations,
    k_neighbors = k_neighbors,
    nboot = nboot
  )
}