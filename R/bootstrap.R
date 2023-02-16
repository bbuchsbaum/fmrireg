


multiresponse_bootstrap_lm <- function(form, data_env, 
                                       conlist, 
                                       vnames, 
                                       fcon, modmat, 
                                       block_size=30,
                                       boot_rows=FALSE,
                                       nboot=100,
                                       event_indices) {
  
  lm.orig <- if (is.null(modmat)) {
    lm(as.formula(form), data=data_env)
  } else {
    lm.fit(modmat, data_env$.y)
  }
  
  yhat <- fitted(lm.orig)
  rows <- 1:nrow(yhat)
  nblocks <- as.integer(length(rows)/block_size)
  
  first_blocks <- rep(1:nblocks, each=length(rows)/nblocks)
  
  blocks <- split(rows, rep(1:nblocks, each=length(rows)/nblocks, length.out=length(rows)))
  maxind <- max(blocks[[length(blocks)]])
  if (maxind < nrow(yhat)) {
    last_block <- (maxind+1):nrow(yhat)
    blocks <- c(blocks, list(last_block))
  }
  
  
  boot_res <- foreach(1:nboot) %dopar% {
    sam_blocks <- sample(1:length(blocks), replace=TRUE)
    samples <- unlist(blocks[sam_blocks])
    if (length(samples) > nrow(yhat)) {
      samples <- samples[1:nrow(yhat)]
    } else if (length(samples) < nrow(yhat)) {
      delta <- nrow(yhat) - length(samples) 
      samples <- c(samples, sample(rows,delta, replace=TRUE))
    }
    rstar <- lm.orig$residuals[samples,]
    ynew <- yhat + rstar
    lm.boot <- lm.fit(modmat, ynew)
    res <- fit_lm_contrasts(lm.boot, conlist, fcon, vnames, se=FALSE) 
  }
  
  ##beta_cov <- 
  con_cov <- if (length(conlist) > 0) {
    lapply(names(boot_res[[1]]$conres), function(bname) {
      bm <- do.call(rbind, lapply(boot_res, function(br) {
        #br$conres[[bname]]$estimate()
        br$conres[[bname]]$estimate
      }))
      cov(bm)
    })
  }
  
  beta_cov <- lapply(event_indices, function(i) {
    bm <- do.call(rbind, lapply(boot_res, function(br) {
      br$bstats$estimate
    }))
    cov(bm)
  })
  
  orig <- fit_lm_contrasts(lm.orig, conlist, fcon, vnames) %>% 
    purrr::list_modify(con_cov=con_cov, beta_cov=beta_cov, nboot=nboot, bootstrap=TRUE) 
  orig
}

  
  