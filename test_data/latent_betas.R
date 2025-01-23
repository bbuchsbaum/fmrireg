library(neuroim2)
library(dplyr)

mask <- read_vol("~/Dropbox/analysis/lowrank_marie/sub-1008/func/automask.nii")
mask.idx <- which(mask>0)
events <- list.files("~/Dropbox/analysis/lowrank_marie/sub-1008/func", ".*tsv", full.names=TRUE)
des <- do.call(rbind, lapply(events, read.table, header=TRUE))
des <- as_tibble(des %>% mutate(Onset = Onsettime))
des$constant <- factor(rep("cons", nrow(des)))
sframe <- sampling_frame(rep(320, length(events)), TR=1.5)
lvec <- read_vec("~/Dropbox/analysis/lowrank_marie/sub-1008/func/sub-1008_task-encode_bold_space-MNI152NLin2009cAsym_latent.lv.h5")
ldset <- latent_dataset(lvec, TR=1.5, run_length=rep(320, length(events)),
                        event_table=des)

ret <- fmri_latent_lm(Onsettime ~ hrf(Video), block=~run, dataset=ldset, autocor="arma")


bdes <- gen_beta_design(Onset ~ hrf(constant), Onset ~ trialwise(), block= ~ run, 5, ldset)

X <- as.matrix(cbind(bdes$dmat_ran,bdes$dmat_fixed))
Y <- as.matrix(ldset$datamat)
dx <- data.frame(X=X, Y=Y)
ncomp=78
fit <- pls::plsr(Y ~ X, data=dx, ncomp=ncomp, method="simpls", scale=TRUE)
betas <- coef(fit, ncomp=ncomp)[,,1]
recon = betas %*% t(ldset$lvec@loadings)

options('future.globals.maxSize' = 4014*1024^2)

rmat <- lvec@basis %*% t(lvec@loadings)
#svec <- SparseNeuroVec(as.matrix(rmat), space=add_dim(space(lvec@mask), 2240), mask=lvec@mask)
dset <- matrix_dataset(as.matrix(rmat[,1:20]), TR=1.5, run_length=rep(320,7), event_table=ldset$event_table )  
  
recon2 <- estimate_betas(dset, fixed = Onset ~ hrf(constant), ran = Onset ~ trialwise(), block = ~ run, 
                       method="pls", ncomp=8)


