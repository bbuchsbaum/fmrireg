library(devtools)
library(neuroim2)

library(dplyr)
library(broom)
library(stringr)
load_all()

mask <- read_vol("test_data/phoneme/sub-1006/automask.nii")
mask.idx <- which(mask>0)
events <- list.files("test_data/phoneme/sub-1006/", ".*tsv", full.names=TRUE)
des <- do.call(rbind, lapply(events, read.table, header=TRUE))
des <- as_tibble(des %>% mutate(Onset = Onset/1000))
des$constant <- factor(rep("cons", nrow(des)))
sframe <- sampling_frame(rep(244, length(events)), TR=1.5)
evmod <- event_model(Onset ~ hrf(constant, basis="spmg3"), block = ~ block_num, sampling_frame=sframe, data=des)
bmod <- baseline_model(basis="bs", degree=3, sframe=sframe)
fmod <- fmri_model(evmod, bmod)

dmat <- design_matrix(fmod)
svec <- read_vec("test_data/phoneme/sub-1006/sub-1006_task-phoneme_bold_space-MNI152NLin2009cAsym_latent.lv.h5")
basis <- svec@basis

Xmat <- as.matrix(dmat)
ret <- lapply(1:ncol(basis), function(i) {
  y <- basis[,i]
  lm.1 <- lm(y ~ Xmat)
  tidy(lm.1) %>% filter(str_detect(term, "Xmatconstant:basis")) %>% mutate(component=i, basis=1:3)
}) %>% bind_rows()

ret1 <- ret %>% filter(basis==1)
bvol1 <- ret1$estimate %*% t(svec@loadings)
bvol1 <- NeuroVol(bvol1[1,], space(mask), indices=mask.idx)

ret2 <- ret %>% filter(basis==2)
bvol2 <- ret2$estimate %*% t(svec@loadings)
bvol2 <- NeuroVol(bvol2[1,], space(mask), indices=mask.idx)

ret3 <- ret %>% filter(basis==3)
bvol3 <- ret3$estimate %*% t(svec@loadings)
bvol3 <- NeuroVol(bvol3[1,], space(mask), indices=mask.idx)

dset <- matrix_dataset(as.matrix(basis), TR=1.5, run_length=rep(244, length(events)),event_table=des)
con1  <- pair_contrast(~ voicing == "Vls", ~ voicing == "Vcd", name="Vls_min_Vcd")
con2  <- pair_contrast(~ manner == "Fric", ~ manner == "Stop", name="Fric_min_Stop")
con3  <- pair_contrast(~ place == "Lab", ~ place != "Lab", name="Lab_min_Other")
cset <- contrast_set(con1,con2)
lm.1 <- fmri_lm(Onset ~ hrf(voicing, manner, contrasts=con1, subset=voicing!="none"), block= ~ block_num, 
                dataset=dset)
lm.2 <- fmri_lm(Onset ~ hrf(manner, contrasts=con2), block= ~ block_num, dataset=dset)
lm.3 <- fmri_lm(Onset ~ hrf(place, contrasts=con3), block= ~ block_num, dataset=dset)

cmap = RColorBrewer::brewer.pal(11, "Spectral")
est1 <- lm.1$result$contrasts$estimate()
tvcd <- t(est1) %*% t(svec@loadings)
tvcdvol <- NeuroVol(tvcd[1,], space(mask), indices=mask.idx)
plot(tvcdvol, zlevel=seq(16,40,by=3), cmap=rev(cmap), thresh=c(-.03, .03))

logp= -log(lm.2$result$contrasts$prob())
est1 <- lm.2$result$contrasts$estimate()
#stat= lm.2$result$contrasts$stat()
tmanner <- t(est1) %*% t(svec@loadings)
tmannervol <- NeuroVol(tmanner[1,], space(mask), indices=mask.idx)
plot(tmannervol, zlevel=seq(16,40,by=3), cmap=rev(cmap), thresh=c(-.03, .03))

se1 <- lm.3$result$contrasts$se()
est1 <- lm.3$result$contrasts$estimate()
#stat= lm.2$result$contrasts$stat()
place_est <- t(est1) %*% t(svec@loadings)
place_se <-  t(se1) %*% t(svec@loadings)
tplacevol <- NeuroVol(place_est[1,], space(mask), indices=mask.idx)
plot(tplacevol, zlevel=seq(16,40,by=3), cmap=rev(cmap), thresh=c(-.03, .03))


ldset <- latent_dataset(svec, TR=1.5, run_length=rep(244, length(events)),event_table=des)

con1  <- pair_contrast(~ syllable %in% c("ma", "na"), ~ !(syllable %in% c("ma", "na")), name="nasals")
res <- try(fmri_latent_lm(Onset ~ hrf(syllable, contrasts=con1, subset=syllable!="null"), 
                      block= ~ block_num, dataset=ldset))


