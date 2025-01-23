library(fmrireg)
library(dplyr)
library(mixmeta)
library(multivarious)
library(ggplot2)

atlas <- neuroatlas::get_schaefer_atlas("400", "17")
sids <- 1001:1028
hrfb <- gen_hrf(HRF_BSPLINE, N=9)
#hrfb <- readRDS("empirical_hrf.rds")
cset <- NULL

vals <- hrfb(0:22)

dfhrf <- data.frame(amplitude=c(vals[,2], vals[,1]), time=rep(0:22, 2), basis=rep(c("basis_1","basis_2"),each=23))
qplot(time, amplitude, colour=basis, data=dfhrf) + geom_line() + theme_bw(base_size=17)

elab_freq <- lapply(sids, function(sid) {
  dat <- readRDS(paste0("~/Dropbox/analysis/SOD/fmri/schaefer/", sid, "_schaefer_mean_404.rds"))
  des <- dat$design %>% mutate(elab_onset = elab_onset_ms/1000)
  des$elab <- des$elaboration
  if (sid == "1015" || sid == "1014") {
    des$elab = case_when(
      des$elab == 1 ~ 4,
      des$elab == 2 ~ 3,
      des$elab == 3 ~ 2,
      des$elab == 4 ~ 1, 
      des$elab == 0 ~ 0)
  }
  
  
  des$felab <- factor(des$elab)
  data.frame(subject=sid, elabfreq=table(des$elab), felab=names(table(des$elab)))
  
}) %>% bind_rows()

run_subject <- function(sid) {
  print(sid)
  dat <- readRDS(paste0("~/Dropbox/analysis/SOD/fmri/schaefer/", sid, "_schaefer_mean_404.rds"))
  des <- dat$design %>% mutate(elab_onset = elab_onset_ms/1000)
  des$elab <- des$elaboration
  if (sid == "1015" || sid == "1014") {
    des$elab = case_when(
                         des$elab == 1 ~ 4,
                         des$elab == 2 ~ 3,
                         des$elab == 3 ~ 2,
                         des$elab == 4 ~ 1, 
                         des$elab == 0 ~ 0)
  }
  
  
  des$felab <- factor(des$elab)
  print(table(des$elab))
  des$elab[des$elab==0] <- NA
  des$elab[is.na(des$elab)] <- mean(des$elab, na.rm=TRUE)
  des$elab <- scale(des$elab)[,1]

  dset <- fmrireg::matrix_dataset(as.matrix(dat$data), TR=1.3, run_length=rep(528,4), event_table=des)

  #con <- fmrireg::pair_contrast(~ attention == "FA", ~ attention == "DA", name="FA_min_DA")
  #con2 <- oneway_contrast(~ basis, name="basis", where=NULL)
  #con3 <- fmrireg::pair_contrast(~ basis == "basis1", ~ basis == "basis4", name="b4_min_b1")
  #cset <- contrast_set(con2,con3)
  #con_elabr <- fmrireg:::unit_contrast(~ attention:resid_elab, name="main_elabr")

  #hrfb - gen_hrf(HRF_BSPLINE, N=7)
  # ret <- fmrireg::fmri_lm(elab_onset ~ hrf(attention, felab, basis=hrfb, subset=phase=="encoding" & felab!=0), #+ 
  #                           #hrf(attention,elab, basis=hrfb, subset=phase=="encoding"), 
  #                         strategy="chunkwise",
  #                       dataset=dset, block = ~ Scan)
  
  ret <- fmrireg::fmri_lm(elab_onset ~ hrf(attention,basis=hrfb, subset=phase=="encoding" & felab!=0), #+ 
                                                     #hrf(attention,elab, basis=hrfb, subset=phase=="encoding"), 
                                                   strategy="chunkwise",
                                                 dataset=dset, block = ~ Scan)
  
  
  ret
}

res <- lapply(sids, run_subject)

get_fitted_hrfs <- function(res) {
  out <- lapply(1:length(res), function(i) {
    print(i)
    fit <- res[[i]]
    sid <- sids[i]
    tmp = fitted_hrf(fit)[[1]]
    des <- replicate(ncol(tmp$pred), tmp$design, simplify=FALSE) %>% bind_rows() %>%
      mutate(roi=rep(1:404, each=nrow(tmp$design)), subject=sid, estimate=as.vector(tmp$pred))
  }) %>% bind_rows()
}

library(purrr)
fitted <- get_fitted_hrfs(res)
mean_hrfs_attn <- dplyr::filter(fitted) %>% group_by(time, subject, attention, roi) %>% summarize(estimate=mean(estimate))
mean_hrfs_attn_elab <- dplyr::filter(fitted) %>% group_by(time, attention, felab, roi) %>% summarize(estimate=mean(estimate))
write.table(mean_hrfs_attn, "~/Dropbox/analysis/SOD/fmri/glm_bspline_mean_hrf_attn.txt", row.names=FALSE)
write.table(mean_hrfs_attn_elab, "~/Dropbox/analysis/SOD/fmri/glm_bspline_mean_hrf_attn_elab.txt", row.names=FALSE)

allhrfs <- mean_hrfs_attn %>% group_by(roi, attention) %>% dplyr::group_split() %>% map( ~ .$estimate) %>% bind_cols()

pres <- prcomp(allhrfs)
pca_hrf <- multivarious::pca(as.matrix(allhrfs), ncomp=5, preproc=standardize())
pca_hrf <- multivarious::pca(as.matrix(allhrfs), ncomp=5, preproc=pass())
h1=gen_empirical_hrf(seq(0,24,by=1), pca_hrf$u[,1])
h2=gen_empirical_hrf(seq(0,24,by=1), pca_hrf$u[,2])
h3=gen_empirical_hrf(seq(0,24,by=1), pca_hrf$u[,3])
hset <- fmrireg::gen_hrf_set(h1,h2,h3)

betas <- lapply(res, coef)
se <- lapply(res, standard_error)
des <- cells(res[[1]]$model$event_model$terms[[1]])
df2 <- lapply(1:length(sids), function(i) {
  print(i)
  sid = sids[i]
  cf1=reshape_coef(betas[[i]], des) %>% rename(beta=value)
  se1=reshape_coef(se[[i]], des) %>% rename(se=value)
  cf1 <- cf1 %>% mutate(se=se1$se, subject=sid)
  cf1$attention <- as.character(cf1$attention)
  cf1$basis <- as.character(cf1$basis)
  cf1
}) %>% bind_rows()

medse <- df2 %>% group_by(subject, row_id) %>% summarize(medse=median(se))
medse <- medse %>% mutate(flag=ifelse(medse < .001, "bad", "ok"))
df2 <- inner_join(df2, medse, by=c("subject", "row_id"))

df2=inner_join(df2, elab_freq, by=c("subject"="subject", "felab"="felab")) 
df2 <- df2 %>% rename(elabfreq.Var1="elabfreq")
df2 <- df2 %>% rename(elabfreq="elabfreq.Freq")
df2 <- df2 %>% mutate(row_order = row_number())

df2_wide <- df2 %>%
  tidyr::pivot_wider(names_from = row_id, values_from = beta, 
              id_cols = c(felab, attention, basis, elabfreq, subject)) 

write.table(df2_wide, "~/Dropbox/analysis/SOD/fmri/glm_SOD_betas_wide.txt", row.names=FALSE)



options(mc.cores = 1)
library(brms)

df2_roi1 <- dplyr::filter(df2, row_id==310 & flag == "ok")

new_priors <- c(
  prior(normal(0, 12), class = "b"),
  prior(cauchy(0, 4), class = "sd")
)

bres2 <- brm(beta | se(se) ~ felab*attention*basis + (attention*basis|subject),
  data = df2_roi1, 
  family = gaussian(),
  backend = "cmdstanr", 
  #algorithm = "meanfield",
  normalize=FALSE,
  iter=500,
  chains=4,
  #decomp = "QR",
  prior = new_priors
)

library(metafor)
library(lme4)


main_basis <- lapply(1:404, function(i) {
  print(i)
  df2_roi<- dplyr::filter(df2, row_id==i & flag == "ok")
  lme1 <- lmer(beta ~ basis -1 + (1 | subject), data=df2_roi)
  tvals <- coef(summary(lme1))[,4]
  pvals <- coef(summary(lme1))[,5]
  
  data.frame(roi=i, t_basis1=tvals[1], t_basis2=tvals[2], p_basis1=pvals[1], p_basis2=pvals[2])

}) %>% bind_rows()

write.table(main_basis, "~/Dropbox/analysis/SOD/fmri/glm_SOD_main_basis.txt", row.names=FALSE)

attn_by_basis <- lapply(1:404, function(i) {
  print(i)
  df2_roi<- dplyr::filter(df2, row_id==i & flag == "ok")
  es <- rma.mv(yi = beta, 
               V = se^2, 
               mods = ~ attention * basis-1,
               random = list(~ attention | subject, ~ basis | subject), 
               data = df2_roi,
               control=list(iter.max=1000,rel.tol=1e-7))
  
  inter <- anova(es, btt=9:14)
  print(paste("ROI:", i, " pval = ", inter$QMp))
  data.frame(roi=i, logp=-log(inter$QMp), QM=inter$QM, mean_beta=mean(es$beta[10:13,1]))
}) %>% bind_rows()


elab_by_basis <- lapply(1:404, function(i) {
  print(i)
  df2_roi<- dplyr::filter(df2, row_id==i & flag == "ok")
  
  lme2 <- lmer(beta ~ felab*basis-1 + (1 | subject), data=df2_roi)
  #lme1 <- lmer(beta ~ felab+basis + (1 | subject), data=df2_roi)
  lme0 <- lmer(beta ~ basis + (1 | subject), data=df2_roi)
  av <- anova(lme2)
  av2 <- anova(lme2,lme0)
  pval2 <- av2$`Pr(>Chisq)`[2]
  data.frame(roi=i, log_elab=-log(av$`Pr(>F)`[1]), 
             log_basis=-log(av$`Pr(>F)`[2]), 
             log_basis_elab=-log(av$`Pr(>F)`[3]),
             log_either=-log(pval2))
}) %>% bind_rows()

library(lme4)
library(lmerTest)
elab_attention_by_basis <- lapply(1:404, function(i) {
  print(i)
  df2_roi<- dplyr::filter(df2, row_id==i & flag == "ok")
  df2_roi <- df2_roi %>% mutate(elab=scale(as.numeric(as.character(felab))))
  
  
  lme3 <- lmer(beta ~ attention*elab*basis  + (1 | subject), data=df2_roi)
  
  em2 <- emmeans(lme3, "attention", by="basis")
  pem2 <- test(pairs(em2))
  
  em4 <- emtrends(lme3, pairwise ~ attention, var = "elab", by="basis")
  
  elab_basis <- test(emtrends(lme3, ~ 1, var="elab", by="basis"))
  elab_attn_basis <- test(emtrends(lme3, var="elab"))
  elab_attn_conn_basis <- test(em4)
  

  av <- anova(lme3)
  
  data.frame(roi=i, 
             t_elab_b1=elab_basis$t.ratio[1],
             t_elab_b2=elab_basis$t.ratio[2],
             t_elab_b3=elab_basis$t.ratio[3],
             t_elab_da_b1=elab_attn_basis$t.ratio[1],
             t_elab_fa_b1=elab_attn_basis$t.ratio[2],
             t_elab_da_b2=elab_attn_basis$t.ratio[3],
             t_elab_fa_b2=elab_attn_basis$t.ratio[4],
             t_elab_da_b3=elab_attn_basis$t.ratio[5],
             t_elab_fa_b3=elab_attn_basis$t.ratio[6],
             t_da_fa_b1=pem2$t.ratio[1],
             t_da_fa_b2=pem2$t.ratio[2],
             t_da_fa_b3=pem2$t.ratio[3],
             t_elab_da_fa_b1=elab_attn_conn_basis$contrasts$t.ratio[1],
             t_elab_da_fa_b2=elab_attn_conn_basis$contrasts$t.ratio[2],
             t_elab_da_fa_b3=elab_attn_conn_basis$contrasts$t.ratio[3],
             
             log_elab_attention=-log(av$`Pr(>F)`[4]),
             log_attention_basis=-log(av$`Pr(>F)`[6]),
             log_elab_attention_basis=-log(av$`Pr(>F)`[7]))
             
}) %>% bind_rows()
  
write.table(elab_attention_by_basis, "~/Dropbox/analysis/SOD/fmri/lmer_elab_attention_by_basis.txt", row.names=FALSE)

  
  # mm <- mixmeta(beta ~ felab*basis,
  #               S=se^2,
  #               random = ~ felab | subject,
  #               data=df2_roi,
  #               method="ml")
  # es <- rma.mv(yi = beta, 
  #              V = se^2, 
  #              mods = ~ felab * basis-1,
  #              random = list(~ felab | subject, ~ basis | subject), 
  #              data = df2_roi,
  #              control=list(iter.max=200,rel.tol=1e-7))
  # 
  # inter <- anova(es, btt=12:32)
  # print(paste("ROI:", i, " pval = ", inter$QMp))
  # data.frame(roi=i, logp=-log(inter$QMp), QM=inter$QM)
})

library(metafor)
#res2 <- rma(yi = beta, sei = se, data = df2_roi1, mods = ~ attention + basis)
res <- rma.mv(yi = beta, 
              V = se^2, 
              mods = ~ attention * basis, 
              random = list(~ attention | subject, ~ basis | subject), 
              data = df2_roi1, struc="CS")
res0 <- rma.mv(yi = beta, 
              V = se^2, 
              mods = ~ attention + basis, 
              random = list(~ attention | subject, ~ basis | subject), 
              data = df2_roi1)

res01 <- rma.mv(yi = beta, 
               V = se^2, 
               mods = ~ attention*basis, 
               random = list(~ 1 | subject), 
               data = df2_roi1)


# Standard errors of the predictions
SE <- sqrt(diag(newX %*% vcov(model) %*% t(newX)))

# Confidence intervals
CI_lower <- Yhat - 1.96 * SE
CI_upper <- Yhat + 1.96 * SE


get_predict.rma <- function(model, newdata, ...) {
  #browser()
  newX <- model.matrix(model$formula.mods, data = newdata)
  
  Yhat <- newX %*% model$b
  out <- data.frame(
    rowid = seq_len(nrow(Yhat)),
    estimate = as.vector(Yhat))
  return(out)
}

get_coef.rma <- function(model, ...) {
  b <- coef(model)
  b <- setNames(as.vector(b), names(b))
  return(b)
}

set_coef.rma <- function(model, coefs, ...) {
  out <- model
  out$b = as.matrix(coefs)
  out$beta = as.matrix(coefs)
  return(out)
}

get_vcov.rma <- function(model, ...) {
  vcov(model)
}

