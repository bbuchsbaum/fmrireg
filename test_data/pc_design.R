library(tibble)

des <- as_tibble(read.table("test_data/1007_test_design_file.txt", header=TRUE))
des$Input <- factor(des$Saliency)


sframe <- sampling_frame(rep(169,3), TR=1.77)


## Saliency as categorical factor

sal_lin_main <- poly_contrast(~ Saliency, "lin_saliency")
sal_lin_lure <- poly_contrast(~ Saliency, "lin_lure_saliency", where = Repetition == "lure")
sal_lin_new <- poly_contrast(~ Saliency, "lin_new_saliency", where = Repetition == "newtest")
sal_lin_old <- poly_contrast(~ Saliency, "lin_old_saliency", where = Repetition == "old")

sal_lin_old_min_new <- sal_lin_old - sal_lin_new
conlist <- list(
  sal_lin_main,
  sal_lin_lure,
  sal_lin_new,
  sal_lin_old)
  
emodel <- event_model(Onset ~ hrf(Repetition, Input), block = ~ Run, sampling_frame=sframe, data=des, contrasts=conlist)