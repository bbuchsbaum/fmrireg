library(tibble)

des <- as_tibble(read.table("test_data/1007_test_design_file.txt", header=TRUE))
des$Input <- factor(des$Saliency)


sframe <- sampling_frame(rep(169,3), TR=1.77)


## Saliency as categorical factor

input_lin_main <- poly_contrast(~ Input, "lin_input")
input_lin_lure <- poly_contrast(~ Input, "lin_lure_input", where = Repetition == "lure")
input_lin_new <- poly_contrast(~ Input, "lin_new_input", where = Repetition == "newtest")
input_lin_old <- poly_contrast(~ Input, "lin_old_input", where = Repetition == "old")

input_lin_old_min_new <- input_lin_old - input_lin_new
conlist <- contrast_set(
  input_lin_main,
  input_lin_lure,
  input_lin_new,
  input_lin_old)
  
emodel <- event_model(Onset ~ hrf(Repetition, Input, contrasts=conlist), block = ~ Run, sampling_frame=sframe, data=des)
bmodel <- baseline_model(basis="bs", degree=5, sframe=sframe)
fmodel <- fmri_model(emodel, bmodel)
