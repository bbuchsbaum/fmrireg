library(tibble)

### read in design file
des <- as_tibble(read.table("test_data/1007_test_design_file.txt", header=TRUE))

## convert numeric Saliency to a factor
des$Input <- factor(des$Saliency)


## There are three 'runs' in the recognition phase. Each has 169 scans and a TR = 1.77.
## We construct a 'sampling_frame' which represents the way images are samplid in time and over runs.
sframe <- sampling_frame(rep(169,3), TR=1.77)


## set up linear contrasts
input_lin_main <- poly_contrast(~ Input, degree=2, "lin_input")
input_lin_lure <- poly_contrast(~ Input, "lin_lure_input", where = Repetition == "lure")
input_lin_new <-  poly_contrast(~ Input, "lin_new_input", where = Repetition == "newtest")
input_lin_old <-  poly_contrast(~ Input, "lin_old_input", where = Repetition == "old")

## differences if linear contrasts
input_lin_old_min_new <- input_lin_old - input_lin_new
input_lin_old_min_lure <- input_lin_old - input_lin_lure
input_lin_lure_min_new <- input_lin_lure - input_lin_new

## pairwise contrasts of probe type, averaging over 'input'
old_new  <-  pair_contrast(~ Repetition == "old", ~ Repetition == "newtest", name="old_min_new")
new_lure <-  pair_contrast(~ Repetition == "newtest", ~ Repetition=="lure", name="new_min_lure")
old_lure <-  pair_contrast(~ Repetition == "old",  ~ Repetition == "lure", name="old_min_lure")

conlist <- contrast_set(
  input_lin_main,
  input_lin_lure,
  input_lin_new,
  input_lin_old,
  input_lin_old_min_new,
  input_lin_old_min_lure,
  input_lin_lure_min_new,
  old_new,
  new_lure,
  old_lure)

## construct event model
emodel <- event_model(Onset ~ hrf(Repetition, Input, contrasts=conlist), block = ~ Run, sampling_frame=sframe, data=des)

## construct baseline model
bmodel <- baseline_model(basis="bs", degree=5, sframe=sframe)

## put them together
fmodel <- fmri_model(emodel, bmodel)

## contruct fmri_dataset
dset <- fmri_dataset(
  scans=c("../epi/rscan02.nii", "../epi/rscan04.nii","../epi/rscan06.nii"),
  mask="../epi/global_mask.nii",
  TR=1.77,
  run_length=169,
  event_table=des,
  base_path="."
)

## convert to an AFNI model
alm <- afni_lm(fmodel, dset)

## execute the AFNI command
run(alm, outdir="glm_input_reptype", execute=FALSE)
          