
base_path="/Users/brad/code/fmrireg/test_data/images_study"

scans = paste0("epi/",
               c("rscan01.nii",
                 "rscan02.nii",
                 "rscan03.nii",
                 "rscan04.nii",
                 "rscan05.nii",
                 "rscan06.nii"))




design = "behavior/design.txt"

aux_data = "epi/aux_data.txt"

#block_column = "run"

event_model = onsetTime ~ hrf(imageName) | run

baseline_model = ~ bs(time,4) + PC1 + PC2 + PC3

output_dir = "glm_out"

mask = "epi/global_mask.nii"

TR = 1.5

run_length = 348
                      