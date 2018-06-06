
base_path="/Users/bbuchsbaum/code/fmrireg/test_data/images_study"

scans = paste0("epi/",
               c("rscan01.nii",
                 "rscan02.nii",
                 "rscan03.nii",
                 "rscan04.nii",
                 "rscan05.nii",
                 "rscan06.nii"))


event_table = "behavior/design.txt"

aux_table = "epi/aux_data.txt"

block_column = "run"

event_model = onsetTime ~ hrf(imageName, id="iname", subset=!is.na(imageName)) 

baseline_model = ~ splines::bs(.scan_time,4) 

output_dir = "glm_out"

mask = "epi/global_mask.nii"

TR = 1.5

run_length = 348
                      