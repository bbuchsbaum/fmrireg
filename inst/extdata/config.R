
base_path="."

scans = paste0("epi/",
               c("rscan01.nii",
                 "rscan02.nii",
                 "rscan03.nii",
                 "rscan04.nii",
                 "rscan05.nii",
                 "rscan06.nii"))


event_table = "face_design.txt"

block_column = "run"

event_model = onsetTime ~ hrf(imageName, id="iname", subset=!is.na(imageName)) 

baseline_model = list(
  basis="bs",
  degree=5,
  nuisance_list="nuisance.txt"
)
  
output_dir = "glm_out"

mask = "epi/global_mask.nii"

TR = 1.5

run_length = 348
                      