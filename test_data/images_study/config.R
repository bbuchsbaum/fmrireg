
base_path="/Users/brad/code/fmrireg/test_data/images_study"

scans = paste0("epi/",
               c("rscan01.nii",
               "rscan02.nii",
               "rscan03.nii",
               "rscan04.nii",
               "rscan05.nii",
               "rscan06.nii"))

nuisance_reg = paste0("epi/",
                      c("nuisance_rscan01.txt",
                        "nuisance_rscan02.txt",
                        "nuisance_rscan03.txt",
                        "nuisance_rscan04.txt",
                        "nuisance_rscan05.txti",
                        "nuisance_rscan06.txt"))

design = "behavior/design.txt"

model = ~ onsetTime ~ hrf(imageName) | block
                      