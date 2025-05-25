## Code to prepare benchmark datasets for the fmrireg package

# Source the benchmark dataset generation script
source("data-raw/generate_benchmark_datasets.R")

# The script will create and save fmri_benchmark_datasets
# No additional usethis::use_data() call needed as it's handled in the generation script
