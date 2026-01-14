require(here)

# rproj root - folder containing scripts
scripts_dir <- here()

# data and results folders are at the same level as the scripts folder
data_dir <- file.path(dirname(scripts_dir), "data")
results_dir <- file.path(dirname(scripts_dir), "results")
