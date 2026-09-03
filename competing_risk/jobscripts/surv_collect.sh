#!/bin/bash
# mem=32gb, well above continuous's 10gb: surv_analysis.R attaches the full
# simulated dataset and truth to every per-run result (results$data,
# results$truth), not just model outputs, and this study has 7,000 runs (7
# scenarios x 2 censoring x 500) vs. continuous's 4,000. Locally, the
# collected surv_all.RDS is already ~729MB serialized, which get_results()
# has to hold unserialized (larger) in memory alongside every per-run file
# it reads - this figure is a reasoned estimate, not a profiled one (this
# study has no *_profile.R step like continuous/crossfitting do). Check
# actual job stats after the first real submission and tighten it.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=32gb
#PBS -N surv_collect


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript surv_collect.R
