#!/bin/bash
#PBS -l walltime=24:00:00  
#PBS -l select=1:ncpus=10:ompthreads=10:mem=20gb
#PBS -J 1-500%43
#PBS -N cts_miss_ci
#PBS -o logs/
#PBS -e logs/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters. Trailing args are workers[1]/workers[2]/
# grf_threads, overwritten by cts_miss_ci_profile_summary.R alongside the
# #PBS -l lines above once jobscripts/cts_miss_ci_profile.sh's sweep has run -
# do not hand-edit one without the other, run the sweep instead. Below is
# still the pre-profiling placeholder (workers unchanged at 3/3; grf_threads=1
# is a guess - today's true behaviour is unset/grf-default, which isn't
# expressible as a CLI arg).
Rscript cts_miss_ci_analysis.R "$PBS_ARRAY_INDEX" 3 3 1
