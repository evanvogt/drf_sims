#!/bin/bash
# PLACEHOLDER resources below - this study has never completed a run, so
# unlike every other jobscript in this repo these numbers are not measured.
# Run the profiling workflow first, which overwrites the #PBS -l lines and
# the trailing Rscript args together (see me_profile_summary.R):
#
#   Rscript model_evaluation/me_testing.R full        # confirms it runs at all, first
#   Rscript model_evaluation/me_profile.R 1            # smoke-test the profiler locally
#   qsub model_evaluation/jobscripts/me_profile.sh     # 16 profiling jobs
#   Rscript model_evaluation/me_profile_summary.R      # writes measured directives below
#
# The %N array throttle also still needs setting from the HPC queue's real
# memory/fair-share limits before the first real submission - each
# concurrent task starts its own H2O JVM cluster (mem="10G" heap), which
# rules out anything resembling continuous's %190 or crossfitting's %380.
# See README.md.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=10gb
#PBS -J 1-360%4
#PBS -N me_1
#PBS -o logs_1/
#PBS -e logs_1/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters. Trailing args are workers/n_cores -
# placeholders here, overwritten by me_profile_summary.R alongside the
# #PBS -l lines above.
Rscript me_analysis.R "$PBS_ARRAY_INDEX" 2 2
