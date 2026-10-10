#!/bin/bash
# The DGM check (cts_val_dgm_check.R) as one job: its runs are spread over
# NCPUS future workers, one grf thread each. Submit from this folder:
#   qsub cts_val_dgm_check.sh                                  # scenario 2, rho 0, 100 runs
#   qsub -v SCENARIO=3,RHO=0.5 cts_val_dgm_check.sh
#   qsub -v RUNS=50,VIMS=vims -l walltime=12:00:00 cts_val_dgm_check.sh
# Finished runs are cached, so a job that hits its walltime is resumed by
# submitting it again with the same variables. Tables land in the .o file and
# as CSVs under <metrics folder>/investigation/dgm_check/.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=8:ompthreads=8:mem=32gb
#PBS -N cts_val_dgm_check

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script - workers = the cores PBS allocated (detectCores() would see the
# whole node)
Rscript cts_val_dgm_check.R "${RUNS:-100}" "$NCPUS" ${VIMS:-} \
  scenario="${SCENARIO:-2}" rho="${RHO:-0}"
