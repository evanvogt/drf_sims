#!/bin/bash
# Same sizing rationale as me_collect.sh - see README.md.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=20gb
#PBS -N me_metrics


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

# Which result tree to process - pass with qsub, e.g.
#   qsub -v TREE=strategies me_metrics.sh
# Unset means the main study, which is what the bare Rscript call did
# before the derived trees existed.
Rscript me_metrics.R ${TREE}
