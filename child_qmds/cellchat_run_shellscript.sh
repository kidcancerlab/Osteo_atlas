#!/bin/bash
#SBATCH --job-name=cellchat_sharedgenes
#SBATCH --output=slurmout/cellchat_sharedgenes_%j.out
#SBATCH --error=slurmout/cellchat_sharedgenes_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=10
#SBATCH --partition=himem,general
#SBATCH --wait

set -e
#read args
ARGS=("$@")

module load R/4.3.0

Rscript /home/gdrobertslab/lab/Analysis/Yogesh/CellTypeAnnRefs/child_qmds/cellchat_slurm_job.R \
    $SLURM_ARRAY_TASK_ID \
    ${ARGS[1]}

