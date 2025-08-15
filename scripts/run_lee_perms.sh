#!/bin/sh
#SBATCH --account=gdrobertslab
#SBATCH --output=slurmOut/lee_perms-%j.out
#SBATCH --error=slurmOut/lee_perms-%j.out
#SBATCH --job-name=lee_perms
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=general,himem
#SBATCH --wait

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

# using slurm array task id to match the sample number
Rscript scripts/run_lee_perms.R $SLURM_ARRAY_TASK_ID
