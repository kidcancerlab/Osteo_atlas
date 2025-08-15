#!/bin/bash
#SBATCH --account=gdrobertslab
#SBATCH --output=output/spacexr/granular_references/rctd_%j.txt
#SBATCH --error=output/spacexr/granular_references/rctd_%j.txt
#SBATCH --job-name rctd
#SBATCH --array=0-39
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --partition=himem,general
#SBATCH --time=24:00:00
#SBATCH --wait

set -e

# Get arrays holding object name and reference name
while read object_name ref_name
do
    obj_array+=("$object_name")
    ref_array+=("$ref_name")
done < misc/rctd_input.txt

ob_id=${obj_array[$SLURM_ARRAY_TASK_ID]}
ref=${ref_array[$SLURM_ARRAY_TASK_ID]}

echo Annotating $ob_id with $ref

Rscript rctd_scripts/rctd.r $ob_id $ref
