#!/bin/bash
#SBATCH --account=gdrobertslab
#SBATCH --output=output/spacexr/%j.txt
#SBATCH --error=output/spacexr/%j.txt
#SBATCH --job-name rctd
#SBATCH --array=0-40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --partition=himem,general
#SBATCH --time=24:00:00
#SBATCH --wait

# Get arrays holding object name and reference name
while read object_name ref_name
do
    obj_array+=("$object_name")
    ref_array+=("$ref_name")
done < misc/rctd_input.txt

ob_id=${obj_array[$SLURM_ARRAY_TASK_ID]}
ref=${ref_array[$SLURM_ARRAY_TASK_ID]}

echo Annotating $ob_id with $ref

# move and rename output to somewhere findable
mv output/spacexr/${SLURM_JOBID}.txt output/spacexr/granular_references/$ref/${ob_id}_output_multi.txt

Rscript rctd_scripts/rctd_multi.r $ob_id $ref