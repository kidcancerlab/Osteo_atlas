#!/bin/bash
#SBATCH --account=gdrobertslab
#SBATCH --output=DoubletFinder/slurm/%j.txt
#SBATCH --error=DoubletFinder/slurm/%j.txt
#SBATCH --job-name demuxafy
#SBATCH --array=0-153%50
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --partition=himem,general
#SBATCH --time=24:00:00
#SBATCH --wait


#Get arrays holding sample name, path to counts folder, and path to output folder
while read sample_name count_path
do
    sample_array+=("$sample_name")
    count_array+=("$count_path")
done < misc/doubletfinder_input.txt

count_mat=${count_array[$SLURM_ARRAY_TASK_ID]}
sample_id=${sample_array[$SLURM_ARRAY_TASK_ID]}
out_path=DoubletFinder/output/$sample_id


echo Sample is $sample_id. Moving slurm output to slurm/$sample_id.txt

mv DoubletFinder/slurm/${SLURM_JOBID}.txt DoubletFinder/slurm/$sample_id.txt

mkdir -p $out_path/DoubletDetection $out_path/scDblFinder $out_path/scds $out_path/combined

singularity exec --bind /home/gdrobertslab/lab/Analysis/MattGust/projects/Roberts_Lab/Osteo_atlas,/home/gdrobertslab/lab/Counts_2/ \
    /home/gdrobertslab/lab/Tools/Demuxafy.sif \
    DoubletDetection.py \
    -m $count_mat/ \
    -o $out_path/DoubletDetection 

echo Finished DoubletDetection, running scDblFinder

singularity exec --bind /home/gdrobertslab/lab/Analysis/MattGust/projects/Roberts_Lab/Osteo_atlas,/home/gdrobertslab/lab/Counts_2/ \
    /home/gdrobertslab/lab/Tools/Demuxafy.sif \
    scDblFinder.R \
    -t $count_mat/ \
    -o $out_path/scDblFinder


echo Finished scDblFinder, running scds...

singularity exec --bind /home/gdrobertslab/lab/Analysis/MattGust/projects/Roberts_Lab/Osteo_atlas,/home/gdrobertslab/lab/Counts_2/ \
    /home/gdrobertslab/lab/Tools/Demuxafy.sif \
    scds.R \
    -o $out_path/scds \
    -t $count_mat

echo Finished scds, combining results...

singularity exec --bind /home/gdrobertslab/lab/Analysis/MattGust/projects/Roberts_Lab/Osteo_atlas \
    /home/gdrobertslab/lab/Tools/Demuxafy.sif \
    Combine_Results.R \
    -o $out_path/combined/combined_doublets.tsv \
    -t $out_path/DoubletDetection \
    -n $out_path/scDblFinder \
    -c $out_path/scds \
    -m 'MajoritySinglet'

for sample in ${sample_array[@]}; do
    echo \n $sample \n > ls_files.txt
    ls -lht DoubletFinder/output/$sample/combined >> ls_files.txt
done