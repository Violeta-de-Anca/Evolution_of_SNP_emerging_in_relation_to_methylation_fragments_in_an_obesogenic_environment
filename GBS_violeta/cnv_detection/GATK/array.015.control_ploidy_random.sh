#!/bin/bash
#SBATCH -A uppmax2025-2-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00
#SBATCH -J jobarray
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/control_random_case_ploidy_array.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/control_random_case_ploidy_array.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

input_path=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK

for F in $(cat $input_path/control_list.txt); do
        echo $F
        sbatch --export=ALL,a=$F 015.random_control_case_mode_ploidy.sh $F
done

