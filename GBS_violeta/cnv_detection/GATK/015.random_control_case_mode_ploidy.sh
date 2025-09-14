#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core -n 1
#SBATCH -t 10-00:00:00
##SBATCH -t 1:00:00
#SBATCH -J casegatk
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/ploidy_random_case_mode.gatk.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/ploidy_random_case_mode.gatk.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK/4.3.0.0
module load conda/latest
source activate gatk

output_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK

#do it one by one in an array of samples
sample=$1
input_folder=${sample%/*}
echo $input_folder
name_sample=${sample##*/}
echo $name_sample
cleaned=${name_sample#C13}
cleaned_1=${cleaned%_Mouse.unique.sorted.bam.HDF5}

gatk DetermineGermlineContigPloidy --output $output_folder --output-prefix control_random_ploidy_${cleaned_1} --model $output_folder/random_cohort-model \
	--input $output_folder/$sample
