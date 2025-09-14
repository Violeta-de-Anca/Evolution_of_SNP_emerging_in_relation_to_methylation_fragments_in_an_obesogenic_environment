#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core -n 1
#SBATCH -t 10-00:00:00
##SBATCH -t 1:00:00
#SBATCH -J casegatk
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/ploidy_case_mode.gatk.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/ploidy_case_mode.gatk.out
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

gatk DetermineGermlineContigPloidy --output $output_folder --output-prefix overnutrition_${cleaned_1} --model $output_folder/control_cohort-model --input $sample

#############################################################
#THE CODE BELOW DOES NOT WORK!!!!!!!!!!!!!!!
#############################################################
#this does not work with all the individuals at the same time
#gatk DetermineGermlineContigPloidy --output $output_folder --output-prefix overnutrition_ploidy --model $output_folder/control_cohort-model \
#	--input $output_folder/C13F2_13_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F2_14_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F2_15_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F2_16_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F2_17_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F2_18_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F2_19_Mouse.unique.sorted.bam.HDF5 \
#	--input $output_folder/C13F2_20_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F2_21_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F2_23_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F2_24_Mouse.unique.sorted.bam.HDF5 \
#	--input $output_folder/C13F3_14_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F3_15_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F3_16_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F3_17_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F3_18_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F3_19_Mouse.unique.sorted.bam.HDF5 \
#	--input $output_folder/C13F3_20_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F3_21_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F3_22_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F3_24_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F3_25_Mouse.unique.sorted.bam.HDF5 \
#        --input $output_folder/C13F3_27_Mouse.unique.sorted.bam.HDF5

