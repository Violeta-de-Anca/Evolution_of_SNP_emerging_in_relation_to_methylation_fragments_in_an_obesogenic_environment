#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core -n 4
#SBATCH -t 10-00:00:00
##SBATCH -t 1:00:00
#SBATCH -J cnvgatk
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/cnv_control_cohort_mode.gatk.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/cnv_control_cohort_mode.gatk.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK/4.3.0.0
module load conda/latest
source activate gatk

output_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK
bins_1000=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/bin_1000_GAKT.interval_list

#First we need to create the model for the calling of the CNVs, so only inputting the control individuals
gatk GermlineCNVCaller --contig-ploidy-calls $output_folder/control_cohort-calls \
	--run-mode COHORT \
	--intervals filtered_gatk_mm39_bins.interval_list \
	--interval-merging-rule OVERLAPPING_ONLY \
	--output $output_folder \
	--output-prefix control_model_cnv_caller_GATK \
        --input $output_folder/C13F1_11_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F1_1_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F1_2_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F1_3_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F1_4_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F1_5_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F1_6_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F1_7_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F1_8_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F1_9_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F2_2_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F2_3_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F2_4_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F2_5_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F2_6_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F2_7_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F2_8_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F2_9_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F3_10_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F3_11_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F3_12_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F3_13_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F3_1_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F3_2_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F3_3_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F3_4_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F3_5_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F3_7_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F3_8_Mouse.unique.sorted.bam.HDF5 \
        --input $output_folder/C13F3_9_Mouse.unique.sorted.bam.HDF5

