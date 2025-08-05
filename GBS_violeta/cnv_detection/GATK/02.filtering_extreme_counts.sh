#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J filtergatk
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/filter_counts.gatk.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/filter_counts.gatk.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK/4.3.0.0
module load conda/latest
source activate gatk

output_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK
bins_1000=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/bin_1000_GAKT.interval_list

gatk FilterIntervals --intervals $bins_1000 --output $output_folder/filtered_gatk_mm39_bins.interval_list --input C13F1_10_Mouse.unique.sorted.bam.HDF5 \
	--interval-merging-rule OVERLAPPING_ONLY \
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
	--input $output_folder/C13F2_12_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F2_13_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F2_14_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F2_15_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F2_16_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F2_17_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F2_18_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F2_19_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F2_1_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F2_20_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F2_21_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F2_23_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F2_24_Mouse.unique.sorted.bam.HDF5 \
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
	--input $output_folder/C13F3_14_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_15_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_16_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_17_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_18_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_19_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_1_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_20_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_21_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_22_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_24_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_25_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_27_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_2_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_3_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_4_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_5_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_7_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_8_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/C13F3_9_Mouse.unique.sorted.bam.HDF5
