#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J filtergatk
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/filter_random_counts.gatk.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/filter_random_counts.gatk.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK/4.3.0.0
module load conda/latest
source activate gatk

# Check if the GATK conda environment is loaded
[[ -z "$CONDA_DEFAULT_ENV" ]] && echo "No Conda environment is currently active." || echo "The currently active Conda environment is: $CONDA_DEFAULT_ENV"

# Check if gcnvkernel is available
python -c "import gcnvkernel" &> /dev/null && echo "Module gcnvkernel exists." || echo "Module gcnvkernel does not exist."

output_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK
bins_1000=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/bin_1000_GAKT.interval_list

gatk FilterIntervals --intervals $bins_1000 --output $output_folder/filtered_random_control.interval_list \
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
        --input $output_folder/C13F3_9_Mouse.unique.sorted.bam.HDF5 \
	--input $output_folder/random_control_1.HDF5 \
	--input $output_folder/random_control_2.HDF5 \
	--input $output_folder/random_control_3.HDF5 \
	--input $output_folder/random_control_4.HDF5 \
	--input $output_folder/random_control_5.HDF5 \
	--input $output_folder/random_control_6.HDF5 \
	--input $output_folder/random_control_7.HDF5 \
	--input $output_folder/random_control_8.HDF5 \
	--input $output_folder/random_control_9.HDF5 \
	--input $output_folder/random_control_10.HDF5 \
	--input $output_folder/random_control_11.HDF5 \
	--input $output_folder/random_control_12.HDF5 \
	--input $output_folder/random_control_13.HDF5 \
	--input $output_folder/random_control_14.HDF5 \
	--input $output_folder/random_control_15.HDF5 \
	--input $output_folder/random_control_16.HDF5 \
	--input $output_folder/random_control_17.HDF5 \
	--input $output_folder/random_control_18.HDF5 \
	--input $output_folder/random_control_19.HDF5 \
	--input $output_folder/random_control_20.HDF5 \
	--input $output_folder/random_control_21.HDF5 \
	--input $output_folder/random_control_22.HDF5 \
	--input $output_folder/random_control_23.HDF5 \
	--input $output_folder/random_control_24.HDF5 \
	--input $output_folder/random_control_25.HDF5 \
	--input $output_folder/random_control_26.HDF5 \
	--input $output_folder/random_control_27.HDF5
