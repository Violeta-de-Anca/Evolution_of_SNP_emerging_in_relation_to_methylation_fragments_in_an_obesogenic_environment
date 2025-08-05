#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core -n 1
#SBATCH -t 10-00:00:00
##SBATCH -t 1:00:00
#SBATCH -J ploidygatk
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/ploidy_control_cohort_mode.gatk.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/ploidy_control_cohort_mode.gatk.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK/4.3.0.0
module load conda/latest
source activate gatk

output_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK
bins_1000=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/bin_1000_GAKT.interval_list

#make the contig prior table
cut -f 1 /proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/mm39.chromAlias.txt > $output_folder/contig_names
tail -n +2 $output_folder/contig_names > $output_folder/contig_names.1
mv $output_folder/contig_names.1 $output_folder/contig_names_no_header
rm $output_folder/contig_names.1
rm $output_folder/contig_names

while read -r line; do
	if [[ "$line" == *chrX* ]]; then
		echo -e "0.01\t0.49\t0.49\t0.01"
	elif [[ "$line" == *chrY* ]]; then
		echo -e "0.5\t0.5\t0.00\t0.00"
	else
		echo -e "0.01\t0.01\t0.97\t0.01"
	fi
done < $output_folder/contig_names_no_header > $output_folder/prob_lines.txt

paste $output_folder/contig_names_no_header $output_folder/prob_lines.txt > $output_folder/contig_prior_table_mm39.tsv
rm $output_folder/prob_lines.txt
rm $output_folder/contig_names_no_header

(echo -e "CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\tPLOIDY_PRIOR_2\tPLOIDY_PRIOR_3"; cat $output_folder/contig_prior_table_mm39.tsv) > $output_folder/contig_prior_table_mm39_with_header.tsv
mv $output_folder/contig_prior_table_mm39_with_header.tsv $output_folder/contig_prior_table_mm39.tsv

gatk DetermineGermlineContigPloidy --output $output_folder --output-prefix control_cohort --contig-ploidy-priors $output_folder/contig_prior_table_mm39.tsv --interval-merging-rule OVERLAPPING_ONLY \
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
