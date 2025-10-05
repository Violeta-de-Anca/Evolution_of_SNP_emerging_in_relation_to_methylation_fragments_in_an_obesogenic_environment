#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J subgatk
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/clean_up_control_random_for_plotting.gatk.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/clean_up_control_random_for_plotting.gatk.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK/4.3.0.0
module load conda/latest
source activate gatk
source $GATK_HOME/gatk-completion.sh

output_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK
bins_1000=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/bin_1000_GAKT.interval_list
dictionary=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/mm39.dict
ref=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/mm39.fa

#first we need to subset the segments with deletions or duplications
#for i in $output_folder/segment_control_random_F*vcf; do
#	g=${i##*/}
#	a=${g%_cnv.vcf}
#	b=${a#segment_control_random_}
#	grep -E '<DUP>|<DEL>' $i > $output_folder/${b}_control_random_cnv_all_sizes.vcf
#done

#second we need to clean the columns
#for i in $output_folder/segment_control_random_F*vcf; do
#	g=${i##*/}
#	a=${g%_cnv.vcf}
#	b=${a#segment_control_random_}
#	awk 'BEGIN {OFS="\t"} {gsub("END=", "", $8); print}' ${b}_control_random_cnv_all_sizes.vcf > $output_folder/${b}_cnv_control_random_intermediate.vcf
#done

#third is calculate the length of the CNV so we can subset
#for i in $output_folder/segment_control_random_F*vcf; do
#	g=${i##*/}
#	a=${g%_cnv.vcf}
#	b=${a#segment_control_random_}
#	awk 'BEGIN {OFS="\t"} {diff = $8 - $2; print $0, diff}' ${b}_cnv_control_random_intermediate.vcf > $output_folder/${b}_cnv_control_random_lenght.vcf
#done

#subset the CNVs above 50Kbps
for i in $output_folder/segment_control_random_F*vcf; do
	g=${i##*/}
	a=${g%_cnv.vcf}
	b=${a#segment_control_random_}
	awk '$11 >= 50000' ${b}_cnv_control_random_lenght.vcf > $output_folder/${b}_control_filtered_cnvs.vcf
done

#remove intermediate files
for i in $output_folder/segment_control_random_F*vcf; do
	g=${i##*/}
	a=${g%_cnv.vcf}
	b=${a#segment_control_random_}
	rm $output_folder/${b}_control_random_cnv_all_sizes.vcf
	rm $output_folder/${b}_cnv_control_random_intermediate.vcf
	rm $output_folder/${b}_cnv_control_random_lenght.vcf
done
