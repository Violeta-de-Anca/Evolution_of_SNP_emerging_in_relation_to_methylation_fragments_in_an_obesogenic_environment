#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J subgatk
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/clean_up_all_for_plotting.gatk.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/clean_up_all_for_plotting.gatk.out
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
#for i in segment_F*_cnv.vcf; do
#	a=${i%_cnv.vcf}
#	b=${a#segment_}
#	grep -E '<DUP>|<DEL>' $i > $output_folder/${b}_cnv_all_sizes.vcf
#
#done

#then we need to clean the end of the cnv
#for i in segment_F*_cnv.vcf; do
#	a=${i%_cnv.vcf}
#	b=${a#segment_}
#	awk 'BEGIN {OFS="\t"} {gsub("END=", "", $8); print}' ${b}_cnv_all_sizes.vcf > $output_folder/${b}_cnv_intermediate.vcf
#done

#then we need to see how long is the CNV
#for i in segment_F*_cnv.vcf; do
#	a=${i%_cnv.vcf}
#	b=${a#segment_}
#	awk 'BEGIN {OFS="\t"} {diff = $8 - $2; print $0, diff}' ${b}_cnv_intermediate.vcf > $output_folder/${b}_cnv_intermediate_lenght.vcf
#done

#subset the CNV that are at least 50Kbps
#for i in segment_F*_cnv.vcf; do
#	a=${i%_cnv.vcf}
#	b=${a#segment_}
#	awk '$11 >= 50000' ${b}_cnv_intermediate_lenght.vcf.vcf > ${b}_filtered_cnv.vcf
#done

#remove intermediate files
#for i in segment_F*_cnv.vcf; do
#	a=${i%_cnv.vcf}
#	b=${a#segment_}
#	rm ${b}_cnv_intermediate_lenght.vcf.vcf
#	rm ${b}_cnv_intermediate.vcf
#	rm ${b}_cnv_intermediate.vcf
#	rm ${b}_cnv_all_sizes.vcf
#done

#do a simpler bed file for plotting (F2_17_filtered_cnv.vcf)
for i in segment_F*_cnv.vcf; do
	a=${i%_cnv.vcf}
	b=${a#segment_}
	awk '{print $1, $2, $5, $8}' OFS="\t" ${b}_filtered_cnv.vcf > $output_folder/${b}_filtered_cnv.bed
done
