#!/bin/bash -l
#SBATCH -A naiss2023-22-848
#SBATCH -p core -n 1
#SBATCH -t 10:00:00
#SBATCH -J cnv_corr
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/cnv_correlation.on.parental.meth.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files//cnv_correlation.on.parental.meth.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load plink
module load BEDTools
module load bcftools
module load java/OpenJDK_17+35
module load vcftools BEDOPS
module load modkit
module load python3

working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta
input=$working_dir/dynamics_CNVs

for suffix in f0_10 f1_13 f1_14 f1_15 f1_16 f1_1 f1_2; do
	bedtools intersect -a $input/log2_ratios_${suffix}.txt -b $input/${suffix}.normalized.counts.txt -loj > $input/${suffix}_table_cnvs_meth.txt
done

#for suffix in f0_10 f1_13 f1_14 f1_15 f1_16 f1_1 f1_2; do
#	cut -f 1,2,3,4,5,6,11 $input/${suffix}_table_cnvs_meth.txt > $input/${suffix}_final_columns.txt
#	mv $input/${suffix}_final_columns.txt $input/${suffix}_table_cnvs_meth.txt
#done
