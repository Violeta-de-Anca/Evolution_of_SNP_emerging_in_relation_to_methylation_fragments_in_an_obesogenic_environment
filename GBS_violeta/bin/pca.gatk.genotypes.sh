#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 1
##SBATCH --mem=80gb
#SBATCH -t 4-100:00:00
#SBATCH -J PCA
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/PCA_GBS_ICR.gatk.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/PCA_GBS_ICR.gatk.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load plink

working_dir=/proj/naiss2023-23-55/GBS_violeta
variant_output=$working_dir/variant_output

cd $variant_output
# First we need to put from vcf to bed file the genotypes #
#plink --vcf ICR_genotypes_filtered_gatk.vcf --make-bed --out ICR_genotypes_filtered_gatk --allow-extra-chr

# Then do the PCA #
plink --allow-no-sex --bfile ICR_genotypes_filtered_gatk --pca --out PCA_ICR_GBS_gatk --allow-extra-chr

