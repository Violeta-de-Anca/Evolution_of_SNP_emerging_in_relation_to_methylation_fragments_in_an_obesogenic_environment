#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 4
#SBATCH -t 100:00:00
#SBATCH -J pca
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/angsd.pca.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/angsd.pca.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools samtools PCAngsd ANGSD bcftools plink htslib

working_dir=/proj/naiss2023-23-55/GBS_violeta
input_prefix=$working_dir/genotype_likelihood
db_knownsites=$working_dir/snp_indel_databases
pca=$working_dir/pca
reference_genome=/proj/naiss2023-23-55/GBS_violeta/reference_genomes/mus_musculus/uscs_ref/mm39.fa
mkdir -p $pca

cd $pca

pcangsd -b $input_prefix/angsd.GBS.ICR.sperm.beagle.gz -o PCA.angsd.GBS.ICR.sperm
