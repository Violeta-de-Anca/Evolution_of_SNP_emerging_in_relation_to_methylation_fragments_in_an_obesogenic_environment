#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 4
#SBATCH --mem=80gb
#SBATCH -t 100:00:00
#SBATCH -J genotypes_likelihoods
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/genotypes_likelihoods.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/genotypes_likelihoods.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools samtools PCAngsd ANGSD bcftools plink htslib

working_dir=/proj/naiss2023-23-55/GBS_violeta
input_prefix=$working_dir/aligned
genotype_likelihood_output=$working_dir/genotype_likelihood
reference_genome=/proj/naiss2023-23-55/GBS_violeta/reference_genomes/mus_musculus/uscs_ref/mm39.fa
mkdir -p $genotype_likelihood_output

cd $genotype_likelihood_output

angsd -GL 2 -doGlf 2 -b $input_prefix/bamlist.txt -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -nThreads 4
