#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 3
#SBATCH -t 100:00:00
#SBATCH -J checkrefgen
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/checking_refgen.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/checking_refgen.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK
module load picard
reference_genome=/proj/naiss2023-23-55/GBS_violeta/reference_genomes/mus_musculus/uscs_ref/mm39.fa.gz

gatk QCRef
