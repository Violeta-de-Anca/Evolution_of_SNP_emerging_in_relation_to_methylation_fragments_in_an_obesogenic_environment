#!/bin/bash -l
#SBATCH -A naiss2023-22-848
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J vsearch
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/vsearch_GBS_MEDIP.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/vsearch_GBS_MEDIP.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

#load modules
module load bioinfo-tools


./proj/naiss2024-23-57/vsearch-2.28.1-linux-x86_64/bin/vsearch --derep_fulllength [input.fa] --output [output.fa] -sizeout --uc uc_out
