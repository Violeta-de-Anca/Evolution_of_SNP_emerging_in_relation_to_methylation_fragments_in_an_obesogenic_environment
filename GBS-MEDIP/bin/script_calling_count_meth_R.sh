#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 4
#SBATCH --mem=80gb
#SBATCH -t 10:00:00
#SBATCH -J count_matrix
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/meth_matrix_new.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/meth_matrix_new.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load R_packages/4.2.1

working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP

R --no-save --quiet <  $working_dir/bin/R_script_count_matrix.R
