#!/bin/bash -l
#SBATCH -A naiss2024-22-316
#SBATCH -p core -n 4
#SBATCH --mem=80gb
#SBATCH -t 10:00:00
#SBATCH -J ANUBIX_links
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/pathway_enrichment.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/pathway_enrichment.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load R_packages/4.2.1

working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP

R --no-save --quiet <  $working_dir/bin/pathway_enrichment.R
