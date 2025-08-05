#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core -n 4
##SBATCH --mem=80gb
#SBATCH -t 10-00:00:00
#SBATCH -J bingatk
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/bin_1000bp.gatk.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/bin_1000bp.gatk.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK/4.3.0.0
module load conda/latest
source activate gatk

folder=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref

#IF YOU DON'T PUT THE EXTENSION .interval_list YOU WILL GET AN ERROR IN THE NEXT STEP!!!!!!!!

gatk PreprocessIntervals --bin-length 1000 --output $folder/bin_1000_GAKT.interval_list --padding  0 --interval-merging-rule OVERLAPPING_ONLY --reference $folder/mm39.fa
