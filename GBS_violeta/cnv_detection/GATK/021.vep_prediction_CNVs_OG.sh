#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J vep
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/vep_prediction_CNV_OG.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/vep_prediction_CNV_OG.gatk.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools vep/113.0

#variant effect predictor, for the final putative CNVs!

#vep -v -i /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK/F2_cnvs_final.txt -o F2_predictions_VEP.txt --cache --dir $VEP_CACHE --force_overwrite  --assembly GRCm39 --offline --species mus_musculus --format ensembl --stats_text --stats_file F2_stats_file.txt --stats_html

vep -v -i /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK/F3_cnvs_final.ensembl -o F3_predictions_VEP.txt --cache --dir $VEP_CACHE --force_overwrite  --assembly GRCm39 --offline --species mus_musculus --format ensembl --stats_text --stats_file F3_stats_file.txt --stats_html
