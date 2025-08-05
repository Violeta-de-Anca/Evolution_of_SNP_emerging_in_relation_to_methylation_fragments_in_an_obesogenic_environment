#!/bin/bash -l
#SBATCH -A uppmax2025-2-151
#SBATCH -p core -n 1
#SBATCH -t 10-000:00:00
#SBATCH -J flagsta
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/control_stats_cnv.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/control_stats_cnv.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14

output_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/aligned
folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK

while read -r line; do
	a=${line##*/}
	samtools flagstat $line > $folder/$a.flagstats
done < $output_folder/control.bam.txt
