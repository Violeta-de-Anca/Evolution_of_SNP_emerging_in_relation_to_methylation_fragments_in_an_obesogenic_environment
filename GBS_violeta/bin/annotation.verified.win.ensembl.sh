#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 1
##SBATCH --mem=80gb
#SBATCH -t 4-100:00:00
#SBATCH -J intersect
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/win.intersect.MACS3.per.ind.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/win.intersect.MACS3.per.ind.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load plink
module load BEDTools

working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta
veri_input=$working_dir/verification.peaks.MACS3
gtf_input=$working_dir/reference_genomes/mus_musculus/gtf_files/ensembl/Mus_musculus.GRCm39.111.gtf
inter_out=$veri_input/intersect_win_per_indv

# annotate windows with ensembl GFT file
bedtools intersect -a $file -b $gtf_input > annotated.verified.win.MACS3.bed








