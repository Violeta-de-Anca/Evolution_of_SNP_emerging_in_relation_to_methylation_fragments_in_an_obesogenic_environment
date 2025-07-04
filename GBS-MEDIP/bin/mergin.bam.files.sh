#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 2
##SBATCH --mem=80gb
#SBATCH -t 100:00:00
#SBATCH -J sort&ind
#SBATCH --error /proj/ancestry_medi_indiv/log_files/merginbams.MQ20.err
#SBATCH --output /proj/ancestry_medi_indiv/log_files/merginbams.MQ20.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load bamtools

gbs_work_dir=/proj/naiss2023-23-55/GBS_violeta
medip_work_dir=/proj/ancestry_medi_indiv/aligned

#bamtools merge -list $gbs_work_dir/bamuniquelist.txt -out $gbs_work_dir/merged.GBS.uniq.bam
bamtools merge -list $gbs_work_dir/gbs-medip/bam.MQ20.list.txt -out $gbs_work_dir/gbs-medip/merged.GBS_MeDIP.MQ20.bam

#samtools sort $gbs_work_dir/merged.GBS.uniq.bam -o $gbs_work_dir/merged.GBS.uniq.sorted.bam
samtools sort $gbs_work_dir/gbs-medip/merged.GBS_MeDIP.MQ20.bam -o $gbs_work_dir/gbs-medip/merged.GBS_MeDIP.MQ20.sorted.bam

#samtools index $gbs_work_dir/merged.GBS.uniq.sorted.bam $gbs_work_dir/merged.GBS.uniq.sorted.bam.bai
samtools index $gbs_work_dir/gbs-medip/merged.GBS_MeDIP.MQ20.sorted.bam $gbs_work_dir/gbs-medip/merged.GBS_MeDIP.MQ20.sorted.bam.bai

#mv $gbs_work_dir/merged.GBS_MeDIP.sorted.bam $medip_work_dir/.
