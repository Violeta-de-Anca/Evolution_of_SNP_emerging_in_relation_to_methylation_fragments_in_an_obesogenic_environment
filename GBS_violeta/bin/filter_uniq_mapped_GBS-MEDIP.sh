#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 4
#SBATCH --mem=80gb
#SBATCH -t 100:00:00
#SBATCH -J filterMQ20
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/filter_MQ20.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/filter_MQ20.out

module load bioinfo-tools
module load samtools/1.14
module load FastQC/0.11.9
module load MultiQC

working_dir=/proj/ancestry_medi_indiv/aligned
output_dir=/proj/naiss2023-23-55/GBS_violeta/gbs-medip

cd $output_dir
for i in *Mouse.unique.bam; do
  sample=$(basename $i .unique.bam)
#  samtools view -h $i | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > $output_dir/$sample.unique.bam
  samtools view -f 2 -F 524 -q 10 -b -o $output_dir/$sample.unique.MQ20.bam $i
done


