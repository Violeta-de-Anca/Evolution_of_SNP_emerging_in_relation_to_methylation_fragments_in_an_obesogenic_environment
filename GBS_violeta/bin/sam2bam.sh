#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 2
##SBATCH --mem=80gb
#SBATCH -t 100:00:00
#SBATCH -J sam2bam
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/sam2bam.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/sam2bam.out

module load bioinfo-tools
module load samtools/1.14
module load FastQC/0.11.9
module load MultiQC

working_dir=/proj/naiss2023-23-55/GBS_violeta/aligned

cd $working_dir

#for file in *.sam; do
#    samtools view -S -b $file > ${file/%sam/bam}
#done

#for file in *.bam; do
#    samtools sort -m 768M $file > ${file/%bam/sorted.bam}
#done

#for i in *sorted.bam
#do
#  echo "Indexing: "$i
#  samtools index $i $i.bai
#done

#As this is the genomic part, we only want to have the fragments that map uniquely

#for i in *sorted.bam; do
#  sample=$(basename $i .sorted.bam)
#  samtools view -h $i | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > $sample.unique.bam
#done

for file in *.unique.bam; do
    samtools sort -m 768M $file > ${file/%bam/sorted.bam}
done

for i in *unique.sorted.bam
do
  samtools index $i $i.bai
done

#for file in *.sorted.bam; do
#  samtools depth $file | awk '{sum+=$3;cnt++} END {print sum/cnt" "sum}'
#done

echo "finished SAMBAM"
