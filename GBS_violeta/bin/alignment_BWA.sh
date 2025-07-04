#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p node -C mem256GB
#SBATCH -t 100:00:00
#SBATCH -J alignment
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/alignment.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/alignment.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load bwa
#module load bwa-mem2

working_dir=/proj/naiss2023-23-55/GBS_violeta
input_prefix=$working_dir/trimmed
reference_genome=/proj/ancestry_medi_indiv/reference_genome/uscs_mus_musculus_mm39/mm39.fa.gz
output_prefix=$working_dir/aligned
mkdir -p $output_prefix

cd $input_prefix
for i in $(ls *1.trimmed.fq.gz | cut -d "." -f 1); do
	input1=$i.1.trimmed.fq.gz
	input2=$i.2.trimmed.fq.gz
	output=$output_prefix/$i.sam
	bwa mem -t 10 -M -R "@RG\tID:GBSICR\tSM:$input1\tPL:ILLUMINA" $reference_genome $input1 $input2 > $output
done

echo "finished alignment"
