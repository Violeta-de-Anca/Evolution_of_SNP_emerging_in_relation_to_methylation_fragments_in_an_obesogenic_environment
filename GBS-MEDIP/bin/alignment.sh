#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p node -C mem256GB
#SBATCH -t 100:00:00
#SBATCH -J alignment
#SBATCH --error /proj/ancestry_medi_indiv/alignment.err
#SBATCH --output /proj/ancestry_medi_indiv/alignment.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools

module load Stacks/2.62
module load cutadapt/4.0
module load MultiQC
module load samtools/1.14
module load bowtie2/2.3.5.1
module load FastQC/0.11.9

working_dir=/proj/ancestry_medi_indiv
input_prefix=$working_dir/trimmed
reference_genome=/proj/ancestry_medi_indiv/reference_genome/uscs_mus_musculus_mm39/mm39.uscs
# if the output folder has not been created
mkdir -p $working_dir/aligned

################## set the prefix for output files##########################
output_prefix=$working_dir/aligned

#############################################################################

cd $input_prefix

for i in $(ls *1.trimmed.fq.gz | cut -d "." -f 1); do
    # set the input and output files
    input1=$i.1.trimmed.fq.gz
    input2=$i.2.trimmed.fq.gz
    input_rem1=$i.rem.1.trimmed.fq.gz
    input_rem2=$i.rem.2.trimmed.fq.gz
    output=$output_prefix/$i.sam
    bowtie2 --threads 10 --very-sensitive-local -x $reference_genome -1 $input_prefix/$input1 -2 $input_prefix/$input2 -U $input_prefix/$input_rem1 -U $input_prefix/$input_rem2 -S $output

done

echo finished Alignment
