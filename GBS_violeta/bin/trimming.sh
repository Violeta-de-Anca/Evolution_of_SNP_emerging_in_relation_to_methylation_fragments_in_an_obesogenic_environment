#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 4 # -p node -C mem256GB
#SBATCH --mem=80gb
#SBATCH -t 100:00:00
#SBATCH -J trimm
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/Trimming.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/Trimming.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools

module load Stacks/2.62
module load cutadapt/4.0

# Define working directory
working_dir=/proj/naiss2023-23-55/GBS_violeta

# set input directory for reads
reads_dir=/proj/ancestry_medi_indiv_comp/private/GBSMouse/Stacks/samples

# set output directory for trimmed reads
trimmed_dir=$working_dir/trimmed

# create output directory if it doesn't exist
mkdir -p $trimmed_dir

# process paired-end reads
for i in $reads_dir/*.1.fq.gz ; do
    sample=$(basename $i .1.fq.gz)
    cutadapt -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT \
    -g AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG \
    -A CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT \
    -G AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -Q 35,0 --quality-base=30 -m 30 \
    -o $trimmed_dir/$sample.1.trimmed.fq.gz \
    -p $trimmed_dir/$sample.2.trimmed.fq.gz \
    $reads_dir/$sample.1.fq.gz \
    $reads_dir/$sample.2.fq.gz
done

echo finished trimming

