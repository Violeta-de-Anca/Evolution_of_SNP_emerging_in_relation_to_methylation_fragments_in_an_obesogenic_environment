#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 4 # -p node -C mem256GB
#SBATCH --mem=80gb
#SBATCH -t 100:00:00
#SBATCH -J demulti
#SBATCH --error /proj/ancestry_medi_indiv/demultiplexing.err
#SBATCH --output /proj/ancestry_medi_indiv/demultiplexing.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools

module load Stacks/2.62
module load cutadapt/4.0
module load MultiQC
module load samtools/1.14
module load bowtie2/2.3.5.1
module load FastQC/0.11.9

#In the working directory folder should be a folder called barcodes and inside a tab separated file with barcode sequence and name
working_dir=/proj/ancestry_medi_indiv

#Define where your input fastq file directory is located
input_fastq_dir=/proj/ancestry_medi_indiv/samples/Undetermined

#Defining your barcode file location
barcode_file=/proj/ancestry_medi_indiv/barcodes/barcode.txt

#defining your output directory
demultiplexed_dir=$working_dir/demultiplexed

#navegate to the script directory:
cd $working_dir
#create the folder in your script directory, if they do not exist
mkdir -p aligned
mkdir -p demultiplexed

#run the demultiplexing command
process_radtags --threads 4 -P -p $input_fastq_dir  -o $working_dir/demultiplexed -b $barcode_file -e pstI -r -q -c --len-limit 20

# --discards capture discarded reads to a file.
echo finished demultiplexing
