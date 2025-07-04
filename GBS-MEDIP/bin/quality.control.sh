#!/bin/bash -l
#SBATCH -A naiss2023-22-848
#SBATCH -p core -n 1
# -p node -C mem256GB
#SBATCH --mem=80gb
#SBATCH -t 100:00:00
#SBATCH -J QC
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/quality.control.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/quality.control.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools

#module load Stacks/2.62
#module load cutadapt/4.0
#module load MultiQC
module load samtools/1.14
#module load bowtie2/2.3.5.1
#module load FastQC/0.11.9
#module load picard
#module load QualiMap

# Define working directory
working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP
# Define input directory
demulti_dir=$working_dir/demultiplexed
trimmed_dir=$working_dir/trimmed
align_dir=$working_dir/aligned
qc_dir=$working_dir/quality_control
fastqc_output_dir=$qc_dir/fastqc_output
reference_genome=/proj/ancestry_medi_indiv/reference_genome/GCF_000001635.27_GRCm39_genomic.fna
multiqc_output_dir=$qc_dir/multiqc_output
alignedQCoutput=$fastqc_output_dir/aligned_QC
depth_output_dir=$qc_dir/depth
mkdir -p $working_dir/quality_control
mkdir -p $depth_output_dir

# Create the FastQC output directory
mkdir -p $fastqc_output_dir
# Create the MultiQC output directory
mkdir -p $multiqc_output_dir

# Run FastQC on the demultiplexed files
# Retrieve the list of *fq.gz files
#file_list=($demulti_dir/*.fq.gz)
# Run FastQC on each file in the list
#for file in "${file_list[@]}"; do
#  fastqc -o $fastqc_output_dir $file
#done
# Create a file list of the QC of the demultiplexed fastqs
#ls -d "$fastqc_output_dir"/*. > demulti.list.txt

# Run FastQC on the trimmed files
#file_list=($trimmed_dir/*trimmed.fq.gz)
# Run FastQC on each file in the list
#for file in "${file_list[@]}"; do
#  fastqc -o "$fastqc_output_dir" "$file"
#done
# Create a file list of the QC of the trimmed fastqs
#ls -d "$fastqc_output_dir"/*.trimmed_fastqc.zip > trimmed.list.txt

# Run FastQC on the aligned files
#file_list=($align_dir/*.bam)
# Run FastQC on each file in the list
#for file in ${file_list[@]}; do
#	java -jar $PICARD_ROOT/picard.jar CollectAlignmentSummaryMetrics --INPUT $file --REFERENCE_SEQUENCE $reference_genome --OUTPUT $fastqc_output_dir/$file
#	i=${file##*/}
	#samtools flagstat $file > $alignedQCoutput/$i.flagstat
        #samtools idxstats $file > $alignedQCoutput/$i.idxstats
#done

#qualimap multi-bamqc -c -d bamlist.txt -r -outdir ../../multiqc_output/ -outfile multi_qualimap_GBS_ICR -outformat PDF

#cd  $alignedQCoutput
#multiqc -o $multiqc_output_dir -i aligned.GBS_MeDIP .
# Create a file list of the QC of the trimmed fastqs
#ls -d "$fastqc_output_dir"/*.sorted_fastqc.zip > aligned.list.txt

# Run MultiQC on the FastQC output
#multiqc $fastqc_output_dir -o $multiqc_output_dir -n demultiplexed.QC -l demulti.list.txt
#multiqc $fastqc_output_dir -o $multiqc_output_dir -n trimmed.QC -l trimmed.list.txt
#multiqc $fastqc_output_dir -o $multiqc_output_dir -n aligned.QC -l aligned.list.txt

# Run the depth per window in all the individuals

samtools depth -a -H -b /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/merged/combined_regions_sorted.bed -o /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/merged/depth.MACS3.wind.MQ10.ICR.testis -f /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/aligned/bam.mq10


#echo finished fastQC and MultiQC
