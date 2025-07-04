#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 1
##SBATCH --mem=80gb
#SBATCH -t 100:00:00
#SBATCH -J QC
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/quality.control.bam.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/quality.control.bam.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
#module load Stacks/2.62
#module load cutadapt/4.0
module load MultiQC
module load samtools/1.14
#module load bowtie2/2.3.5.1
module load FastQC/0.11.9
#module load picard
module load bcftools/1.17
module load QualiMap/2.2.1
#module unload java

# Define working directory
working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta
trimmed_input=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/trimmed
aligned_input=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/aligned
snp_nofiltered=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/recalibration
snp_call=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/variant_output
gbs_medip_input=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/gbs-medip

# Define output directory
qc_dir=$working_dir/quality_control
fastqc_output_dir=$qc_dir/fastqc_output
demultiplexed_output=$fastqc_output_dir/demultiplexed
trimmed_output=$fastqc_output_dir/trimmed
reference_genome=/proj/ancestry_medi_ind/reference_genome/GCF_000001635.27_GRCm39_genomic.fna
multiqc_output_dir=$qc_dir/multiqc_output
depth_output_dir=$fastqc_output_dir/depth
aligned_output=$fastqc_output_dir/aligned
vcf_output=$fastqc_output_dir/snp_nofiltered

#Create the quality control directory
mkdir -p $qc_dir
# Create the FastQC output directory
mkdir -p $fastqc_output_dir
# Create the MultiQC output directory
mkdir -p $multiqc_output_dir
# Create the depth output directory
mkdir -p $depth_output_dir
# Create the demultiplexed QC output directory
mkdir -p $demultiplexed_output
# Create the trimmed QC output directory
mkdir -p $trimmed_output
# Create the MultiQC output directory
mkdir -p $aligned_output
# Create the variant output directory
mkdir -p $vcf_output

# Run FastQC on the demultiplexed files
#file_list=($samples_dir/*.fq.gz)
# Run FastQC on each file in the list
#for file in "${file_list[@]}"; do
#  fastqc -o $demultiplexed_output $file
#done

# Run FastQC on the trimmed files
#file_list=($trimmed_input/*.fq.gz)
# Run FastQC on each file in the list
#for file in "${file_list[@]}"; do
#  fastqc -o $trimmed_output $file
#done

# Run samtools flagstat on the aligned files
#file_list=($gbs_medip_input/*unique.MQ20.bam)
#for file in "${file_list[@]}"; do
#	#unset DISPLAY
#	name=${file##*/}
#	x=${name%.unique.MQ20.bam}
#	samtools flagstat $file > $aligned_output/$x.MQ10.GBS.MEDIP.flagstat
#	samtools idxstats $file > $aligned_output/$x.MQ10.GBS.MEDIP.idxstats
#	qualimap bamqc -bam $file -c -outdir $aligned_output -outfile $x -outformat PDF:HTML -ip
#done

#for i in $aligned_input/*.unique.sorted.bam
#do
#    echo "Running QualiMap on $i ..."
#    bash qualimap.sh $i
#done

#Now do QC on the vcf file with the non-filtered SNPs:
#cd $snp_call
#for SAMPLE in $(ls *.recalibrated.bam.non.filtered.vcf); do
#	bcftools stats $SAMPLE > $vcf_output/$SAMPLE.bcftools
#	grep "^SN" $vcf_output/$SAMPLE.bcftools > $vcf_output/$SAMPLE.summary.stats
#done

#Now do QC on the vcf file with the filtered SNPs:
#cd $snp_call
#for SAMPLE in $(ls *sorted.recal.filtered.vcf); do
#       bcftools stats $SAMPLE > $vcf_output/$SAMPLE.bcftools
#       grep "^SN" $vcf_output/$SAMPLE.bcftools > $vcf_output/$SAMPLE.summary.stats
#done

#Do last QC of genotyped filtered per invididual #
#cd $snp_call
#for SAMPLE in $(ls *Mouse.genotypes_filtered_gatk.vcf); do
#       bcftools stats $SAMPLE > $vcf_output/$SAMPLE.bcftools
#       grep "^SN" $vcf_output/$SAMPLE.bcftools > $vcf_output/$SAMPLE.summary.stats
#done

# number of mapped reads #
#cd $gbs_medip_input
#for SAMPLE in $(ls *Mouse.unique.MQ20.bam); do
#	echo $SAMPLE >> Num_reads_per_sample_GBS-MEDIP
#	samtools view -f 0x02 $SAMPLE | wc -l >> Num_reads_per_sample_GBS-MEDIP
#done

#cd $aligned_input
#for SAMPLE in $(ls *Mouse.unique.sorted.bam); do
#	echo $SAMPLE >> Num_reads_per_sample_GBS
#	samtools view -f 0x02 $SAMPLE | wc -l >> Num_reads_per_sample_GBS
#done

#bcftools stats combine.nonfiltered.g.vcf.gz > $vcf_output/combine.nonfiltered.bcftools
#grep "^SN" $vcf_output/combine.nonfiltered.bcftools > $vcf_output/combine.nonfiltered.summary.stats
#plot-vcfstats -p $vcf_output $vcf_output/combine.nonfiltered.bcftools

# Run MultiQC on the demultiplexed FastQC output
#multiqc $demultiplexed_output/* -o $multiqc_output_dir -n demultiplexed_GBS.QC

# Run MultiQC on the trimmed fastQC output
#multiqc $trimmed_output/* -o $multiqc_output_dir -n trimmed_GBS.QC

# Run MultiQC on the aligned fastQC output
#cd $aligned_output
#multiqc -o $multiqc_output_dir -n aligned.unique_GBS.QC .

#echo finished fastQC and MultiQC
