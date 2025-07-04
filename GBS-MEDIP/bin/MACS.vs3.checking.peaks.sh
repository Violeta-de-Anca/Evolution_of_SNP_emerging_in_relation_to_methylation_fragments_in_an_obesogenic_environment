#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 4
#SBATCH --mem=80gb
#SBATCH -t 4-100:00:00
#SBATCH -J MACS3
#SBATCH --error /proj/ancestry_medi_indiv/MACS3.merged.verification.peaks.MQ20.err
#SBATCH --output /proj/ancestry_medi_indiv/MACS3.merged.verification.peaks.MQ20.out

module load bioinfo-tools
module load samtools/1.14
module load MACS

working_dir_GBS=/proj/naiss2023-23-55/GBS_violeta
input_prefix_GBS=$working_dir_GBS/aligned
working_dir_MeDIP=/proj/ancestry_medi_indiv
input_MEDIP=$working_dir_GBS/gbs-medip
output_dir=$working_dir_GBS/verification.peaks.MACS3
reference_genome=/proj/naiss2023-23-55/GBS_violeta/reference_genomes/mus_musculus/uscs_ref/mm39.fa
mkdir -p $output_dir

#for i in $working_dir_GBS/gbs-medip/*unique.MQ20.bam ;
#do
#	sample=$(basename $i .unique.MQ20.bam)
#	macs3 callpeak -t $i \
 #       -c $input_prefix_GBS/$sample.unique.sorted.bam \
   #     -f BAMPE \
  #      -B \
 #       -n $sample.MQ20 \
#        --trackline -q 0.05 \
#        --outdir $output_dir
#done

macs3 callpeak -t $input_MEDIP/merged.GBS_MeDIP.MQ20.sorted.bam \
	-c $input_prefix_GBS/merged.GBS.uniq.sorted.bam \
	-f BAMPE \
	-B \
	-n merged.peaks.GBSandGBS-MEDIP.MQ10.ICR.strain \
	--trackline -q 0.05 \
	--outdir $output_dir
