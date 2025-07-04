#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 1
#SBATCH -t 100:00:00
#SBATCH -J QC.bsqr
##SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/QC.nofiltered.gatk.err
##SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/QC.nofiltered.gatk.out
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/QC.bsqr.gatk.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/QC.bsqr.gatk.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK/4.2.0.0
module load bcftools/1.17
module load QualiMap/2.2.1
module load vcftools
module unload java
module load MultiQC

working_dir=/proj/naiss2023-23-55/GBS_violeta
recalibration_output=$working_dir/recalibration
reference_genome=/proj/naiss2023-23-55/GBS_violeta/reference_genomes/mus_musculus/uscs_ref/mm39.fa
mkdir -p $recalibration_output
qc_dir=$working_dir/quality_control
fastqc_output_dir=$qc_dir/fastqc_output
vcf_output=$fastqc_output_dir/snp_nofiltered

#cd $recalibration_output
#vcftools --gzvcf combine.nonfiltered.g.vcf.gz --site-depth --out site.depth.combine.nonfiltered
#vcftools --gzvcf combine.nonfiltered.g.vcf.gz --site-quality --out site.quality.combine.nonfiltered

#module load R
#module load R_packages/4.2.1

#R --no-save --quiet < $working_dir/bin/median.depth.R
#R --no-save --quiet < quality.distribution.vcf.R

#After recalibration QC
cd $recalibration_output
#for SAMPLE in $(ls *.sorted.unique.recalibrated.bam ); do
#	gatk BaseRecalibrator --input $SAMPLE \
# 		--reference $reference_genome --known-sites /proj/naiss2023-23-55/GBS_violeta/snp_indel_databases/mus_musculus.uscs.ensembl.vcf \
#	  	--output $SAMPLE.after.table
#done

#for SAMPLE in $(ls *.sort.dup.recal.table ); do
#	i=$(basename $SAMPLE .sort.dup.recal.table)
#	gatk AnalyzeCovariates -before $recalibration_output/$SAMPLE \
#		-after $i.sorted.unique.recalibrated.bam.after.table \
#		-plots $vcf_output/$SAMPLE.plots.bsqr \
#		-csv $vcf_output/$SAMPLE.csv.bsqr
#done

cd $vcf_output
multiqc .
