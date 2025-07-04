#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 4
#SBATCH --mem=80gb
#SBATCH -t 4-100:00:00
#SBATCH -J ensembl
##SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/joinsnp.nonfilter.gatk.err
##SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/joinsnp.nonfilter.gatk.out
##SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/recalibration.uniq.gatk.err
##SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/recalibration.uniq.gatk.out
##SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/filter.raw.gatk.err
##SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/filter.raw.gatk.out
##SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/actual.recalibration.gatk.err
##SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/actual.recalibration.gatk.out
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/ensembl.recalibration.gatk.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/ensembl.recalibration.gatk.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK/4.2.0.0
module load bcftools/1.17
module load QualiMap/2.2.1
module unload java

working_dir=/proj/naiss2023-23-55/GBS_violeta
input_prefix=$working_dir/aligned
db_knownsites=$working_dir/snp_indel_databases
recalibration_output=$working_dir/recalibration
reference_genome=/proj/naiss2023-23-55/GBS_violeta/reference_genomes/mus_musculus/uscs_ref/mm39.fa
mkdir -p $recalibration_output

#First we need to do a raw file with all the variations
#cd $input_prefix
#There is no mapping quality as we only have unique mappers
#for SAMPLE in $(ls *unique.sorted.bam | cut -d "." -f 1); do
#	gatk HaplotypeCaller -OVI true --emit-ref-confidence GVCF \
#		--annotation FisherStrand -A QualByDepth -A DepthPerAlleleBySample \
#		-G StandardAnnotation --input $SAMPLE.sorted.bam \
#		--output $recalibration_output/$SAMPLE.sort.dup.raw.hc.g.vcf.gz \
#		--reference $reference_genome
#done

# Merge all the vcf files so we can do stats on them to see the filters we need to apply and the thresholds

#cd $recalibration_output

#gatk CombineGVCFs -OVI true --output combine.nonfiltered.g.vcf.gz --reference $reference_genome \
#  --variant  C13F1_10_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F1_11_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant  C13F1_1_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F1_2_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant  C13F1_3_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F1_4_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant  C13F1_5_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F1_6_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F1_7_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F1_8_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F1_9_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F2_12_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F2_13_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F2_14_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F2_15_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F2_16_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F2_17_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F2_18_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F2_19_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F2_1_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F2_20_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F2_21_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F2_23_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F2_24_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F2_2_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F2_3_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F2_4_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F2_5_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F2_6_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F2_7_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F2_8_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F2_9_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F3_10_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F3_11_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F3_12_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F3_13_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F3_14_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F3_15_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F3_16_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F3_17_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F3_18_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F3_19_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F3_1_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F3_20_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F3_21_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F3_22_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F3_24_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F3_25_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F3_27_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F3_2_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F3_3_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F3_4_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F3_5_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F3_7_Mouse.sort.dup.raw.hc.g.vcf.gz \
#  --variant C13F3_8_Mouse.sort.dup.raw.hc.g.vcf.gz --variant C13F3_9_Mouse.sort.dup.raw.hc.g.vcf.gz \

# Now from the combination we got the thresholds for the depth, min: 101.6 and max: 1644 we are going to do this in the combination file
#gatk VariantFiltration -OVI true \
#  --variant combine.nonfiltered.g.vcf.gz \
#  --output combine.raw.filtered.g.vcf.gz \
#  --filter-name minDP --filter-expression "DP > 101" \
#  --filter-name maxDP --filter-expression "DP < 1644"

# Now that we have the quality SNPs, letÃ's do the recalibration
cd $input_prefix

for SAMPLE in $(ls *unique.sorted.bam ); do
	gatk BaseRecalibrator --input $SAMPLE \
  		--reference $reference_genome --known-sites $db_knownsites/mus_musculus.uscs.ensembl.vcf \
  		--output $recalibration_output/$SAMPLE.sort.dup.recal.table
	gatk ApplyBQSR --bqsr-recal-file $recalibration_output/$SAMPLE.sort.dup.recal.table \
  		--input $SAMPLE --output $recalibration_output/$SAMPLE.sorted.unique.recalibrated.bam
done
