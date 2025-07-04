#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 1
##SBATCH --mem=80gb
#SBATCH -t 4-100:00:00
#SBATCH -J QCvarifilt
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/QC_variant_calling_raw.gatk.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/QC_variant_calling_raw.gatk.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load bcftools/1.17

working_dir=/proj/naiss2023-23-55/GBS_violeta
input_prefix=$working_dir/aligned
db_knownsites=$working_dir/snp_indel_databases
recalibration_output=$working_dir/recalibration
reference_genome=/proj/naiss2023-23-55/GBS_violeta/reference_genomes/mus_musculus/uscs_ref/mm39.fa
variant_output=$working_dir/variant_output

for i in $(ls $variant_output/*unique.recalibrated.bam.non.filtered.vcf ); do
	sample=$(basename $i _Mouse.unique.sorted.bam.sorted.unique.recalibrated.bam.non.filtered.vcf)
	bcftools stats $i > $sample.bcftools.stats
done

bcftools stats -d *unique.recalibrated.bam.non.filtered.vcf > merged.bcftools.stats
#| grep ^SN > $sample.resumen.stats
