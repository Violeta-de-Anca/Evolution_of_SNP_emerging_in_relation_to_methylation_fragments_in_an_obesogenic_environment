#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 1
#SBATCH -t 100:00:00
#SBATCH -J filterrecal
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/filter_recalibration.gatk.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/filter_recalibration.gatk.out
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

cd $recalibration_output

#so we want to have at least 20X of coverage and a good quality map to be able to have confidence
for SAMPLE in $(ls *.sort.dup.raw.hc.g.vcf.gz | cut -d "." -f 1); do
        gatk VariantFiltration -OVI true \
                --variant $SAMPLE.sort.dup.raw.hc.g.vcf.gz \
                --output $recalibration_output/$SAMPLE.sort.dup.raw.hc.filtered.g.vcf.gz \
                --filter-expression "MMQ > 20.0" --filter-name "medianmappingqualitymorethan20"  --filter-expression "DP > 20" --filter-name "readdepthmorethan20"
done
