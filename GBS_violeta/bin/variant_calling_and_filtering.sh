#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 4
##SBATCH --mem=80gb
#SBATCH -t 10-00:00:00
#SBATCH -J refiltered
##SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/merged.filtered.gatk.err
##SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/merged.filtered.gatk.out
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/re.filtered.gatk.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/re.filtered.gatk.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK/4.2.0.0
module unload java

working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta
input_prefix=$working_dir/aligned
db_knownsites=$working_dir/snp_indel_databases
recalibration_output=$working_dir/recalibration
reference_genome=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/mm39.fa
variant_output=$working_dir/variant_output
mkdir -p $variant_output

###variant calling per individual###

#cd $recalibration_output
#for SAMPLE in $(ls *.unique.sorted.bam.sorted.unique.recalibrated.bam); do
#	gatk HaplotypeCaller -OVI true --emit-ref-confidence GVCF --annotation FisherStrand -A AlleleFraction \
#	   -A AS_QualByDepth \
#	   -A Coverage -A QualByDepth -A DepthPerAlleleBySample \
#	   --dbsnp /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/snp_indel_databases/mus_musculus.uscs.ensembl.vcf \
#           -G StandardAnnotation --input $SAMPLE \
#           -O $SAMPLE.non.filtered.vcf \
#           --reference $reference_genome
#done

### Merging of all the individuals ####
#cd $variant_output
#gatk CombineGVCFs -OVI true --annotation FisherStrand -A AlleleFraction \
#          -A AS_QualByDepth \
#          -A Coverage -A QualByDepth -A DepthPerAlleleBySample \
#          --dbsnp /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/GBS_violeta/snp_indel_databases/mus_musculus.uscs.ensembl.vcf \
#	  --output merged_ICR_variant_calling_gatk.vcf --reference $reference_genome \
#	  -V variant_calling_ICR_GBS.list

#cd $variant_output
### Filtering per individual ####
#for SAMPLE in $(ls *sorted.unique.recalibrated.bam.non.filtered.vcf); do
#        x=${SAMPLE%.unique.sorted.bam.sorted.unique.recalibrated.bam.non.filtered.vcf}
#	gatk VariantFiltration -OVI true \
#  	   --variant $SAMPLE \
#           --filter-name "FisherStrand" --filter-expression "FS > 60.0" \
#           --filter-name "QualByDepth" --filter-expression "QD < 2.0" \
#           --filter-name "MappingQuality" --filter-expression "MQ < 40.0" \
#	   --filter-name "allelicDepth" --filter-expression "AD > 10" --cluster-size 3 --cluster-window-size 50 \
#	   --filter-name "mindepth" --filter-expression "DP>20" --filter-name "maxdepth" --filter-expression "DP<1644" \
#           -O $x.uniq.sorted.recal.filtered.vcf \
#           --reference $reference_genome
#done

#### Merging with the filtered SNPs vcf file ####
#cd $variant_output
#gatk CombineGVCFs -OVI true --annotation FisherStrand -A AlleleFraction \
#         -A AS_QualByDepth \
#         -A Coverage -A QualByDepth -A DepthPerAlleleBySample \
#         --dbsnp /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/snp_indel_databases/mus_musculus.uscs.ensembl.vcf \
#         --output merged_ICR_filtered_variant_calling_gatk.vcf --reference $reference_genome \
#         -V filtered.vcf.list

### Select only the ones that passed the filters per ind for QC ###
#cd $variant_output
#for SAMPLE in $(ls *uniq.sorted.recal.filtered.vcf); do
#	x=${SAMPLE%.uniq.sorted.recal.filtered.vcf}
#	gatk GenotypeGVCFs -OVI true --reference $reference_genome -V $SAMPLE \
#	  --output $x.genotypes_filtered_gatk.vcf \
#	  --dbsnp /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/snp_indel_databases/mus_musculus.uscs.ensembl.vcf
#done


## Genotyping the merged vcf file ###
#cd $variant_output
#gatk GenotypeGVCFs -OVI true --reference $reference_genome -V merged_ICR_filtered_variant_calling_gatk.vcf \
#  --output ICR_genotypes_filtered_gatk.vcf \
#  --dbsnp /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/snp_indel_databases/mus_musculus.uscs.ensembl.vcf

#cd $variant_output
### Filtering per individual now for analizing dynamics in the methylated regions ####
cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/rerun_filters_vcf
for SAMPLE in $(ls *); do
       x=${SAMPLE%_Mouse.unique.sorted.bam.sorted.unique.recalibrated.bam.non.filtered.vcf}
       gatk VariantFiltration -OVI true \
          --variant $SAMPLE \
          --filter-name "FisherStrand" --filter-expression "FS > 60.0" \
          --filter-name "QualByDepth" --filter-expression "QD < 2.0" \
          --filter-name "MappingQuality" --filter-expression "MQ < 40.0" \
          --filter-name "allelicDepth" --filter-expression "AD > 10" \
          --filter-name "mindepth" --filter-expression "DP>20" --filter-name "maxdepth" --filter-expression "DP<1644" \
          -O $x.uniq.sorted.recal.filt.nocluster.vcf \
          --reference $reference_genome
done

for SAMPLE in $(ls *uniq.sorted.recal.filt.nocluster.vcf); do
       x=${SAMPLE%.uniq.sorted.recal.filt.nocluster.vcf}
       gatk GenotypeGVCFs -OVI true --all-sites --reference $reference_genome -V $SAMPLE \
         --output $x.genotypes_refiltered_gatk.vcf \
         --dbsnp /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/snp_indel_databases/mus_musculus.uscs.ensembl.vcf
done

#now do the fathers #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/rerun_filters_vcf
#for SAMPLE in $(ls C13F2*); do
#       x=${SAMPLE%_Mouse.unique.sorted.bam.sorted.unique.recalibrated.bam.non.filtered.vcf}
#       gatk VariantFiltration -OVI true \
#          --variant $SAMPLE \
#           --filter-name "FisherStrand" --filter-expression "FS > 60.0" \
#           --filter-name "QualByDepth" --filter-expression "QD < 2.0" \
#           --filter-name "MappingQuality" --filter-expression "MQ < 40.0" \
#          --filter-name "allelicDepth" --filter-expression "AD > 10" \
#          --filter-name "mindepth" --filter-expression "DP>20" --filter-name "maxdepth" --filter-expression "DP<1644" \
#           -O $x.uniq.sorted.recal.filt.nocluster.vcf \
#           --reference $reference_genome
#done

#for SAMPLE in $(ls *uniq.sorted.recal.filt.nocluster.vcf); do
#       x=${SAMPLE%.uniq.sorted.recal.filt.nocluster.vcf}
#       gatk GenotypeGVCFs -OVI true --all-sites --reference $reference_genome -V $SAMPLE \
#         --output $x.genotypes_refiltered_gatk.vcf \
#         --dbsnp /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/snp_indel_databases/mus_musculus.uscs.ensembl.vcf
#done


## Genotyping the merged vcf file and getting all the sites to do pop gen test ###
#cd $variant_output
#gatk GenotypeGVCFs -OVI true --all-sites --reference $reference_genome -V merged_ICR_filtered_variant_calling_gatk.vcf \
#  --output ICR_all.sites.genotypes_filtered_gatk.vcf \
#  --dbsnp /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/snp_indel_databases/mus_musculus.uscs.ensembl.vcf

