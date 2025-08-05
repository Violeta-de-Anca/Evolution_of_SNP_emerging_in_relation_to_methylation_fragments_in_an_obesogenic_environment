#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core -n 4
#SBATCH -t 10-00:00:00
##SBATCH -t 1:00:00
#SBATCH -J postgatk
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/postcnv_overnutrition_mode.gatk.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/postcnv_overnutrition_mode.gatk.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK/4.3.0.0
module load conda/latest
source activate gatk
source $GATK_HOME/gatk-completion.sh

output_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK
bins_1000=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/bin_1000_GAKT.interval_list
dictionary=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/mm39.dict
ref=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/mm39.fa

sample=$1
input_folder=${sample%/*}
echo $input_folder
name_sample=${sample##*/}
echo $name_sample
cleaned=${name_sample#C13}
cleaned_1=${cleaned%_Mouse.unique.sorted.bam.HDF5}

#run it in cohort mode with the control samples
gatk PostprocessGermlineCNVCalls --calls-shard-path $output_folder/overnutrition_${cleaned_1}_cnvs-calls \
        --contig-ploidy-calls $output_folder/overnutrition_${cleaned_1}-calls \
        --model-shard-path $output_folder/control_model_cnv_caller_GATK-model \
        --output-denoised-copy-ratios $output_folder/denoised_${cleaned_1}_cnv.vcf \
        --output-genotyped-intervals $output_folder/intervals_${cleaned_1}_cnv.vcf \
        --output-genotyped-segments $output_folder/segment_${cleaned_1}_cnv.vcf \
        --reference $ref \
        --allosomal-contig chrX \
        --allosomal-contig chrY \
        --sequence-dictionary $dictionary

