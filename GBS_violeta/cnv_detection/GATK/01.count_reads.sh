#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J countgatk
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/counts_in_bins.gatk.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/counts_in_bins.gatk.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK/4.3.0.0
module load conda/latest
source activate gatk

sample=$1
input_folder=${sample%/*}
echo $input_folder
name_sample=${sample##*/}
echo $name_sample
name=${name_sample%.*}
echo $name

output_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK
bins_1000=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/bin_1000_GAKT.interval_list

#the bins done by GATK need to have the .interval_list extension!!!!!

#done per sample
#gatk CollectReadCounts -imr OVERLAPPING_ONLY --input ${sample} --intervals ${bins_1000} --output ${output_folder}/${name_sample}.HDF5

#for the random bams for the control
#there is an error with the names, so I'm going to try and do a simplier name input just in case
#gatk CollectReadCounts -imr OVERLAPPING_ONLY --input ${output_folder}/${sample} --intervals ${bins_1000} --output ${output_folder}/${name}.HDF5

cd $output_folder
gatk CollectReadCounts -imr OVERLAPPING_ONLY --input ${sample} --intervals ${bins_1000} --output ${name}.HDF5
