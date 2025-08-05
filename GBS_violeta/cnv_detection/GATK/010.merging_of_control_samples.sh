#!/bin/bash -l
#SBATCH -A uppmax2025-2-151
#SBATCH -p node -n 1
#SBATCH -t 10-000:00:00
#SBATCH -J merge
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/control_merging_cnv.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/control_merging_cnv.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14

output_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/aligned
folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK

samtools merge -o $folder/control_merge_unique.bam $output_folder/C13F1_11_Mouse.unique.sorted.bam \
        $output_folder/C13F1_1_Mouse.unique.sorted.bam \
        $output_folder/C13F1_2_Mouse.unique.sorted.bam \
        $output_folder/C13F1_3_Mouse.unique.sorted.bam \
        $output_folder/C13F1_4_Mouse.unique.sorted.bam \
        $output_folder/C13F1_5_Mouse.unique.sorted.bam \
        $output_folder/C13F1_6_Mouse.unique.sorted.bam \
        $output_folder/C13F1_7_Mouse.unique.sorted.bam \
        $output_folder/C13F1_8_Mouse.unique.sorted.bam \
        $output_folder/C13F1_9_Mouse.unique.sorted.bam \
        $output_folder/C13F2_2_Mouse.unique.sorted.bam \
        $output_folder/C13F2_3_Mouse.unique.sorted.bam \
        $output_folder/C13F2_4_Mouse.unique.sorted.bam \
        $output_folder/C13F2_5_Mouse.unique.sorted.bam \
        $output_folder/C13F2_6_Mouse.unique.sorted.bam \
        $output_folder/C13F2_7_Mouse.unique.sorted.bam \
        $output_folder/C13F2_8_Mouse.unique.sorted.bam \
        $output_folder/C13F2_9_Mouse.unique.sorted.bam \
        $output_folder/C13F3_10_Mouse.unique.sorted.bam \
        $output_folder/C13F3_11_Mouse.unique.sorted.bam \
        $output_folder/C13F3_12_Mouse.unique.sorted.bam \
        $output_folder/C13F3_13_Mouse.unique.sorted.bam \
        $output_folder/C13F3_1_Mouse.unique.sorted.bam \
        $output_folder/C13F3_2_Mouse.unique.sorted.bam \
        $output_folder/C13F3_3_Mouse.unique.sorted.bam \
        $output_folder/C13F3_4_Mouse.unique.sorted.bam \
        $output_folder/C13F3_5_Mouse.unique.sorted.bam \
        $output_folder/C13F3_7_Mouse.unique.sorted.bam \
        $output_folder/C13F3_8_Mouse.unique.sorted.bam \
        $output_folder/C13F3_9_Mouse.unique.sorted.bam

# now also we need to change the header cause there is different group names due to the merging
#get the header!
samtools view -H control_merge_unique.bam > header_merged_controls.txt

#chage the read group
sed 's/SM:[^[:space:]]*/SM:control_sample/g' header_merged_controls.txt > header.controls.txt

#now update the header
samtools reheader header.controls.txt control_merge_unique.bam > control_merge_unique_fixed.bam

#and index the new bam file
samtools index control_merge_unique_fixed.bam
