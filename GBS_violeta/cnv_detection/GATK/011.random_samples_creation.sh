#!/bin/bash -l
#SBATCH -A uppmax2025-2-151
#SBATCH -p core -n 1
#SBATCH -t 10-000:00:00
#SBATCH -J merge
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/control_intermediate_random_cnv.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/control_merging_cnv.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

input_n_seq=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/quality_control/multiqc_output/aligned.unique_GBS.QC_data
module load bioinfo-tools samtools/1.20
input_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK
folder=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref

#rm $input_folder/number_total_seq_controls.txt

#get the number of sequences per individual, this is the total (QC passed and NOT passed)
#for i in $input_folder/*_Mouse.unique.sorted.bam.flagstats; do
#	head -n 1 $i | cut -f 1 -d " " >> $input_folder/number_total_seq_controls.txt
#done

#but samtools merge only keep the passed QC reads so we are going to extract the properly paired
#for i in $input_folder/*_Mouse.unique.sorted.bam.flagstats; do
#	sed -n '12p' $i | cut -f 1 -d " " >> $input_folder/number_total_seq_controls.txt
#done

#let's try with the mapped reads
#for i in $input_folder/*_Mouse.unique.sorted.bam.flagstats; do
#	sed -n '7p' $i | cut -f 1 -d " " >> $input_folder/number_total_seq_controls.txt
#done

#number of sequences in the merged bam file
#samtools flagstat $input_folder/control_merge_unique.bam > $input_folder/number_seq_merged_control.txt

#do the % per row to input into samtools
#head -n 1 $input_folder/number_seq_merged_control.txt | cut -f 1 -d " " > $input_folder/total_n_seq_merged.txt
#rm $input_folder/number_seq_merged_control.txt

#with that list iterate per line and extract randomly that number from the merged bam file
#awk -v total=$(cat $input_folder/total_n_seq_merged.txt) '{printf "%.2f\n", $1 / total}' $input_folder/number_total_seq_controls.txt > $input_folder/normalized_by_total.txt

#see if the percentages are good, they should sum up 1
#echo 'awk '{sum += $1} END {printf "%.2f\n", sum}' $input_folder/normalized_by_total.txt'

#now create the random samples with the % calculated before
	#so we need a unique sample name per random bam otherwise GATK gives an error
#i=1
#while read -r line; do
#	samtools view -s $line -b $input_folder/control_merge_unique_fixed.bam > $input_folder/random_control_${i}.sam
#	samtools view -bS $input_folder/random_control_${i}.sam > $input_folder/random_control_${i}.bam
#	samtools view -H $input_folder/random_control_${i}.bam | sed "s/SM:[^[:space:]]*/SM:control_sample_${i}/g" > tem_header.txt
#	samtools reheader tem_header.txt $input_folder/random_control_${i}.bam > $input_folder/random_control_good_header_${i}.bam
#	mv $input_folder/random_control_good_header_${i}.bam $input_folder/random_control_${i}.bam
#	((i++))
#done < $input_folder/normalized_by_total.txt

#that solution is not enough, so we need to change all the read groups for all the bam files
i=1
while read -r line; do
	bam_out=$input_folder/random_control_${i}.bam
	sample_name=control_sample_${i}
	rgid=rg${i}
	# Subsample from merged BAM (reference required for correct header)
	samtools view -s $line -b -T $folder/mm39.fa  $input_folder/control_merge_unique_fixed.bam > $bam_out
	samtools view -H $bam_out > $input_folder/temp
	grep -v "^@PG" $input_folder/temp | grep -v "^@RG" > $input_folder/temp.1
	echo -e "@RG\tID:${i}\tSM:control_${i}\tLB:lib1\tPL:ILLUMINA\tPU:unit1" >> $input_folder/temp.1
	samtools reheader $input_folder/temp.1 $bam_out > $input_folder/random_control_x_${i}.bam
	mv $input_folder/random_control_x_${i}.bam $input_folder/random_control_${i}.bam
	samtools index $bam_out
	((i++))
done < $input_folder/normalized_by_total.txt

# we need to create the bai!
#for i in $input_folder/random_control*bam; do
#	samtools index $i
#done
