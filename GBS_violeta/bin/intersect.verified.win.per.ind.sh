#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 1
##SBATCH --mem=80gb
#SBATCH -t 4-100:00:00
#SBATCH -J intersect
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/win.intersect.MACS3.per.ind.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/win.intersect.MACS3.per.ind.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load plink
module load BEDTools
module load bcftools

working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta
veri_input=$working_dir/verification.peaks.MACS3
ind_win=$veri_input/peaks_MQ10_MACS3
inter_out=$veri_input/intersect_win_per_indv/GBS_with_MACS3.win
aligned_input=$working_dir/aligned
parents=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/no_cluster_filters
mkdir -p $inter_out

# intersect windows called invidually in GBS-MEDIP data using MACS3 with merged windows also in MACS3 #
#cd $ind_win
#file_list=($ind_win/*Mouse.MQ20_summits.bed)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%_Mouse.MQ20_summits.bed}
#	bedtools intersect -a $file -b $veri_input/veri.wind.MACS3.summit.location.MQ10.bed > $inter_out/$x.MQ10.intersect
#done

# Intersect windows called individually in MACS3 with the GBS bam files #
#cd $aligned_input
#file_list=($aligned_input/*Mouse.unique.sorted.bam)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%_Mouse.unique.sorted.bam}
#	y=${name%.unique.sorted.bam}
#        bedtools intersect -a $ind_win/$y.MQ20_peaks.narrowPeak -u -wa -b $file > $inter_out/$x.GBS.with.MACS3.win.intersect.bed
#done

# Count the number of windows and put it in a file
#file_list=($inter_out/*GBS.with.MACS3.win.intersect.bed)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.GBS.with.MACS3.win.intersect.bed}
#	echo $file >> MACS3.win.with.GBS.summary
#	wc -l $file >> MACS3.win.with.GBS.summary
#done

# Now intersect all the files per gen and group? #
#cd $inter_out
#F1 obese
#bedtools multiinter -header -i \
#	C13F2_12.GBS.with.MACS3.win.intersect.bed\
#        C13F2_13.GBS.with.MACS3.win.intersect.bed\
#        C13F2_14.GBS.with.MACS3.win.intersect.bed\
#        C13F2_15.GBS.with.MACS3.win.intersect.bed\
#        C13F2_16.GBS.with.MACS3.win.intersect.bed\
#        C13F2_17.GBS.with.MACS3.win.intersect.bed\
#        C13F2_18.GBS.with.MACS3.win.intersect.bed\
#        C13F2_19.GBS.with.MACS3.win.intersect.bed\
#	C13F2_20.GBS.with.MACS3.win.intersect.bed\
#        C13F2_21.GBS.with.MACS3.win.intersect.bed\
#        C13F2_23.GBS.with.MACS3.win.intersect.bed\
#        C13F2_24.GBS.with.MACS3.win.intersect.bed > win.MACS3.obese.F1.bed
#count windows with cover in 80% indv
#awk '$4>10' win.MACS3.obese.F1.bed > windows.MACS3_80percent.obese.F1.bed
#awk '$4>6' win.MACS3.obese.F1.bed > windows.MACS3_50percent.obese.F1.bed
#F1 Control
#bedtools multiinter -header -i \
#	C13F2_1.GBS.with.MACS3.win.intersect.bed\
#	C13F2_2.GBS.with.MACS3.win.intersect.bed\
#        C13F2_3.GBS.with.MACS3.win.intersect.bed\
#        C13F2_4.GBS.with.MACS3.win.intersect.bed\
#        C13F2_5.GBS.with.MACS3.win.intersect.bed\
#        C13F2_6.GBS.with.MACS3.win.intersect.bed\
#        C13F2_7.GBS.with.MACS3.win.intersect.bed\
#        C13F2_8.GBS.with.MACS3.win.intersect.bed\
#        C13F2_9.GBS.with.MACS3.win.intersect.bed > win.MACS3.F1.control.bed
#awk '$4>8' win.MACS3.F1.control.bed > windows.MACS3_80percent.control.F1.bed
#awk '$4>5' win.MACS3.F1.control.bed > windows.MACS3_50percent.control.F1.bed
#F2 obese
#bedtools multiinter -header -i \
#	C13F3_14.GBS.with.MACS3.win.intersect.bed\
#        C13F3_15.GBS.with.MACS3.win.intersect.bed\
#        C13F3_16.GBS.with.MACS3.win.intersect.bed\
#        C13F3_17.GBS.with.MACS3.win.intersect.bed\
#        C13F3_18.GBS.with.MACS3.win.intersect.bed\
#        C13F3_19.GBS.with.MACS3.win.intersect.bed\
#	C13F3_20.GBS.with.MACS3.win.intersect.bed\
#        C13F3_21.GBS.with.MACS3.win.intersect.bed\
#        C13F3_22.GBS.with.MACS3.win.intersect.bed\
#        C13F3_24.GBS.with.MACS3.win.intersect.bed\
#        C13F3_25.GBS.with.MACS3.win.intersect.bed\
#        C13F3_27.GBS.with.MACS3.win.intersect.bed > win.MACS3.F2.obese.bed
#awk '$4>10' win.MACS3.F2.obese.bed > windows.MACS3_80percent.obese.F2.bed
#awk '$4>6' win.MACS3.F2.obese.bed > windows.MACS3_50percent.obese.F2.bed
#F2 control
#bedtools multiinter -header -i \
#	C13F3_2.GBS.with.MACS3.win.intersect.bed\
#        C13F3_3.GBS.with.MACS3.win.intersect.bed\
#        C13F3_4.GBS.with.MACS3.win.intersect.bed\
#        C13F3_5.GBS.with.MACS3.win.intersect.bed\
#        C13F3_7.GBS.with.MACS3.win.intersect.bed\
#        C13F3_8.GBS.with.MACS3.win.intersect.bed\
#        C13F3_9.GBS.with.MACS3.win.intersect.bed\
#	C13F3_1.GBS.with.MACS3.win.intersect.bed > win.MACS3.F2.control.bed
#awk '$4>10' win.MACS3.F2.control.bed > windows.MACS3_80percent.control.F2.bed
#awk '$4>6' win.MACS3.F2.control.bed > windows.MACS3_50percent.control.F2.bed

#cd $inter_out
#bedtools multiinter -header -i C13F1_10.GBS.with.MACS3.win.intersect.bed \
#	C13F1_11.GBS.with.MACS3.win.intersect.bed \
#	C13F1_1.GBS.with.MACS3.win.intersect.bed \
#	C13F1_2.GBS.with.MACS3.win.intersect.bed\
#	C13F1_3.GBS.with.MACS3.win.intersect.bed\
#	C13F1_4.GBS.with.MACS3.win.intersect.bed\
#	C13F1_5.GBS.with.MACS3.win.intersect.bed\
#	C13F1_6.GBS.with.MACS3.win.intersect.bed\
#	C13F1_7.GBS.with.MACS3.win.intersect.bed\
#	C13F1_8.GBS.with.MACS3.win.intersect.bed\
#	C13F1_9.GBS.with.MACS3.win.intersect.bed\
#	C13F2_12.GBS.with.MACS3.win.intersect.bed\
#	C13F2_13.GBS.with.MACS3.win.intersect.bed\
#	C13F2_14.GBS.with.MACS3.win.intersect.bed\
#	C13F2_15.GBS.with.MACS3.win.intersect.bed\
#	C13F2_16.GBS.with.MACS3.win.intersect.bed\
#	C13F2_17.GBS.with.MACS3.win.intersect.bed\
#	C13F2_18.GBS.with.MACS3.win.intersect.bed\
#	C13F2_19.GBS.with.MACS3.win.intersect.bed\
#	C13F2_1.GBS.with.MACS3.win.intersect.bed\
#	C13F2_20.GBS.with.MACS3.win.intersect.bed\
#	C13F2_21.GBS.with.MACS3.win.intersect.bed\
#	C13F2_23.GBS.with.MACS3.win.intersect.bed\
#	C13F2_24.GBS.with.MACS3.win.intersect.bed\
#	C13F2_2.GBS.with.MACS3.win.intersect.bed\
#	C13F2_3.GBS.with.MACS3.win.intersect.bed\
#	C13F2_4.GBS.with.MACS3.win.intersect.bed\
#	C13F2_5.GBS.with.MACS3.win.intersect.bed\
#	C13F2_6.GBS.with.MACS3.win.intersect.bed\
#	C13F2_7.GBS.with.MACS3.win.intersect.bed\
#	C13F2_8.GBS.with.MACS3.win.intersect.bed\
#	C13F2_9.GBS.with.MACS3.win.intersect.bed\
#	C13F3_10.GBS.with.MACS3.win.intersect.bed\
#	C13F3_11.GBS.with.MACS3.win.intersect.bed\
#	C13F3_12.GBS.with.MACS3.win.intersect.bed\
#	C13F3_13.GBS.with.MACS3.win.intersect.bed\
#	C13F3_14.GBS.with.MACS3.win.intersect.bed\
#	C13F3_15.GBS.with.MACS3.win.intersect.bed\
#	C13F3_16.GBS.with.MACS3.win.intersect.bed\
#	C13F3_17.GBS.with.MACS3.win.intersect.bed\
#	C13F3_18.GBS.with.MACS3.win.intersect.bed\
#	C13F3_19.GBS.with.MACS3.win.intersect.bed\
#	C13F3_1.GBS.with.MACS3.win.intersect.bed\
#	C13F3_20.GBS.with.MACS3.win.intersect.bed\
#	C13F3_21.GBS.with.MACS3.win.intersect.bed\
#	C13F3_22.GBS.with.MACS3.win.intersect.bed\
#	C13F3_24.GBS.with.MACS3.win.intersect.bed\
#	C13F3_25.GBS.with.MACS3.win.intersect.bed\
#	C13F3_27.GBS.with.MACS3.win.intersect.bed\
#	C13F3_2.GBS.with.MACS3.win.intersect.bed\
#	C13F3_3.GBS.with.MACS3.win.intersect.bed\
#	C13F3_4.GBS.with.MACS3.win.intersect.bed\
#	C13F3_5.GBS.with.MACS3.win.intersect.bed\
#	C13F3_7.GBS.with.MACS3.win.intersect.bed\
#	C13F3_8.GBS.with.MACS3.win.intersect.bed\
#	C13F3_9.GBS.with.MACS3.win.intersect.bed > $inter_out/joint_win_MACS3_all_individuals.bed

# Now do for each of the windows see how repetitions we have #
#cd $inter_out
#file_list=($inter_out/*GBS.with.verfied.win.intersect)
#for file in "${file_list[@]}"; do
#	cut $file -f 1-3 | uniq > $file.uniq.intersect
#	fileA=$file
#	fileB=$file.uniq.intersect

	# Loop through each line in file A
#	while IFS= read -r line
#	do
	    # Use grep to count the occurrences of the line in file B
	    #count=$(grep -c $line $fileB)
#	    echo "$line" >> $file.summary.freq.GBS.intersect.veri.win
#	    grep -c $line $fileA >> $file.summary.freq.GBS.intersect.veri.win
#	done < $fileB
#done

# Do the intersect of F2_6 GBS with its own indvidually called MACS3 wind #
# bedtools intersect -bed -abam /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/aligned/C13F2_6_Mouse.unique.sorted.bam -b /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/peaks_MQ10_MACS3/C13F2_6_Mouse.MQ20_peaks.narrowPeak > /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/F2_6_intersect_GBS_and_MACS3_win.bed

# Now we understand that MACS3 what does with multiple files is that it just merge them,
#but it also has a certain false positive calling when GBS data is nearby GBS-MEDIP data,
#so let's first do the GBS vs MACS3 inv called and after intersect of at least 80% indv#
#that's already done, yaaay, so let's do intersect with all the rest of the files
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/GBS_with_MACS3.win
# in joint_win_MACS3_all_individuals.bed, third column there is the number of ocurrences for each win, it shows the number of individuals that has it #

# See how many in all individuals #
#awk '$4>56' joint_win_MACS3_all_individuals.bed > windows.MACS3_in_all_invd.bed

# see how many in 80% ind (45 individuos)#
#awk '$4>45' joint_win_MACS3_all_individuals.bed > windows.MACS3_in_80.percent_invd.bed
#is only 1 window in 1 satellite #

# see how many in 70% ind (39 individuos)#
#awk '$4>39' joint_win_MACS3_all_individuals.bed > windows.MACS3_in_70.percent_invd.bed
#8 wind

# see how many in 50% ind (28 individuos)#
#awk '$4>28' joint_win_MACS3_all_individuals.bed > windows.MACS3_in_50.percent_invd.bed
#here there are only 15 windows

# get the regions where the fathers have methylation #
#for F2_15
#out_progeny=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.15
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.15/*_Mouse.genotypes_filtered_gatk.vcf.gz)

#for file in "${file_list[@]}"; do
#	name=${file##*/}
#        x=${name%_Mouse.genotypes_filtered_gatk.vcf}
#	bgzip $file
#	tabix $file.gz
#	bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/C13F2_15.GBS.with.MACS3.win.intersect.bed -o $out_progeny/$x.genotypes.gatk.filtered.by.parental.f2_15 $file.gz
#done

#Now let's restrict more the regions as we have now identified several regions where we now there is differential methylation in at least 50% of the groups
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#	x=${name%_Mouse.genotypes_filtered_gatk.vcf.gz}
#        bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/windows.MACS3.veri.GBS.hyper-hypo.bed -o $out_progeny/$x.genotypes.gatk.filtered.by.DMRs.bed $file
#done

#let's also filter the father
#out_parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win
#parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/no_cluster_filters/C13F2_15.genotypes_refiltered_gatk.vcf
#name=${parental##*/}
#bgzip $parental
#tabix $parental.gz
#x=${name%_Mouse.genotypes_filtered_gatk.vcf.gz}
#bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/windows.MACS3.veri.GBS.hyper-hypo.bed -o $out_parental/$x.genotypes.gatk.filtered.by.DMRs.vcf $parental
#bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/C13F2_15.GBS.with.MACS3.win.intersect.bed -o $out_parental/$x.genotypes.gatk.filtered.by.own.win.vcf $parental.gz

#for F2_16
#out_progeny=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.16
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.16/*_Mouse.genotypes_filtered_gatk.vcf.gz)

#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%_Mouse.genotypes_filtered_gatk.vcf}
        #bgzip $file
        #tabix $file.gz
#        bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/C13F2_16.GBS.with.MACS3.win.intersect.bed \
#	-o $out_progeny/$x.genotypes.gatk.filtered.by.parental.f2_16.vcf $file.gz
#done

#Now let's restrict more the regions as we have now identified several regions where we now there is differential methylation in at least 50% of the groups
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%_Mouse.genotypes_filtered_gatk.vcf.gz}
#        bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/windows.MACS3.veri.GBS.hyper-hypo.bed -o $out_progeny/$x.genotypes.gatk.filtered.by.DMRs.bed $file
#done

#let's also filter the fathers
#out_parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win
#parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/variant_output/C13F2_16_Mouse.genotypes_filtered_gatk.vcf.gz
#name=${parental##*/}
#bgzip $parental
#tabix $parental.gz
#x=${name%_Mouse.genotypes_filtered_gatk.vcf.gz}
#bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/windows.MACS3.veri.GBS.hyper-hypo.bed -o $out_parental/$x.genotypes.gatk.filtered.by.DMRs.vcf $parental
#bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/C13F2_16.GBS.with.MACS3.win.intersect.bed -o $out_parental/$x.genotypes.gatk.filtered.by.own.win.vcf $parental

#for F2_14
#out_progeny=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.14
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.14/*_Mouse.genotypes_filtered_gatk.vcf.gz)

#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%_Mouse.genotypes_filtered_gatk.vcf}
#        bgzip $file
#        tabix $file.gz
#        bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/C13F2_14.GBS.with.MACS3.win.intersect.bed \
#        -o $out_progeny/$x.genotypes.gatk.filtered.by.parental.f2_14.vcf $file.gz
#done

#Now let's restrict more the regions as we have now identified several regions where we now there is differential methylation in at least 50% of the groups
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%_Mouse.genotypes_filtered_gatk.vcf.gz}
#        bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/windows.MACS3.veri.GBS.hyper-hypo.bed -o $out_progeny/$x.genotypes.gatk.filtered.by.DMRs.bed $file
#done

#let's also filter the fathers
#out_parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win
#parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/rerun_filters_vcf/C13F2_14.genotypes_refiltered_gatk.vcf
#name=${parental##*/}
#bgzip $parental
#tabix $parental.gz
#x=${name%.genotypes_refiltered_gatk.vcf.gz}
#bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/windows.MACS3.veri.GBS.hyper-hypo.bed -o $out_parental/$x.genotypes.gatk.filtered.by.DMRs.vcf $parental
#bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/C13F2_14.GBS.with.MACS3.win.intersect.bed -o $out_parental/$x.genotypes.gatk.filtered.by.own.win.vcf $parental.gz

#for F2_2
#out_progeny=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.2
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.2/*_Mouse.genotypes_filtered_gatk.vcf.gz)

#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%_Mouse.genotypes_filtered_gatk.vcf}
#        bgzip $file
#        tabix $file.gz
#        bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/C13F2_2.GBS.with.MACS3.win.intersect.bed \
#        -o $out_progeny/$x.genotypes.gatk.filtered.by.parental.f2_2.vcf $file.gz
#done

#Now let's restrict more the regions as we have now identified several regions where we now there is differential methylation in at least 50% of the groups
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%_Mouse.genotypes_filtered_gatk.vcf.gz}
#        bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/windows.MACS3.veri.GBS.hyper-hypo.bed -o $out_progeny/$x.genotypes.gatk.filtered.by.DMRs.bed $file
#done

#let's also filter the fathers
#out_parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win
#parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/variant_output/C13F2_2_Mouse.genotypes_filtered_gatk.vcf.gz
#name=${parental##*/}
#bgzip $parental
#tabix $parental.gz
#x=${name%_Mouse.genotypes_filtered_gatk.vcf.gz}
#bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/windows.MACS3.veri.GBS.hyper-hypo.bed -o $out_parental/$x.genotypes.gatk.filtered.by.DMRs.vcf $parental
#bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/C13F2_2.GBS.with.MACS3.win.intersect.bed -o $out_parental/$x.genotypes.gatk.filtered.by.own.win.vcf $parental

#for F2_1
#out_progeny=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.1
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.1/*_Mouse.genotypes_filtered_gatk.vcf.gz)

#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%_Mouse.genotypes_filtered_gatk.vcf}
#        bgzip $file
#        tabix $file.gz
#        bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/C13F2_1.GBS.with.MACS3.win.intersect.bed \
#        -o $out_progeny/$x.genotypes.gatk.filtered.by.parental.f2_1.vcf $file.gz
#done

#Now let's restrict more the regions as we have now identified several regions where we now there is differential methylation in at least 50% of the groups
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%_Mouse.genotypes_filtered_gatk.vcf.gz}
#        bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/windows.MACS3.veri.GBS.hyper-hypo.bed -o $out_progeny/$x.genotypes.gatk.filtered.by.DMRs.bed $file
#done

#let's also filter the fathers
#out_parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win
##parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/rerun_filters_vcf/C13F2_1.genotypes_refiltered_gatk.vcf
#name=${parental##*/}
#bgzip $parental
#tabix $parental.gz
#x=${name%.genotypes_refiltered_gatk.vcf.gz}
#bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/windows.MACS3.veri.GBS.hyper-hypo.bed -o $out_parental/$x.genotypes.gatk.filtered.by.DMRs.vcf $parental
#bcftools view -R /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/C13F2_1.GBS.with.MACS3.win.intersect.bed -o $out_parental/$x.genotypes.gatk.filtered.by.own.win.vcf $parental.gz


#see if the windows in obese are in the GBS of the controls
#win=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/GBS_with_MACS3.win

#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/)
###for file in  $(cat "$aligned_input/controlF1.bam.txt"); do
##	bedtools intersect -a $file -u -wa -bed -b $win/windows.MACS3_50percent.obese.F1.bed >> $win/win.in.obese.F1.inGBS.F1.control.bed
#done

#for file in  $(cat "$aligned_input/obeseF1.bam.txt"); do
#        bedtools intersect -a $file -u -wa -bed -b $win/windows.MACS3_50percent.control.F1.bed >> $win/win.in.control.F1.inGBS.F1.obese.bed
#done

#for file in  $(cat "$aligned_input/bamF3Clist.txt"); do
#        bedtools intersect -a $file -u -wa -bed -b $win/windows.MACS3_50percent.obese.F2.bed >> $win/win.in.obese.F2.inGBS.F2.control.bed
#done

#for file in  $(cat "$aligned_input/bamF3Olist.txt"); do
#        bedtools intersect -a $file -bed -u -wa -b $win/windows.MACS3_50percent.control.F2.bed >> $win/win.in.control.F2.inGBS.F2.obese.bed
#done

#bedtools intersect -a /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/merged/mouse_meth_countmatrix.txt -wa -b /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/windows.MACS3.veri.GBS.hyper-hypo.bed > /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/counts.verified.wind.MACS3.bed

####################################
# get all the methylated regions from the fathers ##
#file_list=($parents/*genotypes_refiltered_gatk.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.genotypes_refiltered_gatk.vcf}
#	bgzip $file
#	tabix $file.gz
#	bcftools view -R $inter_out/$x.GBS.with.MACS3.win.intersect.bed -o $parents/$x.refilt.own.win.vcf $file.gz
#done

# Cause we already had them gz samples F2_1 and F2_14 did not filtered, so doing it again but without the gzip step
file_list=($parents/*genotypes_refiltered_gatk.vcf.gz)
for file in "${file_list[@]}"; do
        name=${file##*/}
        x=${name%.genotypes_refiltered_gatk.vcf.gz}
        bcftools view -R $inter_out/$x.GBS.with.MACS3.win.intersect.bed -o $parents/$x.refilt.own.win.vcf $file
done
