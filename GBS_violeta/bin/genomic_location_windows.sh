#!/bin/bash -l
#SBATCH -A naiss2023-22-848
#SBATCH -p core -n 1
##SBATCH --mem=80gb
#SBATCH -t 10:00:00
#SBATCH -J annotation
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/genomic_features_veri_win.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/genomic_features_veri_win.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load BEDTools

working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta
gtf_folder=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump
indv_win=$working_dir/aligned
output=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/annotation
tss_input=$gtf_folder/promoter_TSS_region_mus_musculus_mm39.bed
exon_input=$gtf_folder/CDS_mus_musculus_mm39.sorted.merged.bed
intron_input=$gtf_folder/intronic_mus_musculus_mm39.bed
three_utr_input=$gtf_folder/3_UTR_mus_musculus_mm39.sorted.merged.bed
five_utr_input=$gtf_folder/5_UTR_mus_musculus_mm39.sorted.merged.bed
start_input=$gtf_folder/start_codon_mus_musculus_mm39.sorted.merged.bed
stop_input=$gtf_folder/stop_codon_mus_musculus_mm39.sorted.merged.bed
inter_input=$gtf_folder/intergenic_mus_musculus_mm39.bed

# From where I downloaded the GTF file
#rsync -av rsync://ftp.ensembl.org/ensembl/pub/current_gtf/mus_musculus/* .

# Now let's do the verfied windows #
# We use the verified_win_ensembl_nomenclature.bed as from there it has been removed the chr from the chromosome name
#bedtools intersect -wa -loj -a $verified_input/verified_win_MACS3_ensembl_nomenclature.bed -b $gtf_file > $verified_input/genomic.features.verfied.DMRs.bed

#with the GFT from ensembl I didn't found anything, so I am going to try with https://remap.univ-amu.fr/download_page #
remap_input=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/remap/remap2022_nr_macs2_mm39_v1_0.bed

#bedtools intersect -wa -loj -a $verified_input/verified_windows_wide_peak.bed -b $remap_input > $verified_input/regulatory_sites_remap_mm39.verfied.DMRs.bed

# With the GTF from the repeated elements we are going to exclude all the RE we have in fragments #
RE_input=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump/UCSC.RepElements.gtf

#let's do also the repeated elements
#file_list=($indv_win/*_Mouse.bam.sorted.merged.bed)
#for i in "${file_list[@]}";do
#        name=${i##*/}
#        x=${name%_Mouse.bam.sorted.merged.bed}
#        bedtools intersect -wa -loj -a $i -b $RE_input > $output/${x}.RE.UCSC.GBS.bed
#done

#cuantify for each ind how much of RE #
#file_list=($output/*.RE.UCSC.GBS.bed)
#for i in "${file_list[@]}";do
#       name=${i##*/}
#       x=${name%.RE.UCSC.GBS.bed}
#       cut -f 6,12 $i | sort | uniq > $output/${x}.uniq.RE.UCSC
#done

#file_list=($output/*.uniq.RE.UCSC)
#for i in "${file_list[@]}";do
#        name=${i##*/}
#        x=${name%.uniq.RE.UCSC}
#        cut -f 1 $i | uniq > $output/${x}.types
#        rm $output/$x.ocurrences.type.RE.UCSC
#        while IFS= read -r line; do
#                echo "$line" >> $output/${x}.ocurrences.type.RE.UCSC
#                grep -c "$line" $i >> $output/${x}.ocurrences.type.RE.UCSC
#        done < $output/${x}.types
#       rm $output/$x.types
#done

#put together all the info from all indv
#tail -n +1 $output/*.ocurrences.type.RE.UCSC > $output/summary.RE.UCSC.GBS

###########################################################

non_repeated=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump/non_repeated_regions_mm39_mus_msuculus.bed
#let's do non-repeated regions
#cd $indv_win
#for i in $(ls *.unique.sorted.bam);do
#       name=${i##*/}
#       x=${name%.unique.sorted.bam}
#       bedtools intersect -wa -bed -wb -a $i -b $non_repeated > $indv_win/${x}.GBS.ucsc.non_repeated.bed
#done

#cd $output
#for i in $(ls *.GBS.ucsc.non_repeated.bed);do
#        echo $i >> summary.ncbi.non-repeated.GBS
#        sort $i | uniq | wc -l >> summary.ncbi.non-repeated.GBS
#done

###########################################################

# we are going to annotate with ncbi database got from ucsc, check script from GBS-MEDIP to know URL

#get bam2bed
#cd $indv_win
#for i in $(ls *.unique.sorted.bam);do
#	name=${i##*/}
#	x=${name%.unique.sorted.bam}
#	bedtools bamtobed -i $i | bedtools sort -i - | bedtools merge -i - > $indv_win/$x.bam.sorted.merged.bed
#done

#EXONS
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/aligned/*.sorted.merged.bed)
#echo "${file_list[@]}"
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%_Mouse.bam.sorted.merged.bed}
#	echo "Output file: $output/${x}_exon_ncbi.GBS.bed"
#	bedtools intersect -bed -wb -wa -a $file -b /proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump/CDS_mus_musculus_mm39.sorted.merged.bed > $output/${x}.exon_ncbi.GBS.bed
#done

#cd $output
#for i in $(ls *.exon_ncbi.GBS.bed);do
#	echo $i >> summary.ncbi.exons.GBS
#	sort $i | uniq | wc -l >> summary.ncbi.exons.GBS
#done

########################################################

#INTERGENIC
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/aligned/*sorted.merged.bed)
#for i in "${file_list[@]}"; do
#        name=${i##*/}
#        x=${name%_Mouse.bam.sorted.merged.bed}
#        bedtools intersect -bed -wb -wa -a $i -b $inter_input > $output/${x}_intergenic_ncbi.GBS.bed
#done

#cd $output
#for i in $(ls *_intergenic_ncbi.GBS.bed);do
#        echo $i >> summary.ncbi.intergenic.GBS
#        sort $i | uniq | wc -l >> summary.ncbi.intergenic.GBS
#done

########################################################

#promotors (2000bps) + TSS
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/aligned/*sorted.merged.bed)
#for i in "${file_list[@]}"; do
#        name=${i##*/}
#        x=${name%_Mouse.bam.sorted.merged.bed}
#        bedtools intersect -bed -wb -wa -a $i -b $tss_input > $output/${x}_promotor.TSS_ncbi.GBS.bed
#done

#cd $output
#for i in $(ls *_promotor.TSS_ncbi.GBS.bed);do
#        echo $i >> summary.ncbi.promotor.TSS.GBS
#	sort $i | uniq | wc -l >> summary.ncbi.promotor.TSS.GBS
#done

########################################################

# Intronic
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/aligned/*sorted.merged.bed)
#for i in "${file_list[@]}"; do
#	name=${i##*/}
#	x=${name%_Mouse.bam.sorted.merged.bed}
#	bedtools intersect -bed -wb -wa -a $i -b $intron_input > $output/${x}_introns_ncbi.GBS.bed
#done

#cd $output
#for i in $(ls *_introns_ncbi.GBS.bed);do
#	echo $i >> summary.ncbi.introns.GBS
#	sort $i | uniq | wc -l >> summary.ncbi.introns.GBS
#done

# 5' UTR
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/aligned/*sorted.merged.bed)
#for i in "${file_list[@]}"; do
#	name=${i##*/}
#	x=${name%_Mouse.bam.sorted.merged.bed}
#	bedtools intersect -bed -wb -wa -a $i -b $five_utr_input > $output/${x}_5.UTR_ncbi.GBS.bed
#done

#cd $output
#for i in $(ls *_5.UTR_ncbi.GBS.bed);do
#        echo $i >> summary.ncbi.5_UTR.GBS
#        sort $i | uniq | wc -l >> summary.ncbi.5_UTR.GBS
#done

# 3' UTR
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/aligned/*sorted.merged.bed)
#for i in "${file_list[@]}"; do
#        name=${i##*/}
#        x=${name%_Mouse.bam.sorted.merged.bed}
#        bedtools intersect -bed -wb -wa -a $i -b $three_utr_input > $output/${x}_3.UTR_ncbi.GBS.bed
#done

#cd $output
#for i in $(ls *_3.UTR_ncbi.GBS.bed);do
#        echo $i >> summary.ncbi.3_UTR.GBS
#        sort $i | uniq | wc -l >> summary.ncbi.3_UTR.GBS
#done

#start codon
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/aligned/*sorted.merged.bed)
#for i in "${file_list[@]}"; do
#        name=${i##*/}
#        x=${name%_Mouse.bam.sorted.merged.bed}
#        bedtools intersect -bed -wb -wa -a $i -b $start_input > $output/${x}_start.codon_ncbi.GBS.bed
#done

#cd $output
#for i in $(ls *_start.codon_ncbi.GBS.bed);do
#        echo $i >> summary.ncbi.start_codon.GBS
#        sort $i | uniq | wc -l >> summary.ncbi.start_codon.GBS
#done

#stop codon
#file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/aligned/*sorted.merged.bed)
#for i in "${file_list[@]}"; do
#        name=${i##*/}
#        x=${name%_Mouse.bam.sorted.merged.bed}
#        bedtools intersect -bed -wb -wa -a $i -b $stop_input > $output/${x}_stop.codon_ncbi.GBS.bed
#done

#cd $output
#for i in $(ls *_stop.codon_ncbi.GBS.bed);do
#        echo $i >> summary.ncbi.stop_codon.GBS
#        sort $i | uniq | wc -l >> summary.ncbi.stop_codon.GBS
#done

###################################################################################################################################

#bedtools intersect -wa -loj -a $verified_input/verified_windows_wide_peak.bed -b $RE_input > $verified_input/rep.elem_veri_win.bed

# let's do the genomic feature part but for every individual #
#cd $indv_win
#for i in $(ls *MQ20_peaks.narrowPeak);do
#	name=${i##*/}
#        x=${name%_Mouse.MQ20_peaks.narrowPeak}
#	bedtools intersect -wa -loj -a $i -b $ensembl_with_chr > $verified_input/intersect_win_per_indv/$x.ensembl.genomic.location.win.MACS3.MQ10.bed
#done

#let's see ind wind if they fall in RE #
#cd $indv_win
#for i in $(ls *MQ20_peaks.narrowPeak);do
#        name=${i##*/}
#        x=${name%_Mouse.MQ20_peaks.narrowPeak}
#        bedtools intersect -wa -loj -a $i -b $RE_input > $verified_input/intersect_win_per_indv/$x.RE.UCSC.win.MACS3.MQ10.bed
#done

# Now let's see how many unique geneIDs we have #
#input_inv=$verified_input/intersect_win_per_indv/genomic_ensembl

#cd $input_inv
#for i in $(ls *.ensembl.genomic.location.win.MACS3.MQ10.bed);do
#	name=${i##*/}
#        x=${name%.ensembl.genomic.location.win.MACS3.MQ10.bed}
#	cut -f 13,19 $i > $x.geno.loc
#	cut -d " " -f 2 $x.geno.loc > $x.ids.geno
#	 cut -f 1 $x.geno.loc > $x.typ.geno
#	 paste $x.typ.geno $x.ids.geno| sort | uniq > $x.uniq.genomic.locations.ensembl
#done

#cuantify for each ind how much of everything#
#cd $input_inv
#for i in $(ls *.uniq.genomic.locations.ensembl);do
#	name=${i##*/}
#	x=${name%.uniq.genomic.locations.ensembl}
#	cut -f 1 $i | uniq > $x.types
#	rm $x.ocurrences.type.genomic.ensembl
#	while IFS= read -r line; do
#		echo "$line" >> $x.ocurrences.type.genomic.ensembl
#		grep -c "$line" $i >> $x.ocurrences.type.genomic.ensembl
#	done < $x.types
#done

#Now let's do regulatory regions with ReMap #
#input_inv=$verified_input/intersect_win_per_indv/*
#cd $indv_win
#for i in $(ls *MQ20_peaks.narrowPeak);do
#        name=${i##*/}
#        x=${name%_Mouse.MQ20_peaks.narrowPeak}
#        bedtools intersect -wa -loj -a $i -b $remap_input > $input_inv/$x.regulatory_sites_remap_mm39.bed
#done


#cuantify for each ind how much of RE #
# Now let's see how many unique geneIDs we have #
#input_inv=$verified_input/intersect_win_per_indv/repeated_elements

#cd $input_inv
#for i in $(ls *.RE.UCSC.win.MACS3.MQ10.bed);do
#       name=${i##*/}
#        x=${name%.RE.UCSC.win.MACS3.MQ10.bed}
#       cut -f 13,19 $i > $x.geno.loc
#       cut -d " " -f 2 $x.geno.loc > $x.ids.geno
#        cut -f 1 $x.geno.loc > $x.typ.geno
#        paste $x.typ.geno $x.ids.geno| sort | uniq > $x.uniq.RE.UCSC
#done

#cd $input_inv
#for i in $(ls *.uniq.RE.UCSC);do
#        name=${i##*/}
#        x=${name%.}
#        cut -f 1 $i | uniq > $x.types
#        rm $x.ocurrences.RE.UCSC.bed
#        while IFS= read -r line; do
#                echo "$line" >> $x.ocurrences.type.RE.UCSC
#                grep -c "$line" $i >> $x.ocurrences.type.RE.UCSC
#        done < $x.types
#done

#####################################################################################

# Let's get the names so we can do pathway analysis

#####################################################################################

## Now do this for each group, so let's do 100% overlap per feature per group ##
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#echo "C13F1_1.uniq.genomic.locations.ensembl" >> control.F1.uniq
#echo "C13F1_2.uniq.genomic.locations.ensembl" >> control.F1.uniq
#echo "C13F1_3.uniq.genomic.locations.ensembl" >> control.F1.uniq
#echo "C13F1_4.uniq.genomic.locations.ensembl" >> control.F1.uniq
#echo "C13F1_5.uniq.genomic.locations.ensembl" >> control.F1.uniq
#echo "C13F1_6.uniq.genomic.locations.ensembl" >> control.F1.uniq
#echo "C13F1_7.uniq.genomic.locations.ensembl" >> obese.F1.uniq
#echo "C13F1_8.uniq.genomic.locations.ensembl" >> obese.F1.uniq
#echo "C13F1_9.uniq.genomic.locations.ensembl" >> obese.F1.uniq
#echo "C13F1_10.uniq.genomic.locations.ensembl" >> obese.F1.uniq
#echo "C13F1_11.uniq.genomic.locations.ensembl" >> obese.F1.uniq
#echo "C13F2_13.uniq.genomic.locations.ensembl" >> obese.F2.uniq
#echo "C13F2_14.uniq.genomic.locations.ensembl" >> obese.F2.uniq
#echo "C13F2_15.uniq.genomic.locations.ensembl" >> obese.F2.uniq
#echo "C13F2_16.uniq.genomic.locations.ensembl" >> obese.F2.uniq
#echo "C13F2_17.uniq.genomic.locations.ensembl" >> obese.F2.uniq
#echo "C13F2_18.uniq.genomic.locations.ensembl" >> obese.F2.uniq
#echo "C13F2_19.uniq.genomic.locations.ensembl" >> obese.F2.uniq
#echo "C13F2_20.uniq.genomic.locations.ensembl" >> obese.F2.uniq
#echo "C13F2_21.uniq.genomic.locations.ensembl" >> obese.F2.uniq
#echo "C13F2_22.uniq.genomic.locations.ensembl" >> obese.F2.uniq
#echo "C13F2_23.uniq.genomic.locations.ensembl" >> obese.F2.uniq
#echo "C13F2_24.uniq.genomic.locations.ensembl" >> obese.F2.uniq
#echo "C13F2_12.uniq.genomic.locations.ensembl" >> control.F2.uniq
#echo "C13F2_1.uniq.genomic.locations.ensembl" >> control.F2.uniq
#echo "C13F2_2.uniq.genomic.locations.ensembl" >> control.F2.uniq
#echo "C13F2_3.uniq.genomic.locations.ensembl" >> control.F2.uniq
#echo "C13F2_4.uniq.genomic.locations.ensembl" >> control.F2.uniq
#echo "C13F2_5.uniq.genomic.locations.ensembl" >> control.F2.uniq
#echo "C13F2_6.uniq.genomic.locations.ensembl" >> control.F2.uniq
#echo "C13F2_7.uniq.genomic.locations.ensembl" >> control.F2.uniq
#echo "C13F2_8.uniq.genomic.locations.ensembl" >> control.F2.uniq
#echo "C13F2_9.uniq.genomic.locations.ensembl" >> control.F2.uniq
#echo "C13F3_10.uniq.genomic.locations.ensembl" >> control.F3.uniq
#echo "C13F3_11.uniq.genomic.locations.ensembl" >> control.F3.uniq
#echo "C13F3_12.uniq.genomic.locations.ensembl" >> control.F3.uniq
#echo "C13F3_13.uniq.genomic.locations.ensembl" >> control.F3.uniq
#echo "C13F3_1.uniq.genomic.locations.ensembl" >> control.F3.uniq
#echo "C13F3_2.uniq.genomic.locations.ensembl" >> control.F3.uniq
#echo "C13F3_3.uniq.genomic.locations.ensembl" >> control.F3.uniq
#echo "C13F3_4.uniq.genomic.locations.ensembl" >> control.F3.uniq
#echo "C13F3_5.uniq.genomic.locations.ensembl" >> control.F3.uniq
#echo "C13F3_7.uniq.genomic.locations.ensembl" >> control.F3.uniq
#echo "C13F3_8.uniq.genomic.locations.ensembl" >> control.F3.uniq
#echo "C13F3_9.uniq.genomic.locations.ensembl" >> control.F3.uniq
#echo "C13F3_14.uniq.genomic.locations.ensembl" >> obese.F3.uniq
#echo "C13F3_15.uniq.genomic.locations.ensembl" >> obese.F3.uniq
#echo "C13F3_16.uniq.genomic.locations.ensembl" >> obese.F3.uniq
#echo "C13F3_17.uniq.genomic.locations.ensembl" >> obese.F3.uniq
#echo "C13F3_18.uniq.genomic.locations.ensembl" >> obese.F3.uniq
#echo "C13F3_19.uniq.genomic.locations.ensembl" >> obese.F3.uniq
#echo "C13F3_20.uniq.genomic.locations.ensembl" >> obese.F3.uniq
#echo "C13F3_21.uniq.genomic.locations.ensembl" >> obese.F3.uniq
#echo "C13F3_22.uniq.genomic.locations.ensembl" >> obese.F3.uniq
#echo "C13F3_24.uniq.genomic.locations.ensembl" >> obese.F3.uniq
#echo "C13F3_25.uniq.genomic.locations.ensembl" >> obese.F3.uniq
#echo "C13F3_27.uniq.genomic.locations.ensembl" >> obese.F3.uniq

# 100% overlapp for control f1 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#	cat $line >> overlapp
#done < control.F1.uniq

#sort overlapp | uniq -c | grep '^ *6 ' | sed 's/^ *6 //' > 100.percent.C.F1.geno.ensembl
#rm overlapp

# 83% overlapp for control f1 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < control.F1.uniq

#sort overlapp | uniq -c | grep '^ *5 ' | sed 's/^ *5 //' > 83.percent.C.F1.geno.ensembl
#rm overlapp

# 100% overlapp for obese f1 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < obese.F1.uniq
#sort overlapp | uniq -c | grep '^ *5 ' | sed 's/^ *5 //' > 100.percent.O.F1.geno.ensembl
#rm overlapp

# 80% overlapp for obese f1 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < obese.F1.uniq
#sort overlapp | uniq -c | grep '^ *4 ' | sed 's/^ *4 //' > 80.percent.O.F1.geno.ensembl
#rm overlapp

# 100% overlapp for control f2 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < control.F2.uniq

#sort overlapp | uniq -c | grep '^ *10 ' | sed 's/^ *10 //' > 100.percent.C.F2.geno.ensembl
#rm overlapp

# 80% overlapp for control f2 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < control.F2.uniq

#sort overlapp | uniq -c | grep '^ *8 ' | sed 's/^ *8 //' > 80.percent.C.F2.geno.ensembl
#rm overlapp

# 100% overlapp for obese f2 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < obese.F2.uniq
#sort overlapp | uniq -c | grep '^ *12 ' | sed 's/^ *12 //' > 100.percent.O.F2.geno.ensembl
#rm overlapp

# 83% overlapp for obese f2 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < obese.F2.uniq
#sort overlapp | uniq -c | grep '^ *10 ' | sed 's/^ *10 //' > 83.percent.O.F2.geno.ensembl
#rm overlapp

# 100% overlapp for control f3 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < control.F3.uniq

#sort overlapp | uniq -c | grep '^ *12 ' | sed 's/^ *12 //' > 100.percent.C.F3.geno.ensembl
#rm overlapp

# 83% overlapp for control f3 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < control.F3.uniq

#sort overlapp | uniq -c | grep '^ *10 ' | sed 's/^ *10 //' > 83.percent.C.F3.geno.ensembl
#rm overlapp

# 100% overlapp for obese f3 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < obese.F3.uniq
#sort overlapp | uniq -c | grep '^ *12 ' | sed 's/^ *12 //' > 100.percent.O.F3.geno.ensembl
#rm overlapp


# 83% overlapp for obese f3 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < obese.F3.uniq
#sort overlapp | uniq -c | grep '^ *10 ' | sed 's/^ *10 //' > 83.percent.O.F3.geno.ensembl
#rm overlapp

# Now let's do the same but with repeated elements #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/repeated_elements
#sed -i 's/[a-z]*\.uniq.genomic.locations.ensembl*$/.uniq.RE.UCSC/' control.F1.uniq
#sed -i 's/[a-z]*\.uniq.genomic.locations.ensembl*$/.uniq.RE.UCSC/' control.F2.uniq
#sed -i 's/[a-z]*\.uniq.genomic.locations.ensembl*$/.uniq.RE.UCSC/' control.F3.uniq
#sed -i 's/[a-z]*\.uniq.genomic.locations.ensembl*$/.uniq.RE.UCSC/' obese.F1.uniq
#sed -i 's/[a-z]*\.uniq.genomic.locations.ensembl*$/.uniq.RE.UCSC/' obese.F1.uniq
#sed -i 's/[a-z]*\.uniq.genomic.locations.ensembl*$/.uniq.RE.UCSC/' obese.F2.uniq
#sed -i 's/[a-z]*\.uniq.genomic.locations.ensembl*$/.uniq.RE.UCSC/' obese.F3.uniq

#### RE ####
# 100% overlapp for control f1 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/repeated_elements
#while IFS= read -r line; do
#       cat $line >> overlapp
#done < control.F1.uniq

#sort overlapp | uniq -c | grep '^ *6 ' | sed 's/^ *6 //' > 100.percent.C.F1.RE.UCSC
#rm overlapp

# 83% overlapp for control f1 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/repeated_elements
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < control.F1.uniq

#sort overlapp | uniq -c | grep '^ *5 ' | sed 's/^ *5 //' > 83.percent.C.F1.RE.UCSC
#rm overlapp

# 100% overlapp for obese f1 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/repeated_elements
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < obese.F1.uniq
#sort overlapp | uniq -c | grep '^ *5 ' | sed 's/^ *5 //' > 100.percent.O.F1.RE.UCSC
#rm overlapp

# 80% overlapp for obese f1 #
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < obese.F1.uniq
#sort overlapp | uniq -c | grep '^ *4 ' | sed 's/^ *4 //' > 80.percent.O.F1.RE.UCSC
#rm overlapp

# 100% overlapp for control f2 #
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < control.F2.uniq

#sort overlapp | uniq -c | grep '^ *10 ' | sed 's/^ *10 //' > 100.percent.C.F2.RE.UCSC
#rm overlapp

# 80% overlapp for control f2 #
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < control.F2.uniq

#sort overlapp | uniq -c | grep '^ *8 ' | sed 's/^ *8 //' > 80.percent.C.F2.RE.UCSC
#rm overlapp

# 100% overlapp for obese f2 #
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < obese.F2.uniq
#sort overlapp | uniq -c | grep '^ *12 ' | sed 's/^ *12 //' > 100.percent.O.F2.RE.UCSC
#rm overlapp

# 83% overlapp for obese f2 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < obese.F2.uniq
#sort overlapp | uniq -c | grep '^ *10 ' | sed 's/^ *10 //' > 83.percent.O.F2.RE.UCSC
#rm overlapp

# 100% overlapp for control f3 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < control.F3.uniq

#sort overlapp | uniq -c | grep '^ *12 ' | sed 's/^ *12 //' > 100.percent.C.F3.RE.UCSC
#rm overlapp

# 83% overlapp for control f3 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < control.F3.uniq

#sort overlapp | uniq -c | grep '^ *10 ' | sed 's/^ *10 //' > 83.percent.C.F3.RE.UCSC
#rm overlapp

# 100% overlapp for obese f3 #
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < obese.F3.uniq
#sort overlapp | uniq -c | grep '^ *12 ' | sed 's/^ *12 //' > 100.percent.O.F3.RE.UCSC
#rm overlapp

# 83% overlapp for obese f3 #
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/genomic_ensembl
#while IFS= read -r line; do
#        cat $line >> overlapp
#done < obese.F3.uniq
#sort overlapp | uniq -c | grep '^ *10 ' | sed 's/^ *10 //' > 83.percent.O.F3.RE.UCSC
#rm overlapp

###################################################################################################################
