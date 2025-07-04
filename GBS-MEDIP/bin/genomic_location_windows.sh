#!/bin/bash -l
#SBATCH -A naiss2023-22-848
#SBATCH -p core -n 1
#SBATCH -t 1:00:00
#SBATCH -J annoRE
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/genomic_features_veri_win.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/genomic_features_veri_win.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load BEDTools

working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP
verified_input=$working_dir/pathway_enrichment
indv_win=$working_dir/merged
gtf_folder=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump
# From where I downloaded the GTF file
#rsync -av rsync://ftp.ensembl.org/ensembl/pub/current_gtf/mus_musculus/* .
ensembl_with_chr=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/ensembl/Mus_musculus.mm39.gtf

############################################################################################################
#I am going to stick to UCSC ref genome and annotation files, from rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/ . I downloaded all the annotation files needed.
ucsc_input=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump/mm39.ncbiRefSeq.gtf.gz

# Let's extract the TSS and the promotor regions from the UCSC gtf file
# Extract TSS and Promoter regions, output in BED format
#zcat $ucsc_input | grep '5UTR' > $gtf_folder/5_UTR_Mus_musculus.mm39.gtf
#awk 'BEGIN{OFS="\t"} {
#    if ($7=="+") {
#        print $1,$4-2001, $4;
#	} else {
#        print $1,$5+1, $5+2001;
#    }
#}' $gtf_folder/5_UTR_Mus_musculus.mm39.gtf | awk 'BEGIN{FS=OFS="\t"} {if ($2 <= 0) $2 = 1; print}' | bedtools sort -i -| bedtools merge -i - > $gtf_folder/promoter_TSS_region_mus_musculus_mm39.bed

# Let's turn the 5' UTR into bed files
#cat $gtf_folder/5_UTR_Mus_musculus.mm39.gtf | awk 'BEGIN{OFS="\t"} {print $1, $4, $5}' | bedtools sort -i - | bedtools merge -i - > $gtf_folder/5_UTR_mus_musculus_mm39.sorted.merged.bed

# Now the 3' UTR
#zcat $ucsc_input| grep '3UTR' > $gtf_folder/3_UTR.mm39.gtf
#cat $gtf_folder/3_UTR.mm39.gtf | awk 'BEGIN{OFS="\t"} {print $1, $4, $5}' | bedtools sort -i - | bedtools merge -i - > $gtf_folder/3_UTR_mus_musculus_mm39.sorted.merged.bed

# Now the exons
#zcat $ucsc_input | grep 'CDS' > $gtf_folder/CDS_mus_musculus_mm39.gtf
#awk 'BEGIN{OFS="\t"} {print $1, $4, $5}' $gtf_folder/CDS_mus_musculus_mm39.gtf | bedtools sort -i - | bedtools merge -i - > $gtf_folder/CDS_mus_musculus_mm39.sorted.merged.bed

#let's get the intronic regions
#cat $gtf_folder/CDS_mus_musculus_mm39.sorted.merged.bed $gtf_folder/5_UTR_mus_musculus_mm39.sorted.merged.bed $gtf_folder/3_UTR_mus_musculus_mm39.sorted.merged.bed $gtf_folder/promoter_TSS_region_mus_musculus_mm39.bed > $gtf_folder/mergedCDS5UTR3UTR_mm39.bed
#awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' $ensembl_with_chr | bedtools sort | bedtools subtract -a stdin -b $gtf_folder/mergedCDS5UTR3UTR_mm39.bed > $gtf_folder/intronic_mus_musculus_mm39.bed

#let's get the start codons regions
#zcat $ucsc_input | grep 'start_codon' | awk 'BEGIN{OFS="\t"} {print $1, $4, $5}' | bedtools sort -i - | bedtools merge -i - > $gtf_folder/start_codon_mus_musculus_mm39.sorted.merged.bed

#do mergin of 3' UTR and stop codon
#cat $gtf_folder/3_UTR_mus_musculus_mm39.sorted.merged.bed $gtf_folder/stop_codon_mus_musculus_mm39.sorted.merged.bed | bedtools sort -i - | bedtools merge -i - -d 5 > 3_UTR_stop_codon_mus_musculus_mm39_ncbi.bed

#do mergin of 5' UTR and start codon
#cat $gtf_folder/5_UTR_mus_musculus_mm39.sorted.merged.bed $gtf_folder/start_codon_mus_musculus_mm39.sorted.merged.bed | bedtools sort -i - | bedtools merge -i - -d 5 > 5_UTR_start_codon_mus_musculus_mm39_ncbi.bed

#let's get the stop codon regions
#zcat $ucsc_input | grep 'stop_codon' | awk 'BEGIN{OFS="\t"} {print $1, $4, $5}' | bedtools sort -i - | bedtools merge -i - > $gtf_folder/stop_codon_mus_musculus_mm39.sorted.merged.bed

#let's do the intergenic regions
#cat $gtf_folder/mergedCDS5UTR3UTR_mm39.bed $gtf_folder/intronic_mus_musculus_mm39.bed | bedtools sort -i - | bedtools merge -d 600 -i - > $gtf_folder/all_genomic_features_mm39.bed
#bedtools complement -i $gtf_folder/all_genomic_features_mm39.bed -g $gtf_folder/chrom_size_mm39_sorted.bed > $gtf_folder/intergenic_mus_musculus_mm39.bed

# do non-repeated regions of the genome
#bedtools complement -i /proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump/UCSC.RepElements.bed -g $gtf_folder/chrom_size_mm39_sorted.bed > $gtf_folder/non_repeated_regions_mm39_mus_msuculus.

tss_input=$gtf_folder/promoter_TSS_region_mus_musculus_mm39.bed
exon_input=$gtf_folder/CDS_mus_musculus_mm39.sorted.merged.bed
intron_input=$gtf_folder/intronic_mus_musculus_mm39.bed
three_utr_input=$gtf_folder/3_UTR_mus_musculus_mm39.sorted.merged.bed
five_utr_input=$gtf_folder/5_UTR_mus_musculus_mm39.sorted.merged.bed
start_input=$gtf_folder/start_codon_mus_musculus_mm39.sorted.merged.bed
stop_input=$gtf_folder/stop_codon_mus_musculus_mm39.sorted.merged.bed
inter_input=$gtf_folder/intergenic_mus_musculus_mm39.bed
5utrandstart=$gtf_folder/5_UTR_start_codon_mus_musculus_mm39_ncbi.bed
3utrandstop=$gtf_folder/3_UTR_stop_codon_mus_musculus_mm39_ncbi.bed

# We use the verified_win_ensembl_nomenclature.bed as from there it has been removed the chr from the chromosome name
# I am going to try also with https://remap.univ-amu.fr/download_page #
remap_input=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/remap/remap2022_nr_macs2_mm39_v1_0.bed
# And with the GTF from the repeated elements we are going to exclude all the RE we have in fragments #
RE_input=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump/UCSC.RepElements.gtf

#################################################################################################
#let's do the start codon regions
#cd $indv_win
#for i in $(ls *.featureCounts.bed);do
#        name=${i##*/}
#        x=${name%.featureCounts.bed}
#        bedtools intersect -wa -bed -wb -a $i -b $start_input > $verified_input/start_codon/$x.ucsc.start_codon.bed
#done

#cd $verified_input/start_codon
#for i in $(ls *.ucsc.start_codon.bed);do
#        echo $i >> summary.ncbi.start_codon.GBS-MEDIP
#        sort $i | uniq | wc -l >> summary.ncbi.start_codon.GBS-MEDIP
#done

#################################################################################################
#let's do the stop codon regions
#cd $indv_win
#for i in $(ls *.featureCounts.bed);do
#        name=${i##*/}
#        x=${name%.featureCounts.bed}
#        bedtools intersect -wa -bed -wb -a $i -b $stop_input > $verified_input/stop_codon/$x.ucsc.stop_codon.bed
#done

#cd $verified_input/stop_codon
#for i in $(ls *.ucsc.stop_codon.bed);do
#	echo $i >> summary.ncbi.stop_codon.GBS-MEDIP
#	sort $i | uniq | wc -l >> summary.ncbi.stop_codon.GBS-MEDIP
#done

#################################################################################################
#let's do the 3 UTR regions
#cd $indv_win
#for i in $(ls *.featureCounts.bed);do
#	name=${i##*/}
#	x=${name%.featureCounts.bed}
#	bedtools intersect -wa -bed -wb -a $i -b $three_utr_input > $verified_input/3_UTR/$x.ucsc.3_UTR.bed
#done

#cd $verified_input/3_UTR
#for i in $(ls *.ucsc.3_UTR.bed);do
#	echo $i >> summary.ncbi.3_UTR.GBS-MEDIP
#	sort $i | uniq | wc -l >> summary.ncbi.3_UTR.GBS-MEDIP
#done

##################################################################################################

#let's do the 5 UTR regions
#cd $indv_win
#for i in $(ls *.featureCounts.bed);do
#	name=${i##*/}
#	x=${name%.featureCounts.bed}
#	bedtools intersect -wa -bed -wb -a $i -b $five_utr_input > $verified_input/5_UTR/$x.ucsc.5_UTR.bed
#done

#cd $verified_input/5_UTR
#for i in $(ls *.ucsc.5_UTR.bed);do
#        echo $i >> summary.ncbi.5_UTR.GBS-MEDIP
#        sort $i | uniq | wc -l >> summary.ncbi.5_UTR.GBS-MEDIP
#done

#######################################################################################

# let's do the TSS + promotor feature part but for every individual #
#cd $indv_win
#for i in $(ls *.featureCounts.bed);do
#	name=${i##*/}
#	x=${name%.featureCounts.bed}
#	bedtools intersect -wa -bed -wb -a $i -b $tss_input > $verified_input/promoter_TSS/$x.ucsc.tss_promotor.MQ10.bed
#done

#cd $verified_input/promoter_TSS
#for i in $(ls *.ucsc.tss_promotor.MQ10.bed);do
#	echo $i >> summary.ncbi.tss_promotor2000bp.GBS-MEDIP
#	sort $i | uniq | wc -l >> summary.ncbi.tss_promotor2000bp.GBS-MEDIP
#done

#################################################################################

#let´s do the exon, start and stop codons and UTR (untranslated regions)
#exon
#cd $indv_win
#for i in $(ls *.featureCounts.bed);do
#       name=${i##*/}
#       x=${name%.featureCounts.bed}
#       bedtools intersect -wa -bed -wb -a $i -b $exon_input > $verified_input/exons/$x.ucsc.exon.bed
#done

#cd $verified_input/exons
#for i in $(ls *.ucsc.exon.bed);do
#        echo $i >> summary.ncbi.exons.GBS-MEDIP
#        sort $i | uniq | wc -l >> summary.ncbi.exons.GBS-MEDIP
#done

##################################################################################

#let's do introns
#cd $indv_win
#for i in $(ls *.featureCounts.bed);do
#       name=${i##*/}
#       x=${name%.featureCounts.bed}
#       bedtools intersect -wa -bed -wb -a $i -b $intron_input > $verified_input/introns/$x.ucsc.intron.bed
#done

#cd $verified_input/introns
#for i in $(ls *.ucsc.intron.bed);do
#	echo $i >> summary.ncbi.introns.GBS-MEDIP
#	sort $i | uniq | wc -l >> summary.ncbi.introns.GBS-MEDIP
#done

##################################################################################

#let's do intergenic regions
#cd $indv_win
#for i in $(ls *.featureCounts.bed);do
#       name=${i##*/}
#       x=${name%.featureCounts.bed}
#       bedtools intersect -wa -bed -wb -a $i -b $inter_input > $verified_input/intergenic/$x.ucsc.intergenic.bed
#done

#cd $verified_input/intergenic
#for i in $(ls *.ucsc.intergenic.bed);do
#        echo $i >> summary.ncbi.intergenic.GBS-MEDIP
#        sort $i | uniq | wc -l >> summary.ncbi.intergenic.GBS-MEDIP
#done

###################################################################################
non_repeated=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump/non_repeated_regions_mm39_mus_msuculus.bed
#let's do non-repeated regions
#cd $indv_win
#for i in $(ls *.featureCounts.bed);do
#       name=${i##*/}
#       x=${name%.featureCounts.bed}
#       bedtools intersect -wa -bed -wb -a $i -b $non_repeated > $verified_input/intergenic/$x.ucsc.non_repeated.bed
#done

#cd $verified_input/intergenic
#for i in $(ls *.ucsc.non_repeated.bed);do
#        echo $i >> summary.ncbi.non-repeated.GBS-MEDIP
#        sort $i | uniq | wc -l >> summary.ncbi.non-repeated.GBS-MEDIP
#done

###################################################################################

# let's see ind wind if they fall in RE #
#cd $indv_win
#for i in $(ls *.featureCounts.bed);do
#        name=${i##*/}
#        x=${name%.featureCounts.bed}
#        bedtools intersect -wa -loj -a $i -b $RE_input > $verified_input/$x.RE.UCSC.win.MACS3.MQ10.bed
#done

###################################################################################

# Now let's see how many promotors we have#
#cd $verified_input
#for i in $(ls *.ensembl.tss_promotor.MQ10.bed);do
#	name=${i##*/}
#	x=${name%.ensembl.tss_promotor.MQ10.bed}
#	cut -f 13,6 $i > $x.geno.loc
#	cut -d " " -f 2 $x.geno.loc > $x.ids.geno
#	cut -f 1 $x.geno.loc > $x.typ.geno
#	paste $x.typ.geno $x.ids.geno| sort | uniq > $x.uniq.tsspromoter.locations.ensembl
#done

#cuantify for each ind how many promotors#
#cd $verified_input
#for i in $(ls *.uniq.promoter.locations.ensembl);do
#	name=${i##*/}
#	x=${name%.uniq.promoter.locations.ensembl}
#	cut -f 1 $i | uniq > $x.types
#	rm $x.ocurrences.promoter.ensembl
#	while IFS= read -r line; do
#		echo "$line" >> $x.ocurrences.promoter.ensembl
#		grep -c "$line" $i >> $x.ocurrences.promoter.ensembl
#	done < $x.types
#done

#put together all the info from all indv
#cd $verified_input
#tail -n +1 *.ocurrences.promoter.ensembl > summary.promoter.ensembl

##############################################################################

# Now let's see how many TSS we have#
#cd $verified_input
#for i in $(ls *.ensembl.tss.MQ10.bed);do
#        name=${i##*/}
#        x=${name%.ensembl.tss.MQ10.bed}
#        cut -f 13,6 $i > $x.geno.loc
#        cut -d " " -f 2 $x.geno.loc > $x.ids.geno
#        cut -f 1 $x.geno.loc > $x.typ.geno
#        paste $x.typ.geno $x.ids.geno| sort | uniq > $x.uniq.tss.locations.ensembl
#done

#cuantify for each ind how many TSS#
#cd $verified_input
#for i in $(ls *.uniq.tss.locations.ensembl);do
#        name=${i##*/}
#        x=${name%.uniq.tss.locations.ensembl}
#        cut -f 1 $i | uniq > $x.types
#        rm $x.ocurrences.tss.ensembl
#        while IFS= read -r line; do
#                echo "$line" >> $x.ocurrences.tss.ensembl
#                grep -c "$line" $i >> $x.ocurrences.tss.ensembl
#        done < $x.types
#done

#put together all the info from all indv
#cd $verified_input
#tail -n +1 *.ocurrences.tss.ensembl > summary.tss.ensembl

##################################################################################

#Now let's do regulatory regions with ReMap #
#cd $indv_win
#for i in $(ls *.featureCounts.bed);do
#        name=${i##*/}
#        x=${name%.featureCounts.bed}
#        bedtools intersect -wa -loj -a $i -b $remap_input > $verified_input/$x.RS.remap_mm39.bed
#done

###################################################################################

#cuantify for each ind how much of RE #
#cd $verified_input
#for i in $(ls *.RE.UCSC.win.MACS3.MQ10.bed);do
#	name=${i##*/}
#	x=${name%.RE.UCSC.win.MACS3.MQ10.bed}
#	cut -f 6,13 $i | sort | uniq > $x.uniq.RE.UCSC
#done

#cd $verified_input
#for i in $(ls *.uniq.RE.UCSC);do
#        name=${i##*/}
#        x=${name%.uniq.RE.UCSC}
#        cut -f 1 $i | uniq > $x.types
#        rm $x.ocurrences.type.RE.UCSC
#        while IFS= read -r line; do
#                echo "$line" >> $x.ocurrences.type.RE.UCSC
#                grep -c "$line" $i >> $x.ocurrences.type.RE.UCSC
#        done < $x.types
#	rm $x.types
#done

#put together all the info from all indv
#cd $verified_input
#tail -n +1 *.ocurrences.type.RE.UCSC > summary.RE.UCSC

#####################################################################################

#let´s extract exons from the ensembl run
#cd $verified_input
#for i in $(ls *.ensembl.exon.StartStopcodon.UTR.MQ10.bed);do
#	name=${i##*/}
#	x=${name%.ensembl.exon.StartStopcodon.UTR.MQ10.bed}
#	awk '$7=="exon"'  $i > $x.exon.ensembl
#	awk '$7 ~ /five_prime_utr/ || $7 ~ /three_prime_utr/' $i > $x.UTR.ensembl
#	awk '$7 ~ /stop_codon/ || $7 ~ /start_codon/' $i > $x.codonStartStop.ensembl
#done

######################################################################################

#let's get the numbers for exons
#cd $verified_input
#for i in $(ls *.exon.ensembl);do
#       name=${i##*/}
#       x=${name%.exon.ensembl}
#       cut -f 13 $i | cut -f 5,6,18 -d ";" | sort -t ';' -k2,2 | uniq > $x.uniq.exons.ensembl
#done

#cd $verified_input
#for i in $(ls *.uniq.exons.ensembl);do
#        name=${i##*/}
#        x=${name%.uniq.exons.ensembl}
#        wc -l $i > $x.ocurrences.exons.ensembl
#done

#put together all the info from all indv
#cd $verified_input
#tail -n +1 *.ocurrences.exons.ensembl > summary.exons.ensembl

######################################################################################

#let's get the numbers for start/stop codons
#cd $verified_input
#for i in $(ls *.codonStartStop.ensembl);do
#	name=${i##*/}
#	x=${name%.codonStartStop.ensembl}
#	cut -f 7,13 $i | grep 'start_codon' > $x.start.codon
#	cut -f 7,13 $i | grep 'stop_codon' > $x.stop.codon
#done
#for i in $(ls *.start.codon);do
#	name=${i##*/}
#	x=${name%.start.codon}
#	cut -f 6 -d ";" $i | sort | uniq > $x.uniq.start.codon
#done

#for i in $(ls *.stop.codon);do
#        name=${i##*/}
#        x=${name%.stop.codon}
#        cut -f 6 -d ";" $i | sort | uniq > $x.uniq.stop.codon
#done


#for i in $(ls *.uniq.stop.codon);do
#	name=${i##*/}
#	x=${name%.uniq.stop.codon}
#	wc -l $i > $x.occurrences.stop.codon
#done

#for i in $(ls *.uniq.start.codon);do
#	name=${i##*/}
#	x=${name%.uniq.start.codon}
#	wc -l $i > $x.occurrences.start.codon
#done

#put together all the info from all indv
#tail -n +1 *.ocurrences.start.codon > summary.start.codons.ensembl
#tail -n +1 *.ocurrences.stop.codon > summary.stop.codons.ensembl

########################################################################################

#let's get the numbers for UTR regions
#cd $verified_input
#for i in $(ls *.UTR.ensembl);do
#       name=${i##*/}
#       x=${name%.UTR.ensembl}
#        cut -f 7,13 $i | grep 'five_prime_utr' > $x.5.UTR
#	cut -f 7,13 $i | grep 'three_prime_utr' > $x.3.UTR
#done

#for i in $(ls *.5.UTR);do
#	name=${i##*/}
#	x=${name%.5.UTR}
#	cut -f 5,15 -d ";" $i | sort -t ";" -k 2,2 | uniq > $x.uniq.5.UTR.ensembl
#done

#for i in $(ls *.3.UTR);do
#        name=${i##*/}
#        x=${name%.3.UTR}
#        cut -f 5,15 -d ";" $i | sort -t ";" -k 2,2 | uniq > $x.uniq.3.UTR.ensembl
#done

#for i in $(ls *.uniq.5.UTR.ensembl);do
#	name=${i##*/}
#        x=${name%.uniq.5.UTR.ensembl}
#	wc -l $i > $x.ocurrences.5.UTR.ensembl
#done

#for i in $(ls *.uniq.3.UTR.ensembl);do
#	name=${i##*/}
#	x=${name%.uniq.3.UTR.ensembl}
#	wc -l $i > $x.ocurrences.3.UTR.ensembl
#done

#put together all the info from all indv
#tail -n +1 *.ocurrences.5.UTR.ensembl > summary.5.UTR.ensembl

#put together all the info from all indv
#tail -n +1 *.ocurrences.3.UTR.ensembl > summary.3.UTR.ensembl

################################################################################
#get the names of the genes so we can do pathway enrichmen
#for prefix in C13F1_10	C13F1_11	C13F1_1	C13F1_2	C13F1_3	C13F1_4	C13F1_5	C13F1_6	C13F1_7	C13F1_8	C13F1_9	C13F2_12	C13F2_13	C13F2_14	C13F2_15	C13F2_16	C13F2_17	C13F2_18	C13F2_19	C13F2_1	C13F2_20	C13F2_21	C13F2_23	C13F2_24	C13F2_2	C13F2_3	C13F2_4	C13F2_5	C13F2_6	C13F2_7	C13F2_8	C13F2_9	C13F3_10	C13F3_11	C13F3_12	C13F3_13	C13F3_14	C13F3_15	C13F3_16	C13F3_17	C13F3_18	C13F3_19	C13F3_1	C13F3_20	C13F3_21	C13F3_22	C13F3_24	C13F3_25	C13F3_27	C13F3_2	C13F3_3	C13F3_4	C13F3_5	C13F3_7	C13F3_8	C13F3_9; do
#	cat $prefix*.ucsc.3_UTR.bed $prefix*.ucsc.5_UTR.bed $prefix*.ucsc.exon.bed $prefix*.ucsc.intron.bed $verified_input/$prefix*.ucsc.tss_promotor.MQ10.bed | cut -f 5,6,7 | bedtools sort -i - | bedtools merge -i - -d 6 > $verified_input/$prefix.all.meth.win.locations
#done

################################################################################
#get the gene names
#for prefix in C13F1_10  C13F1_11        C13F1_1 C13F1_2 C13F1_3 C13F1_4 C13F1_5 C13F1_6 C13F1_7 C13F1_8 C13F1_9 C13F2_12        C13F2_13        C13F2_14        C13F2_15        C13F2_16        C13F2_17        C13F2_18        C13F2_19        C13F2_1 C13F2_20        C13F2_21        C13F2_23        C13F2_24        C13F2_2 C13F2_3 C13F2_4 C13F2_5 C13F2_6 C13F2_7 C13F2_8 C13F2_9 C13F3_10        C13F3_11        C13F3_12        C13F3_13        C13F3_14        C13F3_15        C13F3_16        C13F3_17        C13F3_18        C13F3_19        C13F3_1 C13F3_20        C13F3_21        C13F3_22        C13F3_24        C13F3_25        C13F3_27        C13F3_2 C13F3_3 C13F3_4 C13F3_5 C13F3_7 C13F3_8 C13F3_9; do
#	zcat $ucsc_input | bedtools intersect -a $verified_input/$prefix.all.meth.win.locations -b stdin -wb | cut -f 12 > $verified_input/$prefix.gene.names.bed
#done

################################################################################

# get only the names, do sort and uniq and then put the names in *.gene.names.bed
for prefix in C13F1_10 C13F1_11        C13F1_1 C13F1_2 C13F1_3 C13F1_4 C13F1_5 C13F1_6 C13F1_7 C13F1_8 C13F1_9 C13F2_12        C13F2_13        C13F2_14        C13F2_15        C13F2_16        C13F2_17        C13F2_18        C13F2_19        C13F2_1 C13F2_20        C13F2_21        C13F2_23        C13F2_24        C13F2_2 C13F2_3 C13F2_4 C13F2_5 C13F2_6 C13F2_7 C13F2_8 C13F2_9 C13F3_10        C13F3_11        C13F3_12        C13F3_13        C13F3_14        C13F3_15        C13F3_16        C13F3_17        C13F3_18        C13F3_19        C13F3_1 C13F3_20        C13F3_21        C13F3_22        C13F3_24        C13F3_25        C13F3_27        C13F3_2 C13F3_3 C13F3_4 C13F3_5 C13F3_7 C13F3_8 C13F3_9; do
	awk -F\" '{print $10}' $verified_input/$prefix.gene.names.bed | sort | uniq > $verified_input/$prefix.final.name.txt
	#mv $verified_input/$prefix.final.name.txt $verified_input/$prefix.gene.names.bed
done

################################################################################

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
