#!/bin/bash -l
module load BEDTools

# With the GTF from the repeated elements we are going to exclude all the RE we have in fragments #
RE_input=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump/UCSC.RepElements.gtf

#first putting the cnvs correctly so bedtools can read it
#awk '{ $1 = "chr" $1; print }' F2_cnvs_final.txt > F2_cnvs_final.bed
#awk '{ $1 = "chr" $1; print }' F3_cnvs_final.txt > F3_cnvs_final.bed
cut -d ' ' -f 1 F2_cnvs_final.bed > chr_F2
cut -d ' ' -f 2 F2_cnvs_final.bed > start_F2
cut -d ' ' -f 3 F2_cnvs_final.bed > end_F2
cut -d ' ' -f 4 F2_cnvs_final.bed > type_F2
cut -d ' ' -f 5 F2_cnvs_final.bed > ind_F2
paste chr_F2 start_F2 end_F2 type_F2 ind_F2 > F2_cnvs_final.bed
cut -d ' ' -f 1 F3_cnvs_final.bed > chr_F3
cut -d ' ' -f 2 F3_cnvs_final.bed > start_F3
cut -d ' ' -f 3 F3_cnvs_final.bed > end_F3
cut -d ' ' -f 4 F3_cnvs_final.bed > type_F3
cut -d ' ' -f 5 F3_cnvs_final.bed > ind_F3
paste chr_F3 start_F3 end_F3 type_F3 ind_F3 > F3_cnvs_final.bed

#run the annotation of RE in the CNV regions
bedtools intersect -wa  -loj -a F2_cnvs_final.bed -b $RE_input > F2_cnvs_GATK_annotated_UCSC_RE.bed
bedtools intersect -wa  -loj -a F3_cnvs_final.bed -b $RE_input > F3_cnvs_GATK_annotated_UCSC_RE.bed
