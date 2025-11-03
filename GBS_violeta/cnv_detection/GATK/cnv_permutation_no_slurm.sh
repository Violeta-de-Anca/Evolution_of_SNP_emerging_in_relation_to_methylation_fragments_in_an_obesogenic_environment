#!/bin/bash -l
module load BEDTools SAMtools

mkdir -p /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK/cnv_coverage_check

# first also merge all the CNVs in F2 and F3, then merge all that are in the same coordinates, so we have the distribution
counter=1
#then do the random coordinates of similar size fragments
while IFS= read -r line && [[ $counter -le 8 ]]; do
	printf '%s\n' "$line"
        printf '%s\t' "$line" > ${counter}.bed
        for i in {1..1000}; do
                bedtools shuffle -noOverlapping -i ${counter}.bed -g size_chrom_mm39.bed >> /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK/cnv_coverage_check/cnv_${counter}.bed
        done
        ((counter++))
done < cnvs_all_gen_merged.bed


# do annotations so we can see if there is any enrichtment of any RE
RE_input=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump/UCSC.RepElements.gtf

for i in /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK/cnv_coverage_check/cnv_*.bed; do
        a=${i##*/}
        b=${a%.bed}
        c=${b##*cnv_}
        bedtools intersect -wa -loj -a ${i} -b $RE_input > /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK/cnv_coverage_check/annotation_cnv_random_${c}.bed
done

#now count the amount of RE types and store them
touch cnv_repeat_counts.tsv
for i in /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK/cnv_coverage_check/annotation_cnv_random_*.bed; do
        a=${i##*/}
        b=${a%.bed}
        c=${b##*annotation_cnv_random_}
        echo $a
        echo $c
        awk 'NR==FNR{
          coord[$1":"$2"-"$3]=1; next
        }
        {
          c=$1":"$2"-"$3
          if (c in coord) {
            k=c FS $5
            cnt[k]++
          }
        }
        END{
          OFS="\t"
          for (k in cnt) {
            split(k, a, FS)
            print a[1], a[2], cnt[k]
          }
        }' /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK/cnv_coverage_check/cnv_${c}.bed /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK/cnv_coverage_check/annotation_cnv_random_${c}.bed | sort -k1,1 -k2,2 >> cnv_repeat_counts.tsv
done
