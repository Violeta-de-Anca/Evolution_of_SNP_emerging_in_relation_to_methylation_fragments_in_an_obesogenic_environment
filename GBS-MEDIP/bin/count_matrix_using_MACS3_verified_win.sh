#!/bin/bash -l
#SBATCH -A naiss2023-22-848
#SBATCH -p core -n 1
##SBATCH --mem=80gb
#SBATCH -t 4-100:00:00
#SBATCH -J counts.per.ind
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/counts.per.ind.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/counts.per.ind.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load plink
module load BEDTools
module load bcftools
module load subread

input_gbsmedip_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/aligned
input_win_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/GBS_with_MACS3.win
output=/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/merged

#let's first compare bedtools, with MEDIPS and featureCounts
#bedtools intersect, which counts the whole number of ocurrences, independent of is paired end or anything
#file_list=($input_gbsmedip_dir/*_Mouse.unique.MQ20.bam)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%_Mouse.unique.MQ20.bam}
#	bedtools intersect -a $input_win_dir/$x.GBS.with.MACS3.win.intersect.bed -c -b $file > $output/$x.count.bedtools.intersect.txt
#done

#featureCounts, which counts features, using paired end mode
#First transform the bed files of MACS3 into SAF format: GeneID	Chr	Start	End	Strand
#file_list=($input_win_dir/*GBS.with.MACS3.win.intersect.bed)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.GBS.with.MACS3.win.intersect.bed}
#	awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' $file > $input_win_dir/$x.saf
#	echo -e "GeneID\tChr\tStart\tEnd\tStrand" | cat - $input_win_dir/$x.saf > $input_win_dir/$x.GBS.with.MACS3.win.intersect.saf
#	rm $input_win_dir/$x.saf
#done

# After having the SAF files, let's do the featureCounts
#featureCounts -F SAF -p -B -o $output/F1_2.featurecounts.txt -d 0 -D 1500 -T 1 -a $input_win_dir/C13F1_2.GBS.with.MACS3.win.intersect.saf $input_gbsmedip_dir/C13F1_2_Mouse.unique.MQ20.bam

#Do it counting the fragments, not the reads!
#file_list=($input_gbsmedip_dir/*_Mouse.unique.MQ20.bam)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%_Mouse.unique.MQ20.bam}
#	featureCounts -F SAF -p -B --countReadPairs -o $output/$x.featurecounts.txt -d 0 -D 1500 -T 1 -a $input_win_dir/$x.GBS.with.MACS3.win.intersect.saf $file
#done

# Now let's merge all the individual peaks
#for file in *.featurecounts.txt; do tail -n +3 "$file" | cut -f 2-4 >> combined_regions.bed; done
#bedtools sort -i combined_regions.bed > combined_regions_sorted.bed
#bedtools merge -i combined_regions_sorted.bed > combined_regions.bed

#now get the right format from featureCounts
#file_list=($output/*.featurecounts.txt)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.featurecounts.txt}
#	cut -f 2-4,7 $file | tail -n +3 > $output/$x_featureCounts.bed
#done

# get the counts for all the individuals
#file_list=($output/*.featureCounts.bed)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.featureCounts.bed}
#        bedtools intersect -a $output/combined_regions.bed -b $file -wo > $output/$x.ind.counts.featurecounts.txt
#done

#Now get only the right columns
#file_list=($output/*.ind.counts.featurecounts.txt)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.ind.counts.featurecounts.txt}
#	cut -f 1-3,7 $file > $output/$x.win.featurecounts.bed
#done

#They need to be sorted
#file_list=($output/*.win.featurecounts.bed)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.win.featurecounts.bed}
#	bedtools sort -i $file > $output/$x.sorted.win.featurecounts.bed
#done

# Now the chr, start and stop need to be merged
#file_list=($output/*.sorted.win.featurecounts.bed)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.sorted.win.featurecounts.bed}
#	awk -F $'\t' '{print $1"_"$2"_"$3"\t"$4}' $file | sort -k1,1 > $output/$x.merged.featurecounts.bed
#done

# You also need to tweak the combined regions
#awk -F $'\t' '{print $1"_"$2"_"$3}' $output/combined_regions.bed | sort -k1,1 > $output/nospacecombined_regions.bed

#Now get the final count matrix!
#file_list=($output/*.merged.featurecounts.bed)
cd $output
#for file in "${file_list[@]}"; do
for file in $(ls C13F*.merged.featurecounts.bed);do
	name=${file##*/}
	x=${name%.merged.featurecounts.bed}
	join -a 1 -e '0' -1 1 -o '1.1 2.2' $output/nospacecombined_regions.bed $file > $output/temp1.bed
	echo -e "bad_column\t$x\n$(cat $output/temp1.bed)" > $output/temp.bed
	paste -d '\t' $output/temp.bed $output/MACS3.count.matrix.txt > $output/merged.counts
	mv $output/merged.counts $output/MACS3.count.matrix.txt
done

#file="myfile.txt"
#while IFS= read -r line; do
#	join -a 1 -e '0' -1 1 -o '1.1 2.2' $output/nospacecombined_regions.bed $file > $output/temp1.bed
#	echo -e "bad_column\t$x\n$(cat $output/temp1.bed)" > $output/temp.bed
#	paste -d '\t' $output/temp.bed $output/MACS3.count.matrix.txt > $output/merged.counts
#	mv $output/merged.counts $output/MACS3.count.matrix.txt
#done < "$file"

awk '{ for (i=2; i<=NF; i+=2) printf("%s\t", $i); printf("\n"); }' $output/MACS3.count.matrix.txt > $output/counts.from.ind.MACS3.txt
header="chr\tstart\tend"
echo $header > $output/combined.sorted.regions.bed
awk -F $'_' '{print $1"\t"$2"\t"$3}' $output/nospacecombined_regions.bed >> $output/combined.sorted.regions.bed
paste -d '\t' $output/combined.sorted.regions.bed $output/counts.from.ind.MACS3.txt > $output/MACS3.final.count.matrix.featureCounts.txt
