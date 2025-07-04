#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core
#SBATCH -n 10
#SBATCH --mem=80gb
#SBATCH -t 24:00:00
#SBATCH -J RangeCalling
#SBATCH --error /proj/ancestry_medi_indiv/range.calling.err
#SBATCH --output /proj/ancestry_medi_indiv/range.calling.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

#load the necessary libraries
module load bioinfo-tools
module load BEDOPS
module load BEDTools
module load samtools
module load R/4.2.1

#Get the directory path of the script
working_dir=/proj/ancestry_medi_indiv

#path were the BAM files are located
input_dir=$working_dir/aligned

#path for the output folder
output_dir=$working_dir/merged

# Create the 'merged' folder if it doesn't exist
mkdir -p $output_dir

#FIRST: MERGE THE BAM FILES WITHIN THE INPUT DIRECTORY
samtools merge ${output_dir}/Merged.bam ${input_dir}/*.sorted.bam
samtools sort -m 768M  ${output_dir}/Merged.bam -o ${output_dir}/Merged.sorted.bam
#ls "${output_dir}"*.bam | grep -v .sorted.bam$ | xargs rm
samtools index ${output_dir}/Merged.sorted.bam
samtools depth ${output_dir}/Merged.sorted.bam | awk '{sum+=$3;cnt++}END{print sum/cnt" "sum}' > ${output_dir}/depth.merged.txt

#SECOND: Convert BAM to bed file while keeping the strand information of R1 and R2
bedtools bamtobed -i "${output_dir}/Merged.sorted.bam" > "${output_dir}/aligned_reads.bed"

#THIRD: Creating peaks based only in the paired-end nature of the reads
#Option to print only the first three columns
awk 'NR > 1 {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $4); split($4, a, "/"); if (a[2] == "1") {print $1 "\t" $2 "\t" $3} else if (a[2] == "2") {print $1 "\t" $2 "\t" $3}}' ${output_dir}/aligned_reads.bed | sort -u > ${output_dir}/paired_regions_tab_delimited.bed

#FORTH: MERGING THE PAIRED-END READS BASED IN AN OVERLAP THRESHOULD
#Temporary file
tmp_file="tmp.bed"

# Merge peaks and include merged coordinates
bedtools sort -i ${output_dir}/paired_regions_tab_delimited.bed | bedtools merge -i - -c 1,2,3 -o distinct | awk '{print $1,$2,$3}' > ${output_dir}/$tmp_file

# Sort merged regions by chromosome and start position
sort -k1,1 -k2,2n ${output_dir}/$tmp_file > ${output_dir}/merged_regions_BEDTOOLS.bed

# Clean up temporary file
rm "$tmp_file"

# Calculate statistics on the merged regions
awk 'BEGIN{min=99999999; max=0; count=0; total=0; values="";}{width = $3 - $2; sum += width; sumsq += width^2; count++; if (width < min) min = width; if (width > max) max = width; values = values " " width;}{total += width;}END{avg = sum/count; stddev = sqrt(sumsq/count - (sum/count)^2); print "Count: " count; print "Minimum: " min; print "Maximum: " max; print "Average: " avg; print "Standard Deviation: " stddev; n = split(values, arr); asort(arr); q1 = arr[int(n/4)]; q3 = arr[int(n*3/4)]; print "Q1: " q1; print "Q3: " q3;}' "${output_dir}/merged_regions_BEDTOOLS.bed" > ${output_dir}/statistics.merged.regions.txt

#FIFTH - subset windows above and bellow 300 bps, which is the expected fragment length

input_file="${output_dir}/merged_regions_BEDTOOLS.bed"
output_file_below="${output_dir}/subset_windows_below300.bed"
output_file_above="${output_dir}/subset_windows_above300.bed"
min_width=300

while read -r line; do
    range=($line)
    start=${range[1]}
    end=${range[2]}
    width=$((end - start))

    if ((width <= min_width)); then
        echo -e "${range[0]}\t$start\t$end" >> "$output_file_below"
    else
        echo -e "${range[0]}\t$start\t$end" >> "$output_file_above"
    fi
done < "$input_file"

#input_file_below="subset_windows_below300.bed"

# Calculate average and standard deviation using awk
awk 'BEGIN {sum=0; sumsq=0; count=0} {sum+=$3-$2; sumsq+=($3-$2)*($3-$2); count++} END {avg=sum/count; stdev=sqrt(sumsq/count-(sum/count)^2); print "Average:", avg, "Standard Deviation:", stdev}' "$output_file_below" >> ${output_dir}/statistics.merged.regions.txt

# Save average as variable
desired_average=$(awk 'BEGIN {sum=0; sumsq=0; count=0} {sum+=$3-$2; sumsq+=($3-$2)*($3-$2); count++} END {avg=sum/count; stdev=sqrt(sumsq/count-(sum/count)^2); printf("%.0f", avg)}' "$output_file_below")

###############################
# Input file
input_file="${output_dir}/subset_windows_above300.bed"
# Output file
output_file="${output_dir}/sub_peaks_sub_peaks.bed"

# Split size
split_size=$desired_average

# Temporary file
tmp_file="temp.bed"

# Array to store sub-peak sizes
sub_peak_sizes=()

# Read input file line by line
while IFS=$'\t' read -r chrom start end; do
  length=$((end - start))

  # Calculate the number of sub-peaks to generate
  num_sub_peaks=$((length / split_size))

  # Calculate the size of each sub-peak
  sub_peak_size=$((length / num_sub_peaks))

  # Calculate the remaining bases after splitting into sub-peaks
  remainder=$((length % num_sub_peaks))

  # Adjust the size of the sub-peaks to distribute the remainder
  sub_peak_size=$((sub_peak_size + remainder / num_sub_peaks))

  # Generate the sub-peaks
  for ((i = 0; i < num_sub_peaks; i++)); do
    sub_peak_start=$((start + i * sub_peak_size))
    sub_peak_end=$((sub_peak_start + sub_peak_size))

    echo -e "$chrom\t$sub_peak_start\t$sub_peak_end" >> "$tmp_file"

    # Store the size of each sub-peak
    sub_peak_sizes+=("$sub_peak_size")
  done
done < "${output_dir}/subset_windows_above300.bed"

# Sort the sub-peaks
sort -k1,1 -k2,2n "$tmp_file" > "${output_dir}/sub_peaks_sub_peaks.bed"

# Clean up temporary file
rm "$tmp_file"

# esto no se q es, hay q preguntar a fabio
# Input files
input_file1="${output_dir}/subset_windows_below300.bed"
input_file2="${output_dir}/sub_peaks_sub_peaks.bed"

# Output file
output_file="${output_dir}/joined_sorted.bed"

# Concatenate the input files
cat "$input_file1" "$input_file2" > "$output_file"

# Sort the joined file
sort -k1,1 -k2,2n "$output_file" -o "$output_file"

awk 'BEGIN {sum=0; sumsq=0; count=0} {sum+=$3-$2; sumsq+=($3-$2)*($3-$2); count++} END {avg=sum/count; stdev=sqrt(sumsq/count-(sum/count)^2); print "Average:", avg, "Standard Deviation:", stdev}' "$output_file"

echo "The final file is called ${output_file##*/} and is located at ${output_file%/*}"

#The final file is called joined_sorted.bed

echo Finished
