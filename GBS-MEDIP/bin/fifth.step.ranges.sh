#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core
#SBATCH -n 10
#SBATCH --mem=80gb
#SBATCH -t 24:00:00
#SBATCH -J fifthRC
#SBATCH --error /proj/ancestry_medi_indiv/from.fifth.step.range.call.err
#SBATCH --output /proj/ancestry_medi_indiv/from.fifth.step.range.call.out
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

#FIFTH - subset windows above and bellow 300 bps, which is the expected fragment length

input_file=${output_dir}/merged_regions_BEDTOOLS.bed
output_file_below=${output_dir}/subset_windows_below300.bed
output_file_above=${output_dir}/subset_windows_above300.bed
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
done < $input_file

# Calculate average and standard deviation using awk
awk 'BEGIN {sum=0; sumsq=0; count=0} {sum+=$3-$2; sumsq+=($3-$2)*($3-$2); count++} END {avg=sum/count; stdev=sqrt(sumsq/count-(sum/count)^2); print "Average:", avg, "Standard Deviation:", stdev}' "$output_file_below" >> ${output_dir}/statistics.merged.regions.txt

# Save average as variable
desired_average=$(awk 'BEGIN {sum=0; sumsq=0; count=0} {sum+=$3-$2; sumsq+=($3-$2)*($3-$2); count++} END {avg=sum/count; stdev=sqrt(sumsq/count-(sum/count)^2); printf("%.0f", avg)}' "$output_file_below")

###############################
# Input file
input_file=${output_dir}/subset_windows_above300.bed
# Output file
output_file=${output_dir}/sub_peaks_sub_peaks.bed

# Split size
split_size=$desired_average

# Temporary file
tmp_file=temp.bed

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
done < ${output_dir}/subset_windows_above300.bed

# Sort the sub-peaks
sort -k1,1 -k2,2n "$tmp_file" > "${output_dir}/sub_peaks_sub_peaks.bed"

# Clean up temporary file
rm $tmp_file

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

echo Finished
