import sys
# First, let's define a function to read file B (BED file) and store it in a structured way.
def read_bed_file(bed_file_path):
    bed_ranges = {}
    with open(bed_file_path, 'r') as file:
        for line in file:
            chrom, start, end, count = line.strip().split()
            if chrom not in bed_ranges:
                bed_ranges[chrom] = []
            bed_ranges[chrom].append((int(start), int(end), int(count)))
    return bed_ranges

# Now, let's define a function to process file A and perform the calculations.
def process_file_a(file_a_path, bed_ranges):
    results = []
    with open(file_a_path, 'r') as file:
        for line in file:
            chrom, pos, ref, alt, *individuals = line.strip().split()
            pos = int(pos)
            for start, end, count in bed_ranges.get(chrom, []):
                if start <= pos <= end:
                    for i in range(0, len(individuals), 2):
                        name = individuals[i]
                        genotype = individuals[i+1]
                        # Skip entries without a genotype.
                        if genotype in {'.', './.'}:
                            continue
                        results.append((chrom, start, end, count, name, 1))  # Assuming you want to count entries.
    return results

# Function to sum counts for each unique combination of range and individual.
def summarize_results(results):
    summary = {}
    for chrom, start, end, count, name, value in results:
        key = (chrom, start, end, count, name)
        if key not in summary:
            summary[key] = 0
        summary[key] += value
    return summary

# Main function to tie everything together.
def main(file_a_path, file_b_path):
    bed_ranges = read_bed_file(file_b_path)
    raw_results = process_file_a(file_a_path, bed_ranges)
    summary = summarize_results(raw_results)
    
    # Output the final table
    for (chrom, start, end, count, name), value in summary.items():
        print(f"{chrom}\t{start}\t{end}\t{count}\t{name}\t{value}")

# Replace 'path_to_file_a' and 'path_to_file_b' with the actual paths to your files.
main(sys.argv[1], sys.argv[2])
