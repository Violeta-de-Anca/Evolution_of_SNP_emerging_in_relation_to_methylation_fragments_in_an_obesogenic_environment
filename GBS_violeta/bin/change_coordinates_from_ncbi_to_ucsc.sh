#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 4
#SBATCH -t 100:00:00
#SBATCH -J changecoordi
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/coordinates.change.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/coordinates.change.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load GATK
module load picard

working_dir=/proj/naiss2023-23-55/GBS_violeta
input_prefix=$working_dir/snp_indel_databases

#change the SNP database from ncbi to ucsc, we only need to add chr, the rest of the coordinates is the same
#zcat $input_prefix/mus_musculus.ensembl.vcf.gz | head -n 18 > $input_prefix/header.db.ensembl
#zcat $input_prefix/mus_musculus.ensembl.vcf.gz | awk 'BEGIN{OFS="\t"}$1="chr"$1' > $input_prefix/ensembl.db
#cat $input_prefix/ensembl.db | tail -n 4078812 >> $input_prefix/header.db.ensembl
#mv $input_prefix/header.db.ensembl $input_prefix/mus_musculus.uscs.ensembl.vcf.gz

#change the SNP database from ncbi to ucsc, we only need to add chr, the rest of the coordinates is the same
zcat $input_prefix/mgp_REL2021_snps.vcf.gz | head -n 54 > $input_prefix/header.snp.ncbi
zcat $input_prefix/mgp_REL2021_snps.vcf.gz | cut -f 1-8 > $input_prefix/mgp_REL2021_snps.onlysnps
cat $input_prefix/mgp_REL2021_snps.onlysnps | awk 'BEGIN{OFS="\t"}$1="chr"$1' > $input_prefix/mgp_REL2021_snps.ucsc
cat $input_prefix/mgp_REL2021_snps.ucsc.onlysnps | tail -n 83212163 >> $input_prefix/header.snp.ncbi
mv $input_prefix/header.snp.ncbi $input_prefix/mgp_REL2021_ucsc.snps.vcf.gz

#rm $input_prefix/mgp_REL2021_snps.ucsc
#rm $input_prefix/mgp_REL2021_snps.ucsc.onlysnps

#change the indel database from ncbi to ucsc, we only need to add chr, the rest of the coordinates is t$
#zcat $input_prefix/mgp_REL2021_indels.vcf.gz | head -n 113 > $input_prefix/header.indel.ncbi
#cat $input_prefix/mgp_REL2021_indels.vcf.gz | awk 'BEGIN{OFS="\t"}$1="chr"$1' > $input_prefix/mgp_REL2021_indels
#cat $input_prefix/mgp_REL2021_indels | cut -f 1-8 -d $"\t" > $input_prefix/mgp_REL2021_indels.ucsc.onlyindels
#cat $input_prefix/mgp_REL2021_indels.ucsc.onlyindels | tail -n 23172883 >> $input_prefix/header.indel.ncbi
#mv $input_prefix/header.indel.ncbi $input_prefix/mgp_REL2021_ucsc.indels.vcf.gz

#rm $input_prefix/mgp_REL2021_indels
#rm $input_prefix/mgp_REL2021_indels.ucsc.onlyindels
