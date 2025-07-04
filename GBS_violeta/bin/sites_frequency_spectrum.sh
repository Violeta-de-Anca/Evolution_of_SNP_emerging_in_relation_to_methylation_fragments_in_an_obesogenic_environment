#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 4
#SBATCH -t 100:00:00
#SBATCH -J SFS
#SBATCH --error /proj/naiss2023-23-55/GBS_violeta/log_files/angsd.SFS.err
#SBATCH --output /proj/naiss2023-23-55/GBS_violeta/log_files/angsd.SFS.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools samtools PCAngsd ANGSD bcftools plink htslib

working_dir=/proj/naiss2023-23-55/GBS_violeta
input_prefix=$working_dir/aligned
sfs_output=$working_dir/sites_frequency_spectrum
reference_genome=/proj/naiss2023-23-55/GBS_violeta/reference_genomes/mus_musculus/uscs_ref/mm39.fa
mkdir -p $sfs_output

cd $sfs_output

#we are going to perform a folded SFS as we don't have an outgroup, lets do first the obese group

angsd -bam $input_prefix/bamF3Olist.txt -doSaf 1 -out sfs_gbs_ICR_F3_O -anc $reference_genome -GL 2 -nThreads 4 -minMapQ 30 -minQ 20

#I am doing the last generation as here is where we would expect to have the most differences, but also should be distinguished by group!

#Now there is a process of optimization
realSFS sfs_gbs_ICR_F3_O.saf.idx -maxIter 100 -P 4 > sfs_gbs_ICR_F3_O.sfs

