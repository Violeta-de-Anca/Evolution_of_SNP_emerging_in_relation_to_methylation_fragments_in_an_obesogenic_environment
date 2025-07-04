#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 4
#SBATCH --mem=80gb
#SBATCH -t 100:00:00
#SBATCH -J SFSc
##SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/angsd.SFS.err
##SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/angsd.SFS.out
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/angsd.SFSc.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/angsd.SFSc.out
##SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/angsd.gatk.err
##SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/angsd.gatk.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools samtools PCAngsd ANGSD bcftools plink htslib
#angds v0.940-dirty

working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta
input_prefix=$working_dir/aligned
db_knownsites=$working_dir/snp_indel_databases
allele_freq_output=$working_dir/variant_output
reference_genome=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/mm39.fa
mkdir -p $allele_freq_output

cd $allele_freq_output

# Do the genotype likelihoods, it is set to perform like GATK #
# angsd -b $input_prefix/list.txt -GL 2 -out genotypelikelihood.GATK.ICR.sperm.GBS -doMajorMinor 1 -doMaf 2 -minMapQ 30 -minQ 20 -nThreads 4

# Do the SFS per treatment #
#obese F1
#angsd -bam $input_prefix/obeseF1.bam.txt -doSaf 1 -out sfs_gbs_ICR_obese_F1 -anc $reference_genome -GL 2 -nThreads 4 -minMapQ 30 -minQ 20
#realSFS sfs_gbs_ICR_obese_F1.saf.idx -maxIter 100 -P 4 > ICR.GBS.obese.F1.sfs

#control F1
#angsd -bam $input_prefix/controlF1.bam.txt -doSaf 1 -out sfs_gbs_ICR_control_F1 -anc $reference_genome -GL 2 -nThreads 4 -minMapQ 30 -minQ 20
#realSFS sfs_gbs_ICR_control_F1.saf.idx -maxIter 100 -P 4 > ICR.GBS.control.F1.sfs

# obese F2
angsd -bam $input_prefix/obeseF2.bam.txt -doSaf 1 -out sfs_gbs_ICR_obese_F2 -anc $reference_genome -GL 2 -nThreads 4 -minMapQ 30 -minQ 20
realSFS sfs_gbs_ICR_obese_F2.saf.idx -maxIter 100 -P 4 > ICR.GBS.obese.F2.sfs

#control F1
angsd -bam $input_prefix/controlF2.bam.txt -doSaf 1 -out sfs_gbs_ICR_control_F2 -anc $reference_genome -GL 2 -nThreads 4 -minMapQ 30 -minQ 20
realSFS sfs_gbs_ICR_control_F2.saf.idx -maxIter 100 -P 4 > ICR.GBS.control.F2.sfs
