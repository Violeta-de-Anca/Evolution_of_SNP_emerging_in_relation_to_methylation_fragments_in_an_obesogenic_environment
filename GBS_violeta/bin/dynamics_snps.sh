#!/bin/bash -l
#SBATCH -A naiss2023-22-848
#SBATCH -p core -n 1
##SBATCH --mem=80gb
#SBATCH -t 10:00:00
#SBATCH -J snpeff
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/snpEff.on.parental.meth.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/log_files/snpEff.on.parental.meth.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load samtools/1.14
module load plink
module load BEDTools
module load bcftools
module load java/OpenJDK_17+35
module load vcftools BEDOPS
module load modkit
module load python3

working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta
parental_in=$working_dir/dynamics_SNPs/parental_methylated_win
new_filters=$parental_in/no_cluster_filters
pro_f2_15=$working_dir/dynamics_SNPs/progeny_f2.15/no_cluster_filters
pro_f2_16=$working_dir/dynamics_SNPs/progeny_f2.16/no_cluster_filters
pro_f2_14=$working_dir/dynamics_SNPs/progeny_f2.14/no_cluster_filters
pro_f2_2=$working_dir/dynamics_SNPs/progeny_f2.2/no_cluster_filters
pro_f2_1=$working_dir/dynamics_SNPs/progeny_f2.1/no_cluster_filters
pro_f2_13=$working_dir/dynamics_SNPs/progeny_f2.13/no_cluster_filters
pro_f1_10=$working_dir/dynamics_SNPs/progeny_f1.10/no_cluster_filters
cg_bed=/proj/naiss2024-23-57/reference_genomes/mus_musculus/cg_motif
reference=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref
snpsift=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar

# First filter parents by being homozygous on C in the meth regions

#file_list=($parental_in/*genotypes.gatk.filtered.by.own.win.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.genotypes.gatk.filtered.by.own.win.vcf}
#	cat $file | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "( REF = 'C' ) & isHom( GEN[0] )" > $x.homozygous_C.own.win.vcf
#done

#file_list=($new_filters/*genotypes.gatk.filtered.by.own.win.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes.gatk.filtered.by.own.win.vcf}
#        cat $file | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "( REF = 'C' ) & isHom( GEN[0] ) & isRef( GEN[0] )" > $new_filters/$x.homozygous_C.own.win.vcf
	#there are still ./. snps, so I am testing vcftools --vcf --recode --recode-INFO-all --max-missing --stdout
#	vcftools --vcf $x.homozygous_C.own.win.vcf --max-missing 1 --out $x.homozygous_C.own.win.nomissingdata --recode --recode-INFO-all
#done

#now that we have the homozygous sites for CC in fathers, convert that into a bed file so we can filter the sons
#file_list=($new_filters/*homozygous_C.own.win.nomissingdata.recode.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.homozygous_C.own.win.nomissingdata.recode.vcf}
#	convert2bed -i vcf < $file > $new_filters/$x.cc.in.meth.regions.bed
#done

#with the BED files from the fathers, let's filter the sons
#for F2_1
#for no_cluster_filters
#out_parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win
#parental=$out_parental/no_cluster_filters/F2_1.cc.in.meth.regions.bed
#file_list=($pro_f2_1/*genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
	#bgzip $file
	#tabix $file.gz
#	x=${name%.genotypes_refiltered_gatk.vcf.gz}
#	bcftools view -R $parental -o $pro_f2_1/$x.parental.CC.refilter.vcf $file
#done

#Now that the filtration is done, get rid of the ./. and filter only homozygous DD
#cd $pro_f2_1/no_cluster_filters
#file_list=($pro_f2_1/*.parental.CC.refilter.vcf)
#for file in "${file_list[@]}"; do
#       name=${file##*/}
#       x=${name%.parental.CC.refilter.vcf}
#       cat $file | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $x.DD.SNPs.refilter.parental.win.vcf
#       vcftools --vcf $x.DD.SNPs.refilter.parental.win.vcf --max-missing 1 --out $x.DD.SNPs.refilter.parental.win.nomissingdata --recode --recode-INFO-all
#done

#now for the final stats let's see how many heterozygous and homozygous we have
#for heterozygous
#cat F3_12.DD.SNPs.refilter.parental.win.vcf | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] ) & isVariant( GEN[0] )" --inverse > F3_12.heteroD.SNPs.refilter.parental.win.vcf
#for homozygous
#cat F3_12.DD.SNPs.refilter.parental.win.vcf | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > F3_12.homoDD.SNPs.refilter.parental.win.vcf


#for F2_14
#out_parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win
#parental=$out_parental/no_cluster_filters/F2_14.cc.in.meth.regions.bed
#file_list=($pro_f2_14/*genotypes_refiltered_gatk.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        bgzip $file
#        tabix $file.gz
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $parental -o $pro_f2_14/$x.parental.CC.refilter.vcf $file.gz
#done

#Now that the filtration is done, get rid of the ./. and filter only homozygous DD
#cd $pro_f2_14
#file_list=($pro_f2_14/*genotypes_refiltered_gatk.vcf.parental.CC.refilter.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.genotypes_refiltered_gatk.vcf.parental.CC.refilter.vcf}
#	cat $file | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $x.DD.SNPs.refilter.parental.win.vcf
#	vcftools --vcf $x.DD.SNPs.refilter.parental.win.vcf --max-missing 1 --out $x.DD.SNPs.refilter.parental.win.nomissingdata --recode --recode-INFO-all
#done

#only for F3_17
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.14/no_cluster_filters
#bcftools view -R $parental -o $pro_f2_14/no_cluster_filters/F3_17.parental.CC.refilter.vcf F3_17.genotypes_refiltered_gatk.vcf.gz
#cat $pro_f2_14/no_cluster_filters/F3_17.parental.CC.refilter.vcf | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > F3_17.DD.SNPs.refilter.parental.win.vcf


#for F2_15
#first filter parental
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/no_cluster_filters
#geno_parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/no_cluster_filters/C13F2_15.genotypes.gatk.filtered.by.own.win.vcf
#name=${geno_parental##*/}
#x=${name%.genotypes.gatk.filtered.by.own.win.vcf}
#cat $geno_parental | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "( REF = 'C' ) & isHom( GEN[0] ) & isRef( GEN[0] )" > /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win/no_cluster_filters/$x.homozygous_C.own.win.vcf
#vcftools --vcf $x.homozygous_C.own.win.vcf --max-missing 1 --out $x.homozygous_C.own.win.nomissingdata --recode --recode-INFO-all
#convert2bed -i vcf < C13F2_15.homozygous_C.own.win.nomissingdata.recode.vcf > $new_filters/F2_15.cc.in.meth.regions.bed

#out_parental=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/parental_methylated_win
#parental=$out_parental/no_cluster_filters/F2_15.cc.in.meth.regions.bed
#file_list=($pro_f2_15/no_cluster_filters/*genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $parental -o $pro_f2_15/no_cluster_filters/$x.parental.CC.refilter.vcf $file
#done

##Now that the filtration is done, get rid of the ./. and filter only homozygous DD
#cd $pro_f2_15
#file_list=($pro_f2_15/no_cluster_filters/*.parental.CC.refilter.vcf)
#for file in "${file_list[@]}"; do
#       name=${file##*/}
#       x=${name%.parental.CC.refilter.vcf}
#       cat $file | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $x.DD.SNPs.refilter.parental.win.vcf
       #vcftools --vcf $x.DD.SNPs.refilter.parental.win.vcf --max-missing 1 --out $x.DD.SNPs.refilter.parental.win.nomissingdata --recode --recode-INFO-all
#done


##########################################################################
# Let's repeat everything but this time more orderly
working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta
parental_in=$working_dir/dynamics_SNPs/parental_methylated_win
new_filters=$parental_in/no_cluster_filters
pro_f2_15=$working_dir/dynamics_SNPs/progeny_f2.15/no_cluster_filters
pro_f2_16=$working_dir/dynamics_SNPs/progeny_f2.16/no_cluster_filters
pro_f2_14=$working_dir/dynamics_SNPs/progeny_f2.14/no_cluster_filters
pro_f2_2=$working_dir/dynamics_SNPs/progeny_f2.2/no_cluster_filters
pro_f2_1=$working_dir/dynamics_SNPs/progeny_f2.1/no_cluster_filters
pro_f2_13=$working_dir/dynamics_SNPs/progeny_f2.13/no_cluster_filters
pro_f1_10=$working_dir/dynamics_SNPs/progeny_f1.10/no_cluster_filters
cg_bed=/proj/naiss2024-23-57/reference_genomes/mus_musculus/cg_motif
reference=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref
snpsift=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar
meth=$working_dir/verification.peaks.MACS3/intersect_win_per_indv/GBS_with_MACS3.win

# get all the methylated regions from the fathers ##
#first put everything is gz
#file_list=($new_filters/*genotypes_refiltered_gatk.vcf)
#for file in "${file_list[@]}"; do
#       bgzip $file
#       tabix $file.gz
#done

#####################################################################################

#Then subset the fathers to only their methylated regions
#file_list=($new_filters/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $meth/$x.GBS.with.MACS3.win.intersect.bed -o $new_filters/$x.refilt.own.win.vcf $file
#done

#######################################################################################
# Do the CG bed file
#modkit motif-bed $reference/mm39.fa CG 0 > $cg_bed/CG_motif_mm39.bed

#########################################################################################
# Now subset the CG from the fathers .refilt.own.win.vcf
#file_list=($new_filters/*.refilt.own.win.vcf)
#for file in "${file_list[@]}"; do
#	bgzip $file
#	tabix $file.gz
#done

#file_list=($new_filters/*.refilt.own.win.vcf.gz)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.refilt.own.win.vcf.gz}
#	bcftools view -R $cg_bed/CG_motif_mm39.bed -o $new_filters/$x.refilt.own.win.CG.vcf $file
#done

###########################################################################################
# Now subset homozygous CC in the fathers in CG sites
#file_list=($new_filters/*.refilt.own.win.CG.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.refilt.own.win.CG.vcf}
#	cat $file | java -jar $snpsift filter "( REF = 'C' ) & isHom( GEN[0] ) & isRef( GEN[0] )" > $new_filters/$x.CC.inCG.own.win.vcf
#done

# Now get out all the sites without info ./.
#file_list=($new_filters/*.CC.inCG.own.win.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.CC.inCG.own.win.vcf}
#	vcftools --vcf $file --max-missing 1 --out $new_filters/$x.CC.inCG.own.win.nomissingdata --recode --recode-INFO-all
#done

# Now convert the vcf files from the fathers into a bed file
#file_list=($new_filters/*.CC.inCG.own.win.nomissingdata.recode.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.CC.inCG.own.win.nomissingdata.recode.vcf}
#	convert2bed -i vcf < $file > $new_filters/$x.CC.inCG.own.win.nomissingdata.recode.bed
#done

##################################################################################

#Now is when we need to filter the sons of each parent
# F2_1

#file_list=($pro_f2_1/*genotypes_refiltered_gatk.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	bgzip $file
#	tabix $file.gz
#done

#file_list=($pro_f2_1/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#	x=${name%.genotypes_refiltered_gatk.vcf.gz}
#	bcftools view -R $new_filters/C13F2_1.CC.inCG.own.win.nomissingdata.recode.bed -o $pro_f2_1/$x.parental.CC.refilter.vcf $file
#done

#file_list=($pro_f2_1/*.parental.CC.refilter.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.parental.CC.refilter.vcf}
#	cat $file | java -jar $snpsift filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f2_1/$x.DD.SNPs.refilter.parental.win.vcf
#done

#file_list=($pro_f2_1/*.DD.SNPs.refilter.parental.win.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.DD.SNPs.refilter.parental.win.vcf}
#	cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" --inverse > $pro_f2_1/$x.heteroD.SNPs.refilter.parental.win.vcf
#	cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f2_1/$x.homoDD.SNPs.refilter.parental.win.vcf
#done

# F2_2

#file_list=($pro_f2_2/*genotypes_refiltered_gatk.vcf)
#for file in "${file_list[@]}"; do
#       bgzip $file
#       tabix $file.gz
#done

#file_list=($pro_f2_2/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.genotypes_refiltered_gatk.vcf.gz}
#	bcftools view -R $new_filters/C13F2_2.CC.inCG.own.win.nomissingdata.recode.bed -o $pro_f2_2/$x.parental.CC.refilter.vcf $file
#done

#file_list=($pro_f2_2/*.parental.CC.refilter.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.parental.CC.refilter.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f2_2/$x.DD.SNPs.refilter.parental.win.vcf
#done

#file_list=($pro_f2_2/*.DD.SNPs.refilter.parental.win.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.DD.SNPs.refilter.parental.win.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" --inverse > $pro_f2_2/$x.heteroD.SNPs.refilter.parental.win.vcf
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f2_2/$x.homoDD.SNPs.refilter.parental.win.vcf
#done

# F2_13

#file_list=($pro_f2_13/*genotypes_refiltered_gatk.vcf)
#for file in "${file_list[@]}"; do
#       bgzip $file
#       tabix $file.gz
#done

#file_list=($pro_f2_13/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_13.CC.inCG.own.win.nomissingdata.recode.bed -o $pro_f2_13/$x.parental.CC.refilter.vcf $file
#done

#file_list=($pro_f2_13/*.parental.CC.refilter.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.parental.CC.refilter.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f2_13/$x.DD.SNPs.refilter.parental.win.vcf
#done

#file_list=($pro_f2_13/*.DD.SNPs.refilter.parental.win.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.DD.SNPs.refilter.parental.win.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" --inverse > $pro_f2_13/$x.heteroD.SNPs.refilter.parental.win.vcf
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f2_13/$x.homoDD.SNPs.refilter.parental.win.vcf
#done

# F2_14

#file_list=($pro_f2_14/*genotypes_refiltered_gatk.vcf)
#for file in "${file_list[@]}"; do
#       bgzip $file
#       tabix $file.gz
#done

#file_list=($pro_f2_14/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_14.CC.inCG.own.win.nomissingdata.recode.bed -o $pro_f2_14/$x.parental.CC.refilter.vcf $file
#done

#file_list=($pro_f2_14/*.parental.CC.refilter.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.parental.CC.refilter.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f2_14/$x.DD.SNPs.refilter.parental.win.vcf
#done

#file_list=($pro_f2_14/*.DD.SNPs.refilter.parental.win.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.DD.SNPs.refilter.parental.win.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" --inverse > $pro_f2_14/$x.heteroD.SNPs.refilter.parental.win.vcf
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f2_14/$x.homoDD.SNPs.refilter.parental.win.vcf
#done

# F2_15

#file_list=($pro_f2_15/*genotypes_refiltered_gatk.vcf)
#for file in "${file_list[@]}"; do
#       bgzip $file
#       tabix $file.gz
#done

#file_list=($pro_f2_15/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_15.CC.inCG.own.win.nomissingdata.recode.bed -o $pro_f2_15/$x.parental.CC.refilter.vcf $file
#done

#file_list=($pro_f2_15/*.parental.CC.refilter.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.parental.CC.refilter.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f2_15/$x.DD.SNPs.refilter.parental.win.vcf
#done

#file_list=($pro_f2_15/*.DD.SNPs.refilter.parental.win.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.DD.SNPs.refilter.parental.win.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" --inverse > $pro_f2_15/$x.heteroD.SNPs.refilter.parental.win.vcf
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f2_15/$x.homoDD.SNPs.refilter.parental.win.vcf
#done

# F2_16

#file_list=($pro_f2_16/*genotypes_refiltered_gatk.vcf)
#for file in "${file_list[@]}"; do
#       bgzip $file
#       tabix $file.gz
#done

#file_list=($pro_f2_16/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_16.CC.inCG.own.win.nomissingdata.recode.bed -o $pro_f2_16/$x.parental.CC.refilter.vcf $file
#done

#file_list=($pro_f2_16/*.parental.CC.refilter.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.parental.CC.refilter.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f2_16/$x.DD.SNPs.refilter.parental.win.vcf
#done

#file_list=($pro_f2_16/*.DD.SNPs.refilter.parental.win.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.DD.SNPs.refilter.parental.win.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" --inverse > $pro_f2_16/$x.heteroD.SNPs.refilter.parental.win.vcf
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f2_16/$x.homoDD.SNPs.refilter.parental.win.vcf
#done

# F1_10

#file_list=($pro_f1_10/*genotypes_refiltered_gatk.vcf)
#for file in "${file_list[@]}"; do
#       bgzip $file
#       tabix $file.gz
#done

#file_list=($pro_f1_10/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F1_10.CC.inCG.own.win.vcf.CC.inCG.own.win.nomissingdata.recode.bed -o $pro_f1_10/$x.parental.CC.refilter.vcf $file
#done

#file_list=($pro_f1_10/*.parental.CC.refilter.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.parental.CC.refilter.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f1_10/$x.DD.SNPs.refilter.parental.win.vcf
#done

#file_list=($pro_f1_10/*.DD.SNPs.refilter.parental.win.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.DD.SNPs.refilter.parental.win.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" --inverse > $pro_f1_10/$x.heteroD.SNPs.refilter.parental.win.vcf
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f1_10/$x.homoDD.SNPs.refilter.parental.win.vcf
#done

##################################################################################################
##################################################################################################

#Let's exclude from all the parents the sites which already have a SNP
#file_list=($new_filters/*.refilt.own.win.vcf.gz)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.refilt.own.win.vcf.gz}
#	zcat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isRef( GEN[0] ) " > $new_filters/$x.ref.homo.own.win.refilt.vcf
#done

# take out all the missing data
#file_list=($new_filters/*.ref.homo.own.win.refilt.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.ref.homo.own.win.refilt.vcf}
#	vcftools --vcf $file --max-missing 1 --out $new_filters/$x.ref.homo.own.win.refilt.nomiss --recode --recode-INFO-all
#done

##############################################################################

# Get all the SNPs from all the sons
# F1_10
#sbatfile_list=($pro_f1_10/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.genotypes_refiltered_gatk.vcf.gz}
#	bcftools view -R $new_filters/C13F1_10.ref.homo.own.win.refilt.nomiss.recode.vcf -o $pro_f1_10/$x.parental.meth.win.refilt.vcf $file
#done

# F2_1
#file_list=($pro_f2_1/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_1.ref.homo.own.win.refilt.nomiss.recode.vcf -o $pro_f2_1/$x.parental.meth.win.refilt.vcf $file
#done

# F2_13
#file_list=($pro_f2_13/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_13.ref.homo.own.win.refilt.nomiss.recode.vcf -o $pro_f2_13/$x.parental.meth.win.refilt.vcf $file
#done

#F2_14
#file_list=($pro_f2_14/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_14.ref.homo.own.win.refilt.nomiss.recode.vcf -o $pro_f2_14/$x.parental.meth.win.refilt.vcf $file
#done

#F2_15
#file_list=($pro_f2_15/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_15.ref.homo.own.win.refilt.nomiss.recode.vcf -o $pro_f2_15/$x.parental.meth.win.refilt.vcf $file
#done

#F2_16
#file_list=($pro_f2_16/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_16.ref.homo.own.win.refilt.nomiss.recode.vcf -o $pro_f2_16/$x.parental.meth.win.refilt.vcf $file
#done

#F2_2
#file_list=($pro_f2_2/*.genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_2.ref.homo.own.win.refilt.nomiss.recode.vcf -o $pro_f2_2/$x.parental.meth.win.refilt.vcf $file
#done

#################################################################

#Now we need to extract all the aternative SNPs

#F1_10
#file_list=($pro_f1_10/*.parental.meth.win.refilt.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.parental.meth.win.refilt.vcf}
#	cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f1_10/$x.homo.SNPs.meth.win.father.vcf
#	cat $file | java -jar $snpsift filter "isVariant( GEN[0] )" > $pro_f1_10/$x.all.SNPs.meth.win.father.vcf
#done

#f2_1
#file_list=($pro_f2_1/*.parental.meth.win.refilt.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.parental.meth.win.refilt.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f2_1/$x.homo.SNPs.meth.win.father.vcf
#        cat $file | java -jar $snpsift filter "isVariant( GEN[0] )" > $pro_f2_1/$x.all.SNPs.meth.win.father.vcf
#done

#F2_13
#file_list=($pro_f2_13/*.parental.meth.win.refilt.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.parental.meth.win.refilt.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f2_13/$x.homo.SNPs.meth.win.father.vcf
#        cat $file | java -jar $snpsift filter "isVariant( GEN[0] )" > $pro_f2_13/$x.all.SNPs.meth.win.father.vcf
#done

#F2_14
#file_list=($pro_f2_14/*.parental.meth.win.refilt.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.parental.meth.win.refilt.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f2_14/$x.homo.SNPs.meth.win.father.vcf
#        cat $file | java -jar $snpsift filter "isVariant( GEN[0] )" > $pro_f2_14/$x.all.SNPs.meth.win.father.vcf
#done

#F2_15
#file_list=($pro_f2_15/*.parental.meth.win.refilt.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.parental.meth.win.refilt.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f2_15/$x.homo.SNPs.meth.win.father.vcf
#        cat $file | java -jar $snpsift filter "isVariant( GEN[0] )" > $pro_f2_15/$x.all.SNPs.meth.win.father.vcf
#done

#F2_16
#file_list=($pro_f2_16/*.parental.meth.win.refilt.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.parental.meth.win.refilt.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f2_16/$x.homo.SNPs.meth.win.father.vcf
#        cat $file | java -jar $snpsift filter "isVariant( GEN[0] )" > $pro_f2_16/$x.all.SNPs.meth.win.father.vcf
#done

#F2_2
#file_list=($pro_f2_2/*.parental.meth.win.refilt.vcf)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.parental.meth.win.refilt.vcf}
#        cat $file | java -jar $snpsift filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $pro_f2_2/$x.homo.SNPs.meth.win.father.vcf
#        cat $file | java -jar $snpsift filter "isVariant( GEN[0] )" > $pro_f2_2/$x.all.SNPs.meth.win.father.vcf
#done

###################################################################################

#Merge and extract all the SNPs in a table
#sort and index
#cd $pro_f2_2
#file_list=(*.all.SNPs.meth.win.father.vcf.gz)
#for i in "${file_list[@]}"; do
#	name=${i##*/}
#	x=${name%.all.SNPs.meth.win.father.vcf}
#	#bcftools sort --write-index -o $x.all.SNPs.meth.win.father.sorted.vcf
#	#bgzip $i
#	tabix $i
#done

#cd $pro_f2_16
#file_list=(*.all.SNPs.meth.win.father.vcf.gz)
#for i in "${file_list[@]}"; do
#        name=${i##*/}
#        x=${name%.all.SNPs.meth.win.father.vcf}
	#bcftools sort --write-index -o $x.all.SNPs.meth.win.father.sorted.vcf
	#bgzip $i
	#tabix $i
#done

#cd $pro_f2_15
#file_list=(*.all.SNPs.meth.win.father.vcf.gz)
#for i in "${file_list[@]}"; do
#        name=${i##*/}
#        x=${name%.all.SNPs.meth.win.father.vcf}
        #bcftools sort --write-index -o $x.all.SNPs.meth.win.father.sorted.vcf
	#bgzip $i
	#tabix $i
#done

#cd $pro_f2_14
#file_list=(*.all.SNPs.meth.win.father.vcf.gz)
#for i in "${file_list[@]}"; do
#        name=${i##*/}
#        x=${name%.all.SNPs.meth.win.father.vcf}
        #bcftools sort --write-index -o $x.all.SNPs.meth.win.father.sorted.vcf
	#bgzip $i
	#tabix $i
#done

#cd $pro_f2_13
#file_list=(*.all.SNPs.meth.win.father.vcf.gz)
#for i in "${file_list[@]}"; do
#        name=${i##*/}
#        x=${name%.all.SNPs.meth.win.father.vcf}
        #bcftools sort --write-index -o $x.all.SNPs.meth.win.father.sorted.vcf
	#bgzip $i
	#tabix $i
#done

#cd $pro_f2_1
#file_list=(*.all.SNPs.meth.win.father.vcf.gz)
#for i in "${file_list[@]}"; do
#        name=${i##*/}
#        x=${name%.all.SNPs.meth.win.father.vcf}
        #bcftools sort --write-index -o $x.all.SNPs.meth.win.father.sorted.vcf
	#bgzip $i
        #tabix $i
#done

#cd $pro_f1_10
#file_list=(*.all.SNPs.meth.win.father.vcf.gz)
#for i in "${file_list[@]}"; do
#        name=${i##*/}
#        x=${name%.all.SNPs.meth.win.father.vcf}
        #bcftools sort --write-index -o $x.all.SNPs.meth.win.father.sorted.vcf
	#bgzip $i
	#tabix $i
#done

##############################################################################

#merge
#cd $pro_f2_2
#vcf-merge C13F3_1.all.SNPs.meth.win.father.vcf.gz C13F3_3.all.SNPs.meth.win.father.vcf.gz C13F3_4.all.SNPs.meth.win.father.vcf.gz C13F3_5.all.SNPs.meth.win.father.vcf.gz F3_2.all.SNPs.meth.win.father.vcf.gz > merged.SNPs.sons.from.F2_2.methwind.vcf

#cd $pro_f2_16
#vcf-merge C13F3_18.all.SNPs.meth.win.father.vcf.gz F3_22.all.SNPs.meth.win.father.vcf.gz F3_27.all.SNPs.meth.win.father.vcf.gz > merged.SNPs.sons.from.F2_16.methwind.vcf

#cd $pro_f2_15
#vcf-merge C13F3_20.all.SNPs.meth.win.father.vcf.gz F3_19.all.SNPs.meth.win.father.vcf.gz F3_21.all.SNPs.meth.win.father.vcf.gz > merged.SNPs.sons.from.F2_15.methwind.vcf

#cd $pro_f2_14
#vcf-merge F3_17.all.SNPs.meth.win.father.vcf.gz F3_24.all.SNPs.meth.win.father.vcf.gz F3_25.all.SNPs.meth.win.father.vcf.gz > merged.SNPs.sons.from.F2_14.methwind.vcf

#cd $pro_f2_13
#vcf-merge C13F3_15.all.SNPs.meth.win.father.vcf.gz C13F3_16.all.SNPs.meth.win.father.vcf.gz > merged.SNPs.sons.from.F2_13.methwind.vcf

#cd $pro_f2_1
#vcf-merge C13F3_10.all.SNPs.meth.win.father.vcf.gz F3_11.all.SNPs.meth.win.father.vcf.gz F3_12.all.SNPs.meth.win.father.vcf.gz F3_7.all.SNPs.meth.win.father.vcf.gz F3_8.all.SNPs.meth.win.father.vcf.gz F3_9.all.SNPs.meth.win.father.vcf.gz > merged.SNPs.sons.from.F2_1.methwind.vcf

#cd $pro_f1_10
#vcf-merge C13F2_12.all.SNPs.meth.win.father.vcf.gz C13F2_4.all.SNPs.meth.win.father.vcf.gz C13F2_6.all.SNPs.meth.win.father.vcf.gz > merged.SNPs.sons.from.F1_10.methwind.vcf

######################################################################################

# Index the vcf files to be able to merge them!

#while read line; do
#	bgzip $line
#	tabix $line.gz
#done < /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/list

######################################################################################
gen_out=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs

#Merge all the CG SNPs
#vcf-merge ${pro_f2_2}/C13F3_1.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_2}/C13F3_3.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_2}/C13F3_4.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_2}/C13F3_5.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_2}/F3_2.DD.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.CG.SNPs.sons.from.F2_2.methwind.vcf
#vcf-merge ${pro_f2_16}/C13F3_18.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_16}/F3_22.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_16}/F3_27.DD.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.CG.SNPs.sons.from.F2_16.methwind.vcf
#vcf-merge ${pro_f2_15}/C13F3_20.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_15}/F3_19.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_15}/F3_21.DD.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.CG.SNPs.sons.from.F2_15.methwind.vcf
#vcf-merge ${pro_f2_14}/F3_17.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_14}/F3_24.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_14}/F3_25.DD.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.CG.SNPs.sons.from.F2_14.methwind.vcf
#vcf-merge ${pro_f2_13}/C13F3_15.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_13}/C13F3_16.DD.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.CG.SNPs.sons.from.F2_13.methwind.vcf
#vcf-merge ${pro_f2_1}/C13F3_10.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_1}/C13F3_11.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_1}/C13F3_12.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_1}/C13F3_7.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_1}/C13F3_8.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f2_1}/C13F3_9.DD.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.CG.SNPs.sons.from.F2_1.methwind.vcf
#vcf-merge ${pro_f1_10}/C13F2_12.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f1_10}/C13F2_4.DD.SNPs.refilter.parental.win.vcf.gz ${pro_f1_10}/C13F2_6.DD.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.CG.SNPs.sons.from.F1_10.methwind.vcf

######################################################################################

# Get only the SNPs, NOT APPLICABLE, THEY HAVE BEEN ALREADY SUBSETED
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs
#cat ${gen_out}/merged.CG.SNPs.sons.from.F1_10.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.CG.SNPs.meth.win.father.F1_10.vcf
#cat ${gen_out}/merged.CG.SNPs.sons.from.F2_1.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.CG.SNPs.meth.win.father.F2_1.vcf
#cat ${gen_out}/merged.CG.SNPs.sons.from.F2_2.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.CG.SNPs.meth.win.father.F2_2.vcf
#cat ${gen_out}/merged.CG.SNPs.sons.from.F2_13.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.CG.SNPs.meth.win.father.F2_13.vcf
#cat ${gen_out}/merged.CG.SNPs.sons.from.F2_14.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.CG.SNPs.meth.win.father.F2_14.vcf
#cat ${gen_out}/merged.CG.SNPs.sons.from.F2_15.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.CG.SNPs.meth.win.father.F2_15.vcf
#cat ${gen_out}/merged.CG.SNPs.sons.from.F2_16.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.CG.SNPs.meth.win.father.F2_16.vcf

######################################################################################

# Now extract all SNPs in a table
#cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs

#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.SNPs.sons.from.F1_10.methwind.vcf > table.SNPs.sons.from.F1_10.methwind.txt
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.SNPs.sons.from.F2_1.methwind.vcf > table.SNPs.sons.from.F2_1.methwind.txt
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.SNPs.sons.from.F2_13.methwind.vcf > table.SNPs.sons.from.F2_13.methwind.txt
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.SNPs.sons.from.F2_14.methwind.vcf > table.SNPs.sons.from.F2_14.methwind.txt
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.SNPs.sons.from.F2_15.methwind.vcf > table.SNPs.sons.from.F2_15.methwind.txt
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.SNPs.sons.from.F2_16.methwind.vcf > table.SNPs.sons.from.F2_16.methwind.txt
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.SNPs.sons.from.F2_2.methwind.vcf > table.SNPs.sons.from.F2_2.methwind.txt

#now for only CG
cd /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.CG.SNPs.sons.from.F1_10.methwind.vcf > table.CG.SNPs.sons.from.F1_10.methwind.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.CG.SNPs.sons.from.F2_1.methwind.vcf > table.CG.SNPs.sons.from.F2_1.methwind.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.CG.SNPs.sons.from.F2_2.methwind.vcf > table.CG.SNPs.sons.from.F2_2.methwind.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.CG.SNPs.sons.from.F2_13.methwind.vcf > table.CG.SNPs.sons.from.F2_13.methwind.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.CG.SNPs.sons.from.F2_14.methwind.vcf > table.CG.SNPs.sons.from.F2_14.methwind.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.CG.SNPs.sons.from.F2_15.methwind.vcf > table.CG.SNPs.sons.from.F2_15.methwind.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' merged.CG.SNPs.sons.from.F2_16.methwind.vcf > table.CG.SNPs.sons.from.F2_16.methwind.txt

######################################################################################

#take out all the low quality SNPs and make sure there are no homozygous alleles #
file_list=(/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/table.CG.SNPs.sons.from.*.txt)
for i in "${file_list[@]}"; do
	rm ${i}.nobad
	grep -v '*' ${i} > ${snp_folder}/${i}.nobad
done

# Sum all the SNPs
bin=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin
snp_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs
meth_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/merged
#python3 sum_SNPs.py $snp_folder/table.SNPs.sons.from.F1_10.methwind.txt.nobad $meth_folder/C13F1_10.featureCounts.bed > $snp_folder/finaltable.SNPs.sons.from.F1_10.methwind.txt
#python3 sum_SNPs.py $snp_folder/table.SNPs.sons.from.F2_1.methwind.txt.nobad $meth_folder/C13F2_1.featureCounts.bed > $snp_folder/finaltable.SNPs.sons.from.F2_1.methwind.txt
#python3 sum_SNPs.py $snp_folder/table.SNPs.sons.from.F2_13.methwind.txt.nobad $meth_folder/C13F2_13.featureCounts.bed > $snp_folder/finaltable.SNPs.sons.from.F2_13.methwind.txt
#python3 sum_SNPs.py $snp_folder/table.SNPs.sons.from.F2_14.methwind.txt.nobad $meth_folder/C13F2_14.featureCounts.bed > $snp_folder/finaltable.SNPs.sons.from.F2_14.methwind.txt
#python3 sum_SNPs.py $snp_folder/table.SNPs.sons.from.F2_15.methwind.txt.nobad $meth_folder/C13F2_15.featureCounts.bed > $snp_folder/finaltable.SNPs.sons.from.F2_15.methwind.txt
#python3 sum_SNPs.py $snp_folder/table.SNPs.sons.from.F2_16.methwind.txt.nobad $meth_folder/C13F2_16.featureCounts.bed > $snp_folder/finaltable.SNPs.sons.from.F2_16.methwind.txt
#python3 sum_SNPs.py $snp_folder/table.SNPs.sons.from.F2_2.methwind.txt.nobad $meth_folder/C13F2_2.featureCounts.bed > $snp_folder/finaltable.SNPs.sons.from.F2_2.methwind.txt

python3 $bin/sum_SNPs.py ${snp_folder}/table.CG.SNPs.sons.from.F1_10.methwind.txt.nobad $meth_folder/C13F1_10.featureCounts.bed > $snp_folder/finaltable.CG.SNPs.sons.from.F1_10.methwind.txt
python3 $bin/sum_SNPs.py ${snp_folder}/table.CG.SNPs.sons.from.F2_1.methwind.txt.nobad $meth_folder/C13F2_1.featureCounts.bed > $snp_folder/finaltable.CG.SNPs.sons.from.F2_1.methwind.txt
python3 $bin/sum_SNPs.py ${snp_folder}/table.CG.SNPs.sons.from.F2_2.methwind.txt.nobad $meth_folder/C13F2_2.featureCounts.bed > $snp_folder/finaltable.CG.SNPs.sons.from.F2_2.methwind.txt
python3 $bin/sum_SNPs.py ${snp_folder}/table.CG.SNPs.sons.from.F2_13.methwind.txt.nobad $meth_folder/C13F2_13.featureCounts.bed > $snp_folder/finaltable.CG.SNPs.sons.from.F2_13.methwind.txt
python3 $bin/sum_SNPs.py ${snp_folder}/table.CG.SNPs.sons.from.F2_14.methwind.txt.nobad $meth_folder/C13F2_14.featureCounts.bed > $snp_folder/finaltable.CG.SNPs.sons.from.F2_14.methwind.txt
python3 $bin/sum_SNPs.py ${snp_folder}/table.CG.SNPs.sons.from.F2_15.methwind.txt.nobad $meth_folder/C13F2_15.featureCounts.bed > $snp_folder/finaltable.CG.SNPs.sons.from.F2_15.methwind.txt
python3 $bin/sum_SNPs.py ${snp_folder}/table.CG.SNPs.sons.from.F2_16.methwind.txt.nobad $meth_folder/C13F2_16.featureCounts.bed > $snp_folder/finaltable.CG.SNPs.sons.from.F2_16.methwind.txt

#######################################################################################
# Do the AG bed file
#modkit motif-bed $reference/mm39.fa AG 0 > $cg_bed/AG_motif_mm39.bed

#########################################################################################

#######################################################################################
# Do the TG bed file
#modkit motif-bed $reference/mm39.fa TG 0 > $cg_bed/TG_motif_mm39.bed

#########################################################################################

#######################################################################################
# Do the GG bed file
#modkit motif-bed $reference/mm39.fa GG 0 > $cg_bed/GG_motif_mm39.bed

#########################################################################################

# Now subset the AG/TG/GG from the fathers .refilt.own.win.vcf
#file_list=($new_filters/*.refilt.own.win.vcf.gz)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.refilt.own.win.vcf.gz}
#	bcftools view -R $cg_bed/AG_motif_mm39.bed -o $new_filters/$x.refilt.own.win.AG.vcf $file
#	bcftools view -R $cg_bed/TG_motif_mm39.bed -o $new_filters/$x.refilt.own.win.TG.vcf $file
#	bcftools view -R $cg_bed/GG_motif_mm39.bed -o $new_filters/$x.refilt.own.win.GG.vcf $file
#done

##########################################################################################

# Now subset homozygous AA/TT/GG in the fathers in AG/TG/GG sites
#file_list=($new_filters/*.refilt.own.win.AG.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.refilt.own.win.AG.vcf}
#	cat $file | java -jar $snpsift filter "( REF = 'A' ) & isHom( GEN[0] ) & isRef( GEN[0] )" > $new_filters/$x.AA.inAG.own.win.vcf
#done
#file_list=($new_filters/*.refilt.own.win.GG.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.refilt.own.win.GG.vcf}
#	cat $file | java -jar $snpsift filter "( REF = 'G' ) & isHom( GEN[0] ) & isRef( GEN[0] )" > $new_filters/$x.GG.inGG.own.win.vcf
#done
#file_list=($new_filters/*.refilt.own.win.TG.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.refilt.own.win.TG.vcf}
#	cat $file | java -jar $snpsift filter "( REF = 'T' ) & isHom( GEN[0] ) & isRef( GEN[0] )" > $new_filters/$x.TT.inTG.own.win.vcf
#done

#############################################################################################

# take out the missing values
#file_list=($new_filters/*.TT.inTG.own.win.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.TT.inTG.own.win.vcf}
#	vcftools --vcf $file --max-missing 1 --out $new_filters/${x}_TT.own.win.nomissingdata --recode --recode-INFO-all
#done
#file_list=($new_filters/*.AA.inAG.own.win.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.AA.inAG.own.win.vcf}
#	vcftools --vcf $file --max-missing 1 --out $new_filters/${x}_AA.own.win.nomissingdata --recode --recode-INFO-all
#done
#file_list=($new_filters/*.GG.inGG.own.win.vcf)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.GG.inGG.own.win.vcf}
#	vcftools --vcf $file --max-missing 1 --out $new_filters/${x}_GG.own.win.nomissingdata --recode --recode-INFO-all
#done

###########################################################################

# convert to bed file the vcf files
#file_list=($new_filters/*_GG.own.win.nomissingdata.recode.vcf)
#for file in "${file_list[@]}"; do
#       name=${file##*/}
#       x=${name%_GG.own.win.nomissingdata.recode.vcf}
#       convert2bed -i vcf < $file > $new_filters/${x}.GG.in.meth.regions.bed
#done
#file_list=($new_filters/*_TT.own.win.nomissingdata.recode.vcf)
#for file in "${file_list[@]}"; do
#       name=${file##*/}
#       x=${name%_TT.own.win.nomissingdata.recode.vcf}
#       convert2bed -i vcf < $file > $new_filters/$x.TT.in.meth.regions.bed
#done
#file_list=($new_filters/*_AA.own.win.nomissingdata.recode.vcf)
#for file in "${file_list[@]}"; do
#       name=${file##*/}
#       x=${name%_AA.own.win.nomissingdata.recode.vcf}
#       convert2bed -i vcf < $file > $new_filters/$x.AA.in.meth.regions.bed
#done

#################################################################

# filter the sons
#for F2_1
#file_list=($pro_f2_1/*genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.genotypes_refiltered_gatk.vcf.gz}
#	bcftools view -R $new_filters/C13F2_1_AA.own.win.nomissingdata.recode.vcf -o $pro_f2_1/$x.parental.AA.refilter.vcf $file
#	bcftools view -R $new_filters/C13F2_1_TT.own.win.nomissingdata.recode.vcf -o $pro_f2_1/$x.parental.TT.refilter.vcf $file
#	bcftools view -R $new_filters/C13F2_1_GG.own.win.nomissingdata.recode.vcf -o $pro_f2_1/$x.parental.GG.refilter.vcf $file
#done
#for F2_2
#file_list=($pro_f2_2/*genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_2_AA.own.win.nomissingdata.recode.vcf -o $pro_f2_2/$x.parental.AA.refilter.vcf $file
#        bcftools view -R $new_filters/C13F2_2_TT.own.win.nomissingdata.recode.vcf -o $pro_f2_2/$x.parental.TT.refilter.vcf $file
#        bcftools view -R $new_filters/C13F2_2_GG.own.win.nomissingdata.recode.vcf -o $pro_f2_2/$x.parental.GG.refilter.vcf $file
#done
#for F2_13
#file_list=($pro_f2_13/*genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_13_AA.own.win.nomissingdata.recode.vcf -o $pro_f2_13/$x.parental.AA.refilter.vcf $file
#        bcftools view -R $new_filters/C13F2_13_TT.own.win.nomissingdata.recode.vcf -o $pro_f2_13/$x.parental.TT.refilter.vcf $file
#        bcftools view -R $new_filters/C13F2_13_GG.own.win.nomissingdata.recode.vcf -o $pro_f2_13/$x.parental.GG.refilter.vcf $file
#done
#for F2_14
#file_list=($pro_f2_14/*genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_14_AA.own.win.nomissingdata.recode.vcf -o $pro_f2_14/$x.parental.AA.refilter.vcf $file
#        bcftools view -R $new_filters/C13F2_14_TT.own.win.nomissingdata.recode.vcf -o $pro_f2_14/$x.parental.TT.refilter.vcf $file
#        bcftools view -R $new_filters/C13F2_14_GG.own.win.nomissingdata.recode.vcf -o $pro_f2_14/$x.parental.GG.refilter.vcf $file
#done
#for F2_15
#file_list=($pro_f2_15/*genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_15_AA.own.win.nomissingdata.recode.vcf -o $pro_f2_15/$x.parental.AA.refilter.vcf $file
#        bcftools view -R $new_filters/C13F2_15_TT.own.win.nomissingdata.recode.vcf -o $pro_f2_15/$x.parental.TT.refilter.vcf $file
#        bcftools view -R $new_filters/C13F2_15_GG.own.win.nomissingdata.recode.vcf -o $pro_f2_15/$x.parental.GG.refilter.vcf $file
#done
#for F2_16
#file_list=($pro_f2_16/*genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F2_16_AA.own.win.nomissingdata.recode.vcf -o $pro_f2_16/$x.parental.AA.refilter.vcf $file
#        bcftools view -R $new_filters/C13F2_16_TT.own.win.nomissingdata.recode.vcf -o $pro_f2_16/$x.parental.TT.refilter.vcf $file
#        bcftools view -R $new_filters/C13F2_16_GG.own.win.nomissingdata.recode.vcf -o $pro_f2_16/$x.parental.GG.refilter.vcf $file
#done
#for F1_10
#file_list=($pro_f1_10/*genotypes_refiltered_gatk.vcf.gz)
#for file in "${file_list[@]}"; do
#        name=${file##*/}
#        x=${name%.genotypes_refiltered_gatk.vcf.gz}
#        bcftools view -R $new_filters/C13F1_10_AA.own.win.nomissingdata.recode.vcf -o $pro_f1_10/$x.parental.AA.refilter.vcf $file
#        bcftools view -R $new_filters/C13F1_10_TT.own.win.nomissingdata.recode.vcf -o $pro_f1_10/$x.parental.TT.refilter.vcf $file
#        bcftools view -R $new_filters/C13F1_10_GG.own.win.nomissingdata.recode.vcf -o $pro_f1_10/$x.parental.GG.refilter.vcf $file
#done

##############################################################################################

#Now that the filtration is done, get rid of the ./. and filter only homozygous SNPs
#for f2_1
#for suffix in AA GG TT; do
#	find /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.1/no_cluster_filters -type f -name "*.parental.${suffix}.refilter.vcf" | while read -r file; do
#		name=${file##*/};
#		x=${name%.parental.${suffix}.refilter.vcf};
#		cat ${file} | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $f2_1/${x}.${suffix}.SNPs.refilter.parental.win.vcf ;
#	done;
#done
#for f2_2
#for suffix in AA GG TT; do
#	find /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.2/no_cluster_filters -type f -name "*.parental.${suffix}.refilter.vcf" | while read -r file; do
#		name=${file##*/};
#		x=${name%.parental.${suffix}.refilter.vcf};
#		echo $file
#		cat ${file} | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f2_2/${x}.${suffix}.SNPs.refilter.parental.win.vcf ;
#	done;
#done
#for f2_13
#for suffix in AA GG TT; do
#        find /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.13/no_cluster_filters -type f -name "*.parental.${suffix}.refilter.vcf" | while read -r file; do
#                name=${file##*/};
#		echo $file
#                x=${name%.parental.${suffix}.refilter.vcf};
#                cat ${file} | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f2_13/${x}.${suffix}.SNPs.refilter.parental.win.vcf ;
#        done;
#done
#for f2_14
#for suffix in AA GG TT; do
#        find /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.14/no_cluster_filters -type f -name "*.parental.${suffix}.refilter.vcf" | while read -r file; do
#                name=${file##*/};
#		echo $file
#                x=${name%.parental.${suffix}.refilter.vcf};
#                cat ${file} | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f2_14/${x}.${suffix}.SNPs.refilter.parental.win.vcf ;
#        done;
#done
#for f2_15
#for suffix in AA GG TT; do
#        find /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.15/no_cluster_filters -type f -name "*.parental.${suffix}.refilter.vcf" | while read -r file; do
#                name=${file##*/};
#		echo $file
#                x=${name%.parental.${suffix}.refilter.vcf};
#                cat ${file} | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f2_15/${x}.${suffix}.SNPs.refilter.parental.win.vcf ;
#        done;
#done
#for f2_16
#for suffix in AA GG TT; do
#        find /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f2.16/no_cluster_filters -type f -name "*.parental.${suffix}.refilter.vcf" | while read -r file; do
#                name=${file##*/};
#		echo $file
#                x=${name%.parental.${suffix}.refilter.vcf};
#                cat ${file} | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f2_16/${x}.${suffix}.SNPs.refilter.parental.win.vcf ;
#        done;
#done

#for f1_10
#for suffix in AA GG TT; do
#	find /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/progeny_f1.10/no_cluster_filters -type f -name "*.parental.${suffix}.refilter.vcf" | while read -r file; do
#		name=${file##*/};
#		x=${name%.parental.${suffix}.refilter.vcf};
#		cat ${file} | java -jar /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin/snpEff/SnpSift.jar filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > $pro_f1_10/${x}.${suffix}.SNPs.refilter.parental.win.vcf ;
#	done
#done

##################################################################################

# compressed and index
#while read line; do
#       bgzip $line
#       tabix $line.gz
#done < /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs/list

##################################################################################

# Merged files
#for suffix in AA GG TT; do
	#vcf-merge ${pro_f2_2}/C13F3_1.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_2}/C13F3_3.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_2}/C13F3_4.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_2}/C13F3_5.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_2}/F3_2.${suffix}.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.${suffix}.SNPs.sons.from.F2_2.methwind.vcf
	#vcf-merge ${pro_f2_16}/C13F3_18.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_16}/F3_22.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_16}/F3_27.${suffix}.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.${suffix}.SNPs.sons.from.F2_16.methwind.vcf
	#vcf-merge ${pro_f2_15}/C13F3_20.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_15}/F3_19.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_15}/F3_21.${suffix}.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.${suffix}.SNPs.sons.from.F2_15.methwind.vcf
	#vcf-merge ${pro_f2_14}/F3_17.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_14}/F3_24.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_14}/F3_25.${suffix}.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.${suffix}.SNPs.sons.from.F2_14.methwind.vcf
	#vcf-merge ${pro_f2_13}/C13F3_15.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_13}/C13F3_16.${suffix}.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.${suffix}.SNPs.sons.from.F2_13.methwind.vcf
	#vcf-merge ${pro_f2_1}/C13F3_10.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_1}/F3_11.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_1}/F3_12.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_1}/F3_7.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_1}/F3_8.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f2_1}/F3_9.${suffix}.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.${suffix}.SNPs.sons.from.F2_1.methwind.vcf
	#vcf-merge ${pro_f1_10}/C13F2_12.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f1_10}/C13F2_4.${suffix}.SNPs.refilter.parental.win.vcf.gz ${pro_f1_10}/C13F2_6.${suffix}.SNPs.refilter.parental.win.vcf.gz > ${gen_out}/merged.${suffix}.SNPs.sons.from.F1_10.methwind.vcf
#done

####################################################################################

# Get only the SNPs
#for suffix in AA GG TT; do
#	cat ${gen_out}/merged.${suffix}.SNPs.sons.from.F1_10.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.${suffix}.SNPs.meth.win.father.F1_10.vcf
#	cat ${gen_out}/merged.${suffix}.SNPs.sons.from.F2_1.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.${suffix}.SNPs.meth.win.father.F2_1.vcf
#	cat ${gen_out}/merged.${suffix}.SNPs.sons.from.F2_2.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.${suffix}.SNPs.meth.win.father.F2_2.vcf
#	cat ${gen_out}/merged.${suffix}.SNPs.sons.from.F2_13.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.${suffix}.SNPs.meth.win.father.F2_13.vcf
#	cat ${gen_out}/merged.${suffix}.SNPs.sons.from.F2_14.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.${suffix}.SNPs.meth.win.father.F2_14.vcf
#	cat ${gen_out}/merged.${suffix}.SNPs.sons.from.F2_15.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.${suffix}.SNPs.meth.win.father.F2_15.vcf
#	cat ${gen_out}/merged.${suffix}.SNPs.sons.from.F2_16.methwind.vcf | java -jar ${snpsift} filter "isHom( GEN[0] )  & isRef( GEN[0] ) " --inverse > ${gen_out}/merged.all.${suffix}.SNPs.meth.win.father.F2_16.vcf
#done

######################################################################################

bin=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/bin
snp_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/dynamics_SNPs
meth_folder=/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/merged

# Extract all the SNPs into a table
for suffix in AA GG TT; do
	echo ${snp_folder}/merged.${suffix}.SNPs.meth.win.father.F1_10.vcf
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' ${snp_folder}/merged.${suffix}.SNPs.sons.from.F1_10.methwind.vcf > ${gen_out}/table.${suffix}.SNPs.sons.from.F1_10.methwind.txt
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' ${snp_folder}/merged.${suffix}.SNPs.sons.from.F2_1.methwind.vcf > ${gen_out}/table.${suffix}.SNPs.sons.from.F2_1.methwind.txt
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' ${snp_folder}/merged.${suffix}.SNPs.sons.from.F2_2.methwind.vcf > ${gen_out}/table.${suffix}.SNPs.sons.from.F2_2.methwind.txt
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' ${snp_folder}/merged.${suffix}.SNPs.sons.from.F2_13.methwind.vcf > ${gen_out}/table.${suffix}.SNPs.sons.from.F2_13.methwind.txt
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' ${snp_folder}/merged.${suffix}.SNPs.sons.from.F2_14.methwind.vcf > ${gen_out}/table.${suffix}.SNPs.sons.from.F2_14.methwind.txt
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' ${snp_folder}/merged.${suffix}.SNPs.sons.from.F2_15.methwind.vcf > ${gen_out}/table.${suffix}.SNPs.sons.from.F2_15.methwind.txt
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t]\n' ${snp_folder}/merged.${suffix}.SNPs.sons.from.F2_16.methwind.vcf > ${gen_out}/table.${suffix}.SNPs.sons.from.F2_16.methwind.txt
done

#######################################################################################3

#Sum SNPs
for suffix in AA GG TT; do
	python3 $bin/sum_SNPs.py ${snp_folder}/table.${suffix}.SNPs.sons.from.F1_10.methwind.txt $meth_folder/C13F1_10.featureCounts.bed > $snp_folder/finaltable.${suffix}.SNPs.sons.from.F1_10.methwind.txt
	python3 $bin/sum_SNPs.py ${snp_folder}/table.${suffix}.SNPs.sons.from.F2_1.methwind.txt $meth_folder/C13F2_1.featureCounts.bed > $snp_folder/finaltable.${suffix}.SNPs.sons.from.F2_1.methwind.txt
	python3 $bin/sum_SNPs.py ${snp_folder}/table.${suffix}.SNPs.sons.from.F2_2.methwind.txt $meth_folder/C13F2_2.featureCounts.bed > $snp_folder/finaltable.${suffix}.SNPs.sons.from.F2_2.methwind.txt
	python3 $bin/sum_SNPs.py ${snp_folder}/table.${suffix}.SNPs.sons.from.F2_13.methwind.txt $meth_folder/C13F2_13.featureCounts.bed > $snp_folder/finaltable.${suffix}.SNPs.sons.from.F2_13.methwind.txt
	python3 $bin/sum_SNPs.py ${snp_folder}/table.${suffix}.SNPs.sons.from.F2_14.methwind.txt $meth_folder/C13F2_14.featureCounts.bed > $snp_folder/finaltable.${suffix}.SNPs.sons.from.F2_14.methwind.txt
	python3 $bin/sum_SNPs.py ${snp_folder}/table.${suffix}.SNPs.sons.from.F2_15.methwind.txt $meth_folder/C13F2_15.featureCounts.bed > $snp_folder/finaltable.${suffix}.SNPs.sons.from.F2_15.methwind.txt
	python3 $bin/sum_SNPs.py ${snp_folder}/table.${suffix}.SNPs.sons.from.F2_16.methwind.txt $meth_folder/C13F2_16.featureCounts.bed > $snp_folder/finaltable.${suffix}.SNPs.sons.from.F2_16.methwind.txt
done
