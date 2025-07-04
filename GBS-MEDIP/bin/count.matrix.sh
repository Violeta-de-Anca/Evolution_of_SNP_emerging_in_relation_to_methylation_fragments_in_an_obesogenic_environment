#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 6-00:00:00
#SBATCH -J countmatrix
#SBATCH --error /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/CountMatrix.err
#SBATCH --output /proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/log_files/CountMatrix.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

#load modules
module load bioinfo-tools
#module load R_packages/4.2.1
module load BEDTools
module load bcftools

working_dir=/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP
gbs-medip=/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/gbs-medip

file_list=($gbs-medip/C13F2*_Mouse.unique.MQ20.bam)
for file in "${file_list[@]}"; do
	name=${file##*/}
	x=${name%_Mouse.unique.MQ20.bam}
	bedtools intersect -a $file -c -b /proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/windows.MACS3.veri.GBS.hyper-hypo.bed > 
done
#R code starts here
Rscript - <<EOF

#Capture start time
start_time <- Sys.time()

#Load libraries
#BiocManager::install("MEDIPS", lib="/home/viole/")
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm39", lib="/home/viole/", force = TRUE)
#install.packages("MEDIPS", repos='http://cran.us.r-project.org', lib="/home/viole/")
#install.packages("BSgenome.Mmusculus.UCSC.mm39", repos='http://cran.us.r-project.org', lib="/home/viole/")
#BiocManager::install("BSgenome", lib="/home/viole/", force = TRUE)
BiocManager::install("BiocGenerics", lib="/home/viole/", force = TRUE)
library(BiocGenerics, lib.loc="/home/viole/")
library(BSgenome, lib.loc="/home/viole/")
library(MEDIPS, lib.loc="/home/viole/")
library(BSgenome.Mmusculus.UCSC.mm39, lib.loc="/home/viole/")
library(doParallel)
library(foreach)
library(tidyr)
library(tidyselect)
library(dplyr)

# define working directory
setwd("/proj/ancestry_medi_indiv/aligned")
# Load the list of filenames
filenames <- list.files(pattern = ".sorted.bam$")
# define output directory
output_dir <- file.path("/proj/ancestry_medi_indiv/merged")
# Define the maximum number of samples to run per subset
#max_samples_per_subset <- 10
# Set the number of cores to use
#numCores <- 16
# Divide filenames into subsets
#filename_subsets <- split(filenames, ceiling(seq_along(filenames) / max_samples_per_subset))
# Specify the input file and create the ROI
ROI <- read.delim(file.path(output_dir, "joined_sorted.bed"), sep = "", header = FALSE, stringsAsFactors = FALSE)
newname<-list("chr", "start", "end")
# Specify the input file, ROI, and other parameters
chr_all=unique(ROI$V1)
genome="BSgenome.Mmusculus.UCSC.mm39"
uniq=1e-3
extend=100
ws=300
shift=0
# Get the number of rows in the ROI data frame
n_rows <- nrow(ROI)
# Create a new data frame with the desired column names and numbers from 1 to n_rows
ROI <- data.frame(ROI$V1, ROI$V2, ROI$V3, 1:n_rows, stringsAsFactors = FALSE)
names(ROI)<-newname
newname<-list("chr", "start", "end", "name")
ROI <- subset(ROI, chr %in% chr_all)

# Create a function to process each subset of filenames
#process_subset <- function(subset_filenames) {
# Create a list to store the results for the subset
subset_results <- list()
# Process each file in the subset
for (file in filenames) {
    result <- MEDIPS.createROIset(file = file, BSgenome = genome, uniq = uniq, ROI = ROI, chr.select = chr_all, paired = TRUE)
    subset_results[[file]] <- result
}


# Initialize a parallel cluster
#cl <- makeCluster(numCores)

# Export required variables, packages, and functions to the parallel workers
#clusterExport(cl, varlist = c("genome", "uniq", "ROI", "chr_all", "MEDIPS.createROIset"))

# Load the required package in the parallel workers
#clusterEvalQ(cl, library(BSgenome.Mmusculus.UCSC.mm39))

# Process each subset in parallel
#subset_results <- parLapply(cl, filename_subsets, process_subset)

# Stop the parallel cluster
#stopCluster(cl)

# Merge the results from all subsets in the same order
r <- unlist(subset_results, recursive = FALSE)

# save the objects
save(r, file = file.path(output_dir, "r.rda"))

# remove all objects
rm(list = ls())

###run the necessary libraries and files again
load(file = file.path(output_dir, "r.rda"))
library(MEDIPS, lib.loc="/home/viole/")
library(BSgenome.Mmusculus.UCSC.mm39, lib.loc="/home/viole/")
library(tidyr)
library(tidyselect)
library(dplyr)

# This i don't know what is it
CS = MEDIPS.couplingVector(pattern = "CG", refObj = r[[1]])


methList <- list()  # Create an empty list to store the read assignments

for (i in seq_along(r)) {
  # Perform the MEDIPS.meth analysis
  meth_result <- MEDIPS.meth(r[[i]])

  # Remove the statistics from the result
  meth_result$statistics <- NULL

  # Store the read assignments in the list
  methList[[i]] <- meth_result
}

save(methList, file = file.path(output_dir, "methList.rda"))
load(file="methList.rda")

# count matrix
meth_countmatrix <- methList[[1]][, c("chr", "start", "stop")]

for (i in seq_along(methList)) {
  counts_column <- methList[[i]][[4]]
  column_name <- colnames(methList[[i]])[4]
  meth_countmatrix <- cbind(meth_countmatrix, counts_column)
  colnames(meth_countmatrix)[ncol(meth_countmatrix)] <- column_name
}

save(meth_countmatrix, file = file.path(output_dir, "meth_countmatrix.rda"))

# get a tab separated count matrix
write.table(meth_countmatrix, file = file.path(output_dir, "meth_countmatrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Outputting a message indicating the location of the output files
message <- paste("The output files 'meth_countmatrix.txt and .rda' are located at:", output_dir)
system(paste("echo", shQuote(message)))

# Calculate duration
end_time <- Sys.time()
duration <- end_time - start_time
message <- paste("Execution time:", format(duration, units = "auto"))
system(paste("echo", shQuote(message)))

EOF
