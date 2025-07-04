#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 6-00:00:00
#SBATCH -J ROI
#SBATCH --error /proj/ancestry_medi_indiv/roi.creation.err
#SBATCH --output /proj/ancestry_medi_indiv/roi.creation.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

#load modules
#module load bioinfo-tools
module load R_packages/4.2.1

#R code starts here
Rscript - <<EOF

#Capture start time
start_time <- Sys.time()

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
# Create a list to store the results for the subset
subset_results <- list()
# Process each file in the subset
for (file in filenames) {
  result <- MEDIPS.createROIset(file = file, BSgenome = genome, uniq = uniq, ROI = ROI, chr.select = chr_all, paired = TRUE)
  subset_results[[file]] <- result
}

save(subset_results,file="/proj/ancestry_medi_indiv/subset_results.rda")

EOF
