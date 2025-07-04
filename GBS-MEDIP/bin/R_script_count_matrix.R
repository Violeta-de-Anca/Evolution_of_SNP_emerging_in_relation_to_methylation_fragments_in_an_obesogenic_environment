#library(BiocManager)
#library(S4Vectors)
#library(BiocGenerics)
#library(IRanges)
#library(GenomeInfoDb)
#library(BSgenome)
library(MEDIPS)
library(BSgenome.Mmusculus.UCSC.mm39)
library(tidyr)
library(tidyselect)
library(dplyr)

# define working directory
setwd("/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/aligned")
# Load the list of filenames
filenames <- list.files(pattern = "Mouse.unique.MQ20.bam$")
# define output directory
output_dir <- file.path("/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/merged")
# Specify the input file and create the ROI
ROI <- read.delim(file.path(output_dir, "joined_sorted.bed"), sep = "", header = FALSE, stringsAsFactors = FALSE)
newname<-list("chr", "start", "end")
# Specify the input file, ROI, and other parameters
chr_all=unique(ROI$V1)
genome="BSgenome.Mmusculus.UCSC.mm39"
uniq=0

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
save(subset_results,file="/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/merged/subset_results.rda")
load("/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/merged/subset_results.rda")
r <- unlist(subset_results, recursive = FALSE)

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

meth_countmatrix <- methList[[1]][, c("chr", "start", "stop")]

for (i in seq_along(methList)) {
  counts_column <- methList[[i]][[4]]
  column_name <- colnames(methList[[i]])[4]
  meth_countmatrix <- cbind(meth_countmatrix, counts_column)
  colnames(meth_countmatrix)[ncol(meth_countmatrix)] <- column_name
}
save(meth_countmatrix, file = file.path(output_dir, "meth_countmatrix.rda"))
write.table(meth_countmatrix, file = file.path(output_dir, "meth_countmatrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


######################################################################################
###### trying with MACS3 windows in F1_2 individual ####
# define working directory
setwd("/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/aligned")
gbs_medip_dir=file.path("/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/aligned")
# Load the list of filenames
filename <- file.path(gbs_medip_dir, "C13F1_2_Mouse.unique.MQ20.bam")
# define output directory
output_dir <- file.path("/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/GBS_with_MACS3.win/")
# Specify the input file and create the ROI
ROI <- read.delim(file.path(output_dir, "C13F1_2.GBS.with.MACS3.win.intersect.bed"), sep = "", header = FALSE, stringsAsFactors = FALSE)
newname<-list("chr", "start", "end")
# Specify the input file, ROI, and other parameters
chr_all=unique(ROI$V1)
genome="BSgenome.Mmusculus.UCSC.mm39"
uniq=0

# Get the number of rows in the ROI data frame
n_rows <- nrow(ROI)
# Create a new data frame with the desired column names and numbers from 1 to n_rows
ROI <- data.frame(ROI$V1, ROI$V2, ROI$V3,ROI$V4, stringsAsFactors = FALSE)
newname<-list("chr", "start", "end", "name")
names(ROI)<-newname
#ROI <- subset(ROI, chr %in% chr_all)

# Process each file in the subset
result <- MEDIPS.createROIset(file = filename, BSgenome = genome, uniq = uniq, ROI = ROI, chr.select = chr_all, paired = TRUE)
meth_result <- MEDIPS.meth(result)
write.table(meth_result, file = file.path(output_dir, "F1_2_MEDIPS_countmatrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


######################################################################################
###### trying with MACS3 windows in F1_2 individual - isSecondaryAlignment = TRUE, simpleCigar = TRUE ####
# define working directory
setwd("/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/aligned")
gbs_medip_dir=file.path("/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/aligned")
# Load the list of filenames
filename <- file.path(gbs_medip_dir, "C13F1_2_Mouse.unique.MQ20.bam")
# define output directory
output_dir <- file.path("/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/verification.peaks.MACS3/intersect_win_per_indv/GBS_with_MACS3.win/")
# Specify the input file and create the ROI
ROI <- read.delim(file.path(output_dir, "C13F1_2.GBS.with.MACS3.win.intersect.bed"), sep = "", header = FALSE, stringsAsFactors = FALSE)
newname<-list("chr", "start", "end")
# Specify the input file, ROI, and other parameters
chr_all=unique(ROI$V1)
genome="BSgenome.Mmusculus.UCSC.mm39"
uniq=0

# Get the number of rows in the ROI data frame
n_rows <- nrow(ROI)
# Create a new data frame with the desired column names and numbers from 1 to n_rows
ROI <- data.frame(ROI$V1, ROI$V2, ROI$V3,ROI$V4, stringsAsFactors = FALSE)
newname<-list("chr", "start", "end", "name")
names(ROI)<-newname
#ROI <- subset(ROI, chr %in% chr_all)

# Process each file in the subset
result <- MEDIPS.createROIset(file = filename, BSgenome = genome, uniq = uniq, ROI = ROI, chr.select = chr_all, paired = TRUE,isSecondaryAlignment = TRUE, simpleCigar = TRUE)
meth_result <- MEDIPS.meth(result)
write.table(meth_result, file = file.path(output_dir, "F1_2_MEDIPS_seconalignandcigar_countmatrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
