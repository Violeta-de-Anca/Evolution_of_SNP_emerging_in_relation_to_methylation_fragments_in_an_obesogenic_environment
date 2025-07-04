library(ggplot2)
data <- read.table("/proj/naiss2023-23-55/GBS_violeta/recalibration/site.quality.combine.nonfiltered.lqual", header = TRUE)
#jpeg(file ="/proj/naiss2023-23-55/GBS_violeta/quality_control/fastqc_output/snp_nofiltered/sitequality.nonfiltered.snp.GBS.ICR.jpeg", width = 9000, height = 9000)
tiff("/proj/naiss2023-23-55/GBS_violeta/quality_control/fastqc_output/snp_nofiltered/sitequality.nonfiltered.snp.GBS.ICR.tiff", height = 30, width = 20, units="cm", compression = "lzw", res = 300)
ggplot(subset(data, QUAL < 1000 & QUAL >= 0), aes(x = QUAL)) + geom_histogram(fill = "white",
    color = "black", bins = 50) + xlab("Quality value") + ylab("Count") + geom_vline(xintercept = 30,
    color = "red", size = 1.3) + ggtitle("ICR GBS SNP sperm nonfiltered")
dev.off()
