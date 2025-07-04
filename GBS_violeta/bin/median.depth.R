library(ggplot2)
data <- read.table("/proj/naiss2023-23-55/GBS_violeta/recalibration/site.depth.combine.nonfiltered.ldepth", header = TRUE)
x <- as.data.frame(table(data$SUM_DEPTH))
lower <- 0.8 * median(data$SUM_DEPTH)
lower= median(data$SUM_DEPTH) - 1 * sd(data$SUM_DEPTH)
upper <- median(data$SUM_DEPTH) + 1 * sd(data$SUM_DEPTH)
xupper <- ceiling(upper/100) * 100
#jpeg(file = "/proj/naiss2023-23-55/GBS_violeta/quality_control/fastqc_output/snp_nofiltered/median.depth.snp.GBS.ICR.jpeg", width = 9000, height = 9000)
tiff("/proj/naiss2023-23-55/GBS_violeta/quality_control/fastqc_output/snp_nofiltered/median.depth.snp.GBS.ICR.tiff", height = 30, width = 20, units="cm", compression = "lzw", res = 300)
ggplot(x, aes(x = as.numeric(Var1) , y =as.numeric(Freq))) + geom_line() + xlab("Depth") + ylab("Freq") + xlim(0, xupper)+ylim(0,100000) + geom_vline(xintercept = lower, color = "red", linewidth = 1.3) + geom_vline(xintercept = upper, color = "red", linewidth = 1.3) + ggtitle("ICR GBS sperm SNP nonfiltered median depth + 2sd")
dev.off()                                                                             

