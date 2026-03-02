library(stringr)
# BiocManager::install("edgeR")
#install.packages("rgl")
#library(rgl)
#library(htmltools)
library(edgeR)
library("RColorBrewer")
library(dplyr)
#install.packages("PLNmodels")
#library(PLNmodels)
library(corrplot)
library(factoextra)
library(FactoMineR)
library(ade4)
library(gridExtra)
library(ggrepel)
library(viridis)
library(tidyr)
library(tidyselect)
library(dplyr)
library(ggplot2)
library("viridis") 
library(gplots)
library(qqman)
library(data.table)
library(patchwork)
library(scales)
library(tibble)
library(magrittr)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(grDevices)

# # load the data####
# stats.all.gen=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/meth_countmatrix.txt")
# groups.mus.musculus=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)
# rownames(stats.all.gen) <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# stats.all.gen = stats.all.gen[rowSums(stats.all.gen[,4:59]) >=1,]
# location.stats.all.gen <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# stats.all.gen$location <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# individuals.stats.all.gen=as.data.frame(colnames(stats.all.gen[,4:59]))
# colnames(individuals.stats.all.gen)=c("sample")
# groups.mus.musculus=left_join(individuals.stats.all.gen,groups.mus.musculus,by="sample")
# ##### mann whitney F1 #####
# F1.stats=stats.all.gen[,c(4:14,60)]
# row.names(F1.stats)=location.stats.all.gen
# F1.stats = F1.stats[rowSums(F1.stats[,1:11]) >=10,]
# groups.F1=groups.mus.musculus[1:11,]
# edgeR.F1 <- DGEList(counts=F1.stats[,1:11], genes=F1.stats[,12])
# # rownames(edgeR.F1$counts) <- rownames(edgeR.F1$genes) <- F1.stats[,12]
# edgeR.F1 <- calcNormFactors(edgeR.F1)
# efective.F1= as.data.frame(edgeR.F1$samples$lib.size*edgeR.F1$samples$norm.factors)
# edgeR.F1.MW=data.frame(mapply(`*`,F1.stats[,1:11],t(efective.F1)))
# p.values.F1.mann.whitney=apply(edgeR.F1.MW, 1, function(x){wilcox.test(x~groups.F1$sex)$p.value})
# p.values.F1.mann.whitney=as.data.frame(p.values.F1.mann.whitney)
# colnames(p.values.F1.mann.whitney)=c("p_value")
# p.values.F1.mann.whitney$location=F1.stats$location
# p.values.F1.mann.whitney$corrected_p_value.BH = p.adjust(((p.values.F1.mann.whitney$p_value)), method = "BH")
# p.values.F1.mann.whitney$corrected_p_value = p.adjust(((p.values.F1.mann.whitney$p_value)), method = "bonferroni")
# #p.values.F1.mann.whitney=p.values.F1.mann.whitney[p.values.F1.mann.whitney$p_value<0.05,]
# #write.csv(p.values.F1.mann.whitney,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/pvalues_testis_F1_MW_less_0.05.csv",quote = F,row.names = F)
# # getting the windows that are marked as significant
# MW.bonfe.sig.F1=p.values.F1.mann.whitney[p.values.F1.mann.whitney$p_value<0.01,]
# 
# ##### PCA+BCA approach #####
# #PCA to have the eigthteen vectors and to see how your data looks like
# pca_TMM_CPM.F1<-dudi.pca(edgeR.F1.MW,scannf=F,nf=5) #keep the first 5 axis
# fviz_pca_biplot(pca_TMM_CPM.F1,col.var='contrib',title = "PCA - Biplot, F1 generation")
# 
# #volcano plots, using the log fold change and the p-values F1####
# #get first the fold change, the fold change is calculated respective the control, so we need to change the sign of all the fold changes so is respective the obese group
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F1.rda")
# design.F1=model.matrix(~groups.F1$sex)
# edgeR.F1 <- estimateDisp(edgeR.F1, design.F1)
# x.F1=glmFit(edgeR.F1, design.F1)
# x.F1 = glmLRT(x.F1)
# y.F1=topTags(x.F1,n=nrow(edgeR.F1))
# logFC.F1=data.frame(logFC=y.F1$table$logFC,location=y.F1$table$location)
# z.F1=left_join(p.values.F1.mann.whitney,logFC.F1,by="location")
# z.F1$logFC <- -1 * z.F1$logFC
# z.F1$diffexpressed <- "NoDiff"
# z.F1$diffexpressed[z.F1$logFC > 0.6 & z.F1$p_value < 0.05] <- "Hypermethylated"
# z.F1$diffexpressed[z.F1$logFC < -0.6 & z.F1$p_value < 0.05] <- "Hypomethylated"
# z.F1$delabel <- NA
# z.F1$delabel[z.F1$diffexpressed != "NoDiff"] <- z.F1$location[z.F1$diffexpressed != "NoDiff"]
# save(z.F1,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F1.rda")
# # ggplot(data=z.F1, aes(x=logFC, y=-log10(p_value), col=diffexpressed, label=delabel)) + 
# #   geom_point()+
# #   geom_text_repel()+
# #   theme_minimal()+
# #   scale_color_manual(values=c("blue", "black", "red"))+
# #   geom_vline(xintercept=c(-0.6, 0.6), col="grey") +
# #   geom_hline(yintercept=-log10(0.05), col="grey")+
# #   geom_hline(yintercept=-log10(0.01), col="grey")+
# #   geom_hline(yintercept=-log10(0.001), col="grey")
# ggplot(data=z.F1, aes(x=logFC, y=-log10(p_value), col=diffexpressed)) +
#   geom_point(size=10)+
#   theme(legend.title = element_text(size=25),
#         panel.background=element_blank(),
#         plot.background=element_blank(),
#         axis.text.x=element_text(size=25),
#         axis.text.y =element_text(size=25),
#         legend.text = element_text(size=25),
#         axis.title = element_text(size=25))+
#   scale_color_manual(values=c("#990f4b","#6c88c4","black"))+
#   geom_vline(xintercept=c(-0.6, 0.6), col="grey") +
#   geom_hline(aes(yintercept=-log10(0.05),linetype="dotted"),linewidth=1.5,color="grey",show.legend=TRUE)+
#   geom_hline(aes(yintercept=-log10(0.001),linetype="dashed"),linewidth=1.5,color="grey",show.legend=TRUE)+
#   geom_hline(aes(yintercept=-log10(0.01),linetype="solid"),linewidth=1.5,color="grey",show.legend=TRUE)+
#   scale_linetype(name = "P-values", guide=guide_legend(override.aes = list(colour = "grey"))) +
#   ylab("-log10(p value)")+xlab("logFC methylation obese/control")+labs(col="Differentially\nmethylated")
# 
# 
# ##### mann whitney F2 #####
# F2.stats=stats.all.gen[,c(15:35,60)]
# row.names(F2.stats)=location.stats.all.gen
# F2.stats = F2.stats[rowSums(F2.stats[,1:21]) >=10,]
# groups.F2=groups.mus.musculus[12:32,]
# edgeR.F2 <- DGEList(counts=F2.stats[,1:21], genes=F2.stats[,22])
# # rownames(edgeR.F2$counts) <- rownames(edgeR.F2$genes) <- F2.stats[,22]
# edgeR.F2 <- calcNormFactors(edgeR.F2)
# efective.F2= as.data.frame(edgeR.F2$samples$lib.size*edgeR.F2$samples$norm.factors)
# edgeR.F2.MW=data.frame(mapply(`*`,F2.stats[,1:21],t(efective.F2)))
# p.values.F2.mann.whitney=apply(edgeR.F2.MW, 1, function(x){wilcox.test(x~groups.F2$sex)$p.value})
# p.values.F2.mann.whitney=as.data.frame(p.values.F2.mann.whitney)
# colnames(p.values.F2.mann.whitney)=c("p_value")
# p.values.F2.mann.whitney$location=F2.stats$location
# p.values.F2.mann.whitney$corrected_p_value.BH = p.adjust(((p.values.F2.mann.whitney$p_value)), method = "BH")
# p.values.F2.mann.whitney$corrected_p_value = p.adjust(((p.values.F2.mann.whitney$p_value)), method = "bonferroni")
# #p.values.F2.mann.whitney=p.values.F2.mann.whitney[p.values.F2.mann.whitney$p_value<0.05,]
# #write.csv(p.values.F2.mann.whitney,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/pvalues_testis_F2_MW_less_0.05.csv",quote = F,row.names = F)
# # getting the windows that are marked as significant
# MW.bonfe.sig.F2=p.values.F2.mann.whitney[p.values.F2.mann.whitney$corrected_p_value<0.01,]
# ##### PCA+BCA approach #####
# #PCA to have the eigthteen vectors and to see how your data looks like
# pca_TMM_CPM.F2<-dudi.pca(edgeR.F2.MW,scannf=F,nf=5) #keep the first 5 axis
# fviz_pca_biplot(pca_TMM_CPM.F2,col.var='contrib',title = "PCA - Biplot, F2 generation")
# 
# #volcano plots, using the log fold change and the p-values F2#####
# #get first the fold change
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F2.rda")
# design.F2=model.matrix(~groups.F2$sex)
# edgeR.F2 <- estimateDisp(edgeR.F2, design.F2)
# x.F2=glmFit(edgeR.F2, design.F2)
# x.F2 = glmLRT(x.F2)
# y.F2=topTags(x.F2,n=nrow(edgeR.F2))
# logFC.F2=data.frame(logFC=y.F2$table$logFC,location=y.F2$table$location)
# z.F2=left_join(p.values.F2.mann.whitney,logFC.F2,by="location")
# z.F2$logFC <- -1 * z.F2$logFC
# z.F2$diffexpressed <- "NoDiff"
# z.F2$diffexpressed[z.F2$logFC > 0.6 & z.F2$p_value < 0.05] <- "Hypermethylated"
# z.F2$diffexpressed[z.F2$logFC < -0.6 & z.F2$p_value < 0.05] <- "Hypomethylated"
# z.F2$delabel <- NA
# z.F2$delabel[z.F2$diffexpressed != "NoDiff"] <- z.F2$location[z.F2$diffexpressed != "NoDiff"]
# save(z.F2,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F2.rda")
# # ggplot(data=z.F2, aes(x=logFC, y=-log10(p_value), col=diffexpressed, label=delabel)) + 
# #   geom_point()+
# #   geom_text_repel()+
# #   theme_minimal()+
# #   scale_color_manual(values=c("blue", "black", "red"))+
# #   geom_vline(xintercept=c(-0.6, 0.6), col="grey") +
# #   geom_hline(yintercept=-log10(0.05), col="grey")+
# #   geom_hline(yintercept=-log10(0.01), col="grey")+
# #   geom_hline(yintercept=-log10(0.001), col="grey")
# # ggplot(data=z.F2, aes(x=logFC, y=-log10(p_value), col=diffexpressed)) +
# #   geom_point(size=2.5)+
# #   theme_minimal()+
# #   scale_color_manual(values=c("blue", "black", "red"))+
# #   geom_vline(xintercept=c(-0.6, 0.6), col="grey") +
# #   geom_hline(yintercept=-log10(0.05), col="yellow",size=1.5)+
# #   geom_hline(yintercept=-log10(0.01), col="orange")+
# #   geom_hline(yintercept=-log10(0.001), col="green")+
# #   ggtitle("F2 generation")+ylab("-log10(p value)")+
# #   labs(col = "Differentially\nexpressed")
# ggplot(data=z.F2, aes(x=logFC, y=-log10(p_value), col=diffexpressed)) +
#   geom_point(size=10)+
#   theme(legend.title = element_text(size=25),
#         panel.background=element_blank(),
#         plot.background=element_blank(),
#         axis.text.x=element_text(size=25),
#         axis.text.y =element_text(size=25),
#         legend.text = element_text(size=25),
#         axis.title = element_text(size=25))+
#   scale_color_manual(values=c("#990f4b","#6c88c4","black"))+
#   geom_vline(xintercept=c(-0.6, 0.6), col="grey") +
#   geom_hline(aes(yintercept=-log10(0.05),linetype="dotted"),linewidth=1.5,color="grey",show.legend=TRUE)+
#   geom_hline(aes(yintercept=-log10(0.001),linetype="dashed"),linewidth=1.5,color="grey",show.legend=TRUE)+
#   geom_hline(aes(yintercept=-log10(0.01),linetype="solid"),linewidth=1.5,color="grey",show.legend=TRUE)+
#   scale_linetype(name = "P-values", guide=guide_legend(override.aes = list(colour = "grey"))) +
#   ylab("-log10(p value)")+xlab("logFC methylation obese/control")+labs(col="Differentially\nmethylated")
# 
# 
# ##### mann whitney F3 #####
# F3.stats=stats.all.gen[,c(36:60)]
# row.names(F3.stats)=location.stats.all.gen
# F3.stats = F3.stats[rowSums(F3.stats[,1:24]) >=10,]
# groups.F3=groups.mus.musculus[33:56,]
# edgeR.F3 <- DGEList(counts=F3.stats[,1:24], genes=F3.stats[,25])
# # rownames(edgeR.F3$counts) <- rownames(edgeR.F3$genes) <- F3.stats[,25]
# edgeR.F3 <- calcNormFactors(edgeR.F3)
# efective.F3= as.data.frame(edgeR.F3$samples$lib.size*edgeR.F3$samples$norm.factors)
# edgeR.F3.MW=data.frame(mapply(`*`,F3.stats[,1:24],t(efective.F3)))
# p.values.F3.mann.whitney=apply(edgeR.F3.MW, 1, function(x){wilcox.test(x~groups.F3$sex)$p.value})
# p.values.F3.mann.whitney=as.data.frame(p.values.F3.mann.whitney)
# colnames(p.values.F3.mann.whitney)=c("p_value")
# p.values.F3.mann.whitney$location=F3.stats$location
# p.values.F3.mann.whitney$corrected_p_value.BH = p.adjust(((p.values.F3.mann.whitney$p_value)), method = "BH")
# p.values.F3.mann.whitney$corrected_p_value = p.adjust(((p.values.F3.mann.whitney$p_value)), method = "bonferroni")
# #p.values.F3.mann.whitney=p.values.F3.mann.whitney[p.values.F3.mann.whitney$p_value<0.05,]
# #write.csv(p.values.F3.mann.whitney,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/pvalues_testis_F3_MW_less_0.05.csv",quote = F,row.names = F)
# # getting the windows that are marked as significant
# MW.bonfe.sig.F3=p.values.F3.mann.whitney[p.values.F3.mann.whitney$corrected_p_value<0.01,]
# ##### PCA+BCA approach #####
# #PCA to have the eigthteen vectors and to see how your data looks like
# pca_TMM_CPM.F3<-dudi.pca(edgeR.F3.MW,scannf=F,nf=5) #keep the first 5 axis
# fviz_pca_biplot(pca_TMM_CPM.F3,col.var='contrib',title = "PCA - Biplot, F3 generation")
# 
# #volcano plots, using the log fold change and the p-values F3####
# #get first the fold change####
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F3.rda")
# design.F3=model.matrix(~groups.F3$sex)
# edgeR.F3 <- estimateDisp(edgeR.F3, design.F3)
# x.F3=glmFit(edgeR.F3, design.F3)
# x.F3 = glmLRT(x.F3)
# y.F3=topTags(x.F3,n=nrow(edgeR.F3))
# logFC.F3=data.frame(logFC=y.F3$table$logFC,location=y.F3$table$location)
# z.F3=left_join(p.values.F3.mann.whitney,logFC.F3,by="location")
# z.F3$logFC <- -1 * z.F3$logFC
# z.F3$diffexpressed <- "NoDiff"
# z.F3$diffexpressed[z.F3$logFC > 0.6 & z.F3$p_value < 0.05] <- "Hypermethylation"
# z.F3$diffexpressed[z.F3$logFC < -0.6 & z.F3$p_value < 0.05] <- "Hypomethylation"
# z.F3$delabel <- NA
# z.F3$delabel[z.F3$diffexpressed != "NoDiff"] <- z.F3$location[z.F3$diffexpressed != "NoDiff"]
# save(z.F3,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F3.rda")
# # ggplot(data=z.F3, aes(x=logFC, y=-log10(p_value), col=diffexpressed, label=delabel)) + 
# #   geom_point()+
# #   geom_text_repel()+
# #   theme_minimal()+
# #   scale_color_manual(values=c("blue", "black", "red"))+
# #   geom_vline(xintercept=c(-0.6, 0.6), col="grey") +
# #   geom_hline(yintercept=-log10(0.05), col="grey")+
# #   geom_hline(yintercept=-log10(0.01), col="grey")+
# #   geom_hline(yintercept=-log10(0.001), col="grey")
# # ggplot(data=z.F3, aes(x=logFC, y=-log10(p_value), col=diffexpressed)) +
# #   geom_point(size=2.5)+
# #   theme_minimal()+
# #   scale_color_manual(values=c("blue", "black", "red"))+
# #   geom_vline(xintercept=c(-0.6, 0.6), col="grey") +
# #   geom_hline(yintercept=-log10(0.05), col="yellow",size=1.5)+
# #   geom_hline(yintercept=-log10(0.01), col="orange")+
# #   geom_hline(yintercept=-log10(0.001), col="green")+
# #   ggtitle("F3 generation")+ylab("-log10(p value)")+
# #   labs(col = "Differentially\nexpressed")
# 
# ggplot(data=z.F3, aes(x=logFC, y=-log10(p_value), col=diffexpressed)) +
#   geom_point(size=10)+
#   theme(legend.title = element_text(size=25),
#         panel.background=element_blank(),
#         plot.background=element_blank(),
#         axis.text.x=element_text(size=25),
#         axis.text.y =element_text(size=25),
#         legend.text = element_text(size=25),
#         axis.title = element_text(size=25))+
#   scale_color_manual(values=c("#990f4b","#6c88c4","black"))+
#   geom_vline(xintercept=c(-0.6, 0.6), col="grey") +
#   geom_hline(aes(yintercept=-log10(0.05),linetype="dotted"),linewidth=1.5,color="grey",show.legend=TRUE)+
#   geom_hline(aes(yintercept=-log10(0.001),linetype="dashed"),linewidth=1.5,color="grey",show.legend=TRUE)+
#   geom_hline(aes(yintercept=-log10(0.01),linetype="solid"),linewidth=1.5,color="grey",show.legend=TRUE)+
#   scale_linetype(name = "P-values", guide=guide_legend(override.aes = list(colour = "grey"))) +
#   ylab("-log10(p value)")+xlab("logFC methylation obese/control")+labs(col="Differentially\nmethylated")
# 
# #get a volcano plot with all the DMRs and all the gens#####
# z.F1=z.F1[z.F1$p_value<0.05,]
# z.F1=z.F1[z.F1$logFC > 0.6|z.F1$logFC < -0.6,]
# z.F1$diffexpressed <- "NA"
# z.F1$diffexpressed[z.F1$logFC > 0.6 & z.F1$p_value < 0.05 & z.F1$p_value > 0.01] <- "Hypermethylation with p-value < 0.05"
# z.F1$diffexpressed[z.F1$logFC < -0.6 & z.F1$p_value < 0.05 & z.F1$p_value > 0.01] <- "Hypomethylation with p-value < 0.05"
# z.F1$diffexpressed[z.F1$logFC > 0.6 & z.F1$p_value <= 0.01 & z.F1$p_value > 0.001] <- "Hypermethylation with p-value < 0.01"
# z.F1$diffexpressed[z.F1$logFC < -0.6 & z.F1$p_value <= 0.01 & z.F1$p_value > 0.001] <- "Hypomethylation with p-value < 0.01"
# z.F1$diffexpressed[z.F1$logFC > 0.6 & z.F1$p_value < 0.001] <- "Hypermethylation with p-value < 0.001"
# z.F1$diffexpressed[z.F1$logFC < -0.6 & z.F1$p_value < 0.001] <- "Hypomethylation with p-value < 0.001"
# z.F1$Generation="F1"
# z.F2=z.F2[z.F2$p_value<0.05,]
# z.F2=z.F2[z.F2$logFC > 0.6|z.F2$logFC < -0.6,]
# z.F2$diffexpressed <- "NA"
# z.F2$diffexpressed[z.F2$logFC > 0.6 & z.F2$p_value < 0.05 & z.F2$p_value > 0.01] <- "Hypermethylation with p-value < 0.05"
# z.F2$diffexpressed[z.F2$logFC < -0.6 & z.F2$p_value < 0.05 & z.F2$p_value > 0.01] <- "Hypomethylation with p-value < 0.05"
# z.F2$diffexpressed[z.F2$logFC > 0.6 & z.F2$p_value <= 0.01 & z.F2$p_value > 0.001] <- "Hypermethylation with p-value < 0.01"
# z.F2$diffexpressed[z.F2$logFC < -0.6 & z.F2$p_value <= 0.01 & z.F2$p_value > 0.001] <- "Hypomethylation with p-value < 0.01"
# z.F2$diffexpressed[z.F2$logFC > 0.6 & z.F2$p_value < 0.001] <- "Hypermethylation with p-value < 0.001"
# z.F2$diffexpressed[z.F2$logFC < -0.6 & z.F2$p_value < 0.001] <- "Hypomethylation with p-value < 0.001"
# z.F2$Generation="F2"
# z.F3=z.F3[z.F3$p_value<0.05,]
# z.F3=z.F3[z.F3$logFC > 0.6|z.F3$logFC < -0.6,]
# z.F3$diffexpressed <- "NA"
# z.F3$diffexpressed[z.F3$logFC > 0.6 & z.F3$p_value < 0.05 & z.F3$p_value > 0.01] <- "Hypermethylation with p-value < 0.05"
# z.F3$diffexpressed[z.F3$logFC < -0.6 & z.F3$p_value < 0.05 & z.F3$p_value > 0.01] <- "Hypomethylation with p-value < 0.05"
# z.F3$diffexpressed[z.F3$logFC > 0.6 & z.F3$p_value <= 0.01 & z.F3$p_value > 0.001] <- "Hypermethylation with p-value < 0.01"
# z.F3$diffexpressed[z.F3$logFC < -0.6 & z.F3$p_value <= 0.01 & z.F3$p_value > 0.001] <- "Hypomethylation with p-value < 0.01"
# z.F3$diffexpressed[z.F3$logFC > 0.6 & z.F3$p_value < 0.001] <- "Hypermethylation with p-value < 0.001"
# z.F3$diffexpressed[z.F3$logFC < -0.6 & z.F3$p_value < 0.001] <- "Hypomethylation with p-value < 0.001"
# z.F3$Generation="F3"
# all.z=rbind(z.F1,z.F2,z.F3)
# 
# ggplot(data=all.z, aes(x=logFC, y=-log10(p_value), col=diffexpressed, shape=Generation)) +
#   geom_point(size=10,position = "jitter")+
#   theme(legend.title = element_text(size=25),
#         panel.background=element_blank(),
#         plot.background=element_blank(),
#         axis.text.x=element_text(size=25),
#         axis.text.y =element_text(size=25),
#         legend.text = element_text(size=25),
#         axis.title = element_text(size=25))+
#   scale_color_manual(values=c("#980909", "#d82909", "#e87f6b","#405268","#4e8992","black"))+
#   scale_shape_manual(values=c(19,15,17))+
#   ylab("-log10(p value)")+xlab("logFC methylation obese/control")+labs(col="Differentially\nmethylated")
# 
# ggplot(data=all.z, aes(x=logFC, y=-log10(p_value), col=Generation, shape=diffexpressed)) +
#   geom_point(size=10)+
#   theme(legend.title = element_text(size=25),
#         panel.background=element_blank(),
#         plot.background=element_blank(),
#         axis.text.x=element_text(size=25),
#         axis.text.y =element_text(size=25),
#         legend.text = element_text(size=25),
#         axis.title = element_text(size=25))+
#   scale_color_manual(values=c("#980909", "#d82909", "#e87f6b"))+
#   scale_shape_manual(values=c(19,15,17))+
#   ylab("-log10(p value)")+xlab("logFC methylation obese/control")+labs(col="Differentially\nmethylated")
# 
# # heatmap of the DMRs normalized counts, all individuals ####
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/DMRs.allgens.normalized.counts.rda")
# names(DMRs.allgens.normalized.counts) <- gsub("F1", "F0", names(DMRs.allgens.normalized.counts))
# names(DMRs.allgens.normalized.counts) <- gsub("F2", "F1", names(DMRs.allgens.normalized.counts))
# names(DMRs.allgens.normalized.counts) <- gsub("F3", "F2", names(DMRs.allgens.normalized.counts))
# groups.mus.musculus=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)
#  #### Heatmap F0 ####
# groups.f0 = groups.mus.musculus[groups.mus.musculus$patient == "F0", ]
# group.F0.obese = groups.f0[groups.f0$sex == "Obese", ]
# group.F0.control = groups.f0[groups.f0$sex == "Control", ]
# #columnas DMRs, filas individuos
# DMRs.all.scaled=scale(DMRs.allgens.normalized.counts, center = F, scale = T)
#   #we cannot do all generations at the same time as is too big
# F0.o=group.F0.obese$sample
# F0.c=group.F0.control$sample
# DMRs.f0.obese.scaled <- DMRs.all.scaled[, colnames(DMRs.all.scaled) %in% F0.o]
# DMRs.f0.control.scaled <- DMRs.all.scaled[, colnames(DMRs.all.scaled) %in% F0.c]
# DMRs.f0.scaled=cbind(DMRs.f0.control.scaled,DMRs.f0.obese.scaled)
# tiff("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/heatmap.F0.DMR.GBS.MeDIP.normalizedcounts.scaled.tiff", width = 100000, height = 10000, units = "px", res = 1000)
# par(mar=c(1,1,1,1))
# heatmap.2(DMRs.f0.scaled,dendrogram='none',col=viridis(56) , key=F ,margins=c(5,9), xlab = "DMRs",ylab = "Individuals",Rowv=FALSE, Colv=FALSE,trace='none',labCol = colnames(DMRs.f0.scaled))
# dev.off()
# 
#   #### Heatmap F1 ####
# groups.F1 = groups.mus.musculus[groups.mus.musculus$patient == "F1", ]
# group.F1.obese = groups.F1[groups.F1$sex == "Obese", ]
# group.F1.control = groups.F1[groups.F1$sex == "Control", ]
# #columnas DMRs, filas individuos
# DMRs.all.scaled=scale(DMRs.allgens.normalized.counts, center = F, scale = T)
# #we cannot do all generations at the same time as is too big
# F1.o=group.F1.obese$sample
# F1.c=group.F1.control$sample
# DMRs.F1.obese.scaled <- DMRs.all.scaled[, colnames(DMRs.all.scaled) %in% F1.o]
# DMRs.F1.control.scaled <- DMRs.all.scaled[, colnames(DMRs.all.scaled) %in% F1.c]
# DMRs.F1.scaled=cbind(DMRs.F1.control.scaled,DMRs.F1.obese.scaled)
# tiff("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/heatmap.F1.DMR.GBS.MeDIP.normalizedcounts.scaled.tiff", width = 100000, height = 10000, units = "px", res = 1000)
# heatmap.2(DMRs.F1.scaled,dendrogram='none',col=viridis(56) , key=F ,margins=c(5,9), xlab = "Individuals",ylab = "DMRs",Rowv=FALSE, Colv=FALSE,trace='none',labCol = colnames(DMRs.F1.scaled))
# dev.off()
# 
# #### Heatmap F2 ####
# groups.F2 = groups.mus.musculus[groups.mus.musculus$patient == "F2", ]
# group.F2.obese = groups.F2[groups.F2$sex == "Obese", ]
# group.F2.control = groups.F2[groups.F2$sex == "Control", ]
# #columnas DMRs, filas individuos
# DMRs.all.scaled=scale(DMRs.allgens.normalized.counts, center = F, scale = T)
# #we cannot do all generations at the same time as is too big
# F2.o=group.F2.obese$sample
# F2.c=group.F2.control$sample
# DMRs.F2.obese.scaled <- DMRs.all.scaled[, colnames(DMRs.all.scaled) %in% F2.o]
# DMRs.F2.control.scaled <- DMRs.all.scaled[, colnames(DMRs.all.scaled) %in% F2.c]
# DMRs.F2.scaled=cbind(DMRs.F2.control.scaled,DMRs.F2.obese.scaled)
# tiff("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/heatmap.F2.DMR.GBS.MeDIP.normalizedcounts.scaled.tiff", width = 100000, height = 10000, units = "px", res = 1000)
# heatmap.2(DMRs.F2.scaled,dendrogram='none',col=viridis(56) , key=F ,margins=c(5,9), xlab = "Individuals",ylab = "DMRs",Rowv=FALSE, Colv=FALSE,trace='none',labCol = colnames(DMRs.F2.scaled))
# dev.off()
# 
# #### Heatmap of fold change, across gen ####
# #those were heatmaps with normalized counts for all DMRs, now we do heatmap of 
# #the fold change, so we can see how the DMRs change across generations
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F3.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F2.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F1.rda")
# z.F1$gen = "F0" 
# z.F2$gen = "F1"
# z.F3$gen = "F2"
# f0.hiper=z.F1[ z.F1$diffexpressed=="Hypermethylated",]
# f1.hipe=z.F2[ z.F2$diffexpressed=="UP",]
# f2.hipe=z.F3[ z.F3$diffexpressed=="UP",]
# 
# #NO SCALING, heatmap F0 hipermethylated across gens ####
# fc.f1=z.F2[z.F2$location %in% f0.hiper$location,]
# fc.f2=z.F3[z.F3$location %in% f0.hiper$location,]
# all.hiper=left_join(f0.hiper,fc.f1,by="location")
# all.hiper=data.frame(location=all.hiper$location,logFC.F0=all.hiper$logFC.x,logFC.F1=all.hiper$logFC.y)
# all.hiper=left_join(all.hiper,fc.f2,by="location")
# all.hiper=data.frame(location=all.hiper$location,logFC.F0=all.hiper$logFC.F0,logFC.F1=all.hiper$logFC.F1,logFC.F2=all.hiper$logFC)
# all.hiper=all.hiper[complete.cases(all.hiper),]
# row.names(all.hiper)=all.hiper$location
# all.hiper=all.hiper[,-1]
# all.hiper=as.matrix(all.hiper)
# tiff("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/plots/heatmap.foldchange.F0.across.gens.onlyheatmap.tiff",width = 15000, height = 5000, res = 600, units = "px")
# heatmap.2(all.hiper,dendrogram='none', lhei = c(1,8),lwid = c(0.5,4),col = bluered(100), key=T,key.title ="",key.xlab = "log(Fold Change)",key.ylab = "" ,cexRow =2,cexCol =1.5,notecex=2,margins=c(3,1), xlab = "Generations",ylab = "DMRs",Rowv=FALSE, Colv=FALSE,trace='none',scale ="none",labCol = colnames(all.hiper))
# dev.off()
# 
# #NO SCALING, heatmap F1 hipermethylated across gens ####
# fc.f0=z.F1[z.F1$location %in% f1.hipe$location,]
# fc.f2=z.F3[z.F3$location %in% f1.hipe$location,]
# all.hiper=left_join(f1.hipe,fc.f0,by="location")
# all.hiper=data.frame(location=all.hiper$location,logFC.F1=all.hiper$logFC.x,logFC.F0=all.hiper$logFC.y)
# all.hiper=left_join(all.hiper,fc.f2,by="location")
# all.hiper=data.frame(location=all.hiper$location,logFC.F0=all.hiper$logFC.F0,logFC.F1=all.hiper$logFC.F1,logFC.F2=all.hiper$logFC)
# all.hiper=all.hiper[complete.cases(all.hiper),]
# row.names(all.hiper)=all.hiper$location
# all.hiper=all.hiper[,-1]
# all.hiper=as.matrix(all.hiper)
# tiff("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/plots/heatmap.foldchange.F1.across.gens.tiff",width = 15000, height = 5000, res = 600, units = "px")
# heatmap.2(all.hiper,lhei = c(1,4),lwid = c(0.5,4),dendrogram='none',cex=2,col = bluered(100), key=F,key.title ="",key.xlab = "log(Fold Change)",key.ylab = "" ,cexRow =1.5,cexCol =2,notecex=1,margins=c(9,24), xlab = "",ylab = "",Rowv=FALSE, Colv=FALSE,trace='none',scale ="none",labCol = colnames(all.hiper))
# dev.off()
# 
# #NO SCALING, heatmap F2 hipermethylated across gens ####
# fc.f0=z.F1[z.F1$location %in% f1.hipe$location,]
# fc.f1=z.F2[z.F2$location %in% f0.hiper$location,]
# all.hiper=left_join(f2.hipe,fc.f0,by="location")
# all.hiper=data.frame(location=all.hiper$location,logFC.F2=all.hiper$logFC.x,logFC.F0=all.hiper$logFC.y)
# all.hiper=left_join(all.hiper,fc.f1,by="location")
# all.hiper=data.frame(location=all.hiper$location,logFC.F0=all.hiper$logFC.F0,logFC.F2=all.hiper$logFC.F2,logFC.F1=all.hiper$logFC)
# all.hiper=all.hiper[complete.cases(all.hiper),]
# row.names(all.hiper)=all.hiper$location
# all.hiper=all.hiper[,-1]
# all.hiper=as.matrix(all.hiper)
# tiff("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/plots/heatmap.foldchange.F1.across.gens.tiff",width = 15000, height = 5000, res = 600, units = "px")
# heatmap.2(all.hiper,lhei = c(1,4),lwid = c(0.5,4),dendrogram='none',cex=2,col = bluered(100), key=F,key.title ="",key.xlab = "log(Fold Change)",key.ylab = "" ,cexRow =1.5,cexCol =2,notecex=1,margins=c(9,24), xlab = "",ylab = "",Rowv=FALSE, Colv=FALSE,trace='none',scale ="none",labCol = colnames(all.hiper))
# dev.off()
# 
# f0.hipo=z.F1[ z.F1$diffexpressed=="Hypomethylated",]
# f1.hipo=z.F2[ z.F2$diffexpressed=="DOWN",]
# f2.hipo=z.F3[ z.F3$diffexpressed=="DOWN",]
# 
# #NO SCALING, heatmap F0 hipomethylated across gens ####
# fc.f1=z.F2[z.F2$location %in% f0.hipo$location,]
# fc.f2=z.F3[z.F3$location %in% f0.hipo$location,]
# all.hipo=left_join(f0.hipo,fc.f1,by="location")
# all.hipo=data.frame(location=all.hipo$location,logFC.F0=all.hipo$logFC.x,logFC.F1=all.hipo$logFC.y)
# all.hipo=left_join(all.hipo,fc.f2,by="location")
# all.hipo=data.frame(location=all.hipo$location,logFC.F0=all.hipo$logFC.F0,logFC.F1=all.hipo$logFC.F1,logFC.F2=all.hipo$logFC)
# all.hipo=all.hipo[complete.cases(all.hipo),]
# row.names(all.hipo)=all.hipo$location
# all.hipo=all.hipo[,-1]
# all.hipo=as.matrix(all.hipo)
# tiff("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/plots/heatmap.foldchange.F0DMRhipo.across.gens.tiff",width = 15000, height = 5000, res = 600, units = "px")
# heatmap.2(all.hipo,dendrogram='none', lhei = c(1,8),lwid = c(0.5,4),col = bluered(100), key=F,key.title ="",key.xlab = "log(Fold Change)",key.ylab = "" ,cexRow =1,cexCol =1.5,notecex=2,margins=c(7,50), xlab = "",ylab = "",Rowv=FALSE, Colv=FALSE,trace='none',scale ="none",labCol = colnames(all.hipo))
# dev.off()
# 
# #NO SCALING, heatmap f1 hipomethylated across gens ####
# fc.f0=z.F1[z.F1$location %in% f1.hipo$location,]
# fc.f2=z.F3[z.F3$location %in% f1.hipo$location,]
# all.hipo=left_join(f1.hipo,fc.f0,by="location")
# all.hipo=data.frame(location=all.hipo$location,logFC.f1=all.hipo$logFC.x,logFC.F0=all.hipo$logFC.y)
# all.hipo=left_join(all.hipo,fc.f2,by="location")
# all.hipo=data.frame(location=all.hipo$location,logFC.F0=all.hipo$logFC.F0,logFC.F1=all.hipo$logFC.f1,logFC.F2=all.hipo$logFC)
# all.hipo=all.hipo[complete.cases(all.hipo),]
# row.names(all.hipo)=all.hipo$location
# all.hipo=all.hipo[,-1]
# all.hipo=as.matrix(all.hipo)
# tiff("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/plots/heatmap.foldchange.f1DMRhipo.across.gens.tiff",width = 15000, height = 5000, res = 600, units = "px")
# heatmap.2(all.hipo,dendrogram='none',col = bluered(100), key=T,key.title ="",key.xlab = "log(Fold Change)",key.ylab = "" ,cexRow =2,cexCol =1.5,notecex=2,margins=c(7,50), xlab = "",ylab = "",Rowv=FALSE, Colv=FALSE,trace='none',scale ="none",labCol = colnames(all.hipo))
# dev.off()
# 
# #NO SCALING, heatmap F2 hipomethylated across gens ####
# fc.f0=z.F1[z.F1$location %in% f2.hipo$location,]
# fc.f1=z.F2[z.F2$location %in% f2.hipo$location,]
# all.hipo=left_join(f2.hipo,fc.f0,by="location")
# all.hipo=data.frame(location=all.hipo$location,logFC.F2=all.hipo$logFC.x,logFC.F0=all.hipo$logFC.y)
# all.hipo=left_join(all.hipo,fc.f1,by="location")
# all.hipo=data.frame(location=all.hipo$location,logFC.F0=all.hipo$logFC.F0,logFC.F1=all.hipo$logFC,logFC.F2=all.hipo$logFC.F2)
# all.hipo=all.hipo[complete.cases(all.hipo),]
# row.names(all.hipo)=all.hipo$location
# all.hipo=all.hipo[,-1]
# all.hipo=as.matrix(all.hipo)
# tiff("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/plots/heatmap.foldchange.F2DMRhipo.across.gens.tiff",width = 15000, height = 5000, res = 600, units = "px")
# heatmap.2(all.hipo,dendrogram='none',col = bluered(100),key=T,key.title ="",key.xlab = "log(Fold Change)",key.ylab = "" ,cexRow =2,cexCol =1.5,notecex=2,margins=c(7,50), xlab = "",ylab = "",Rowv=FALSE, Colv=FALSE,trace='none',scale ="none",labCol = colnames(all.hipo))
# dev.off()
# lhei = c(1,8),lwid = c(0.5,4)
# 
# #heatmaps of these windows, but with normalized counts ####
#  #heatmap F0 ####
# f0.hiper
# f0.hipo
# all.F0.FC.win=c(f0.hiper$location,f0.hipo$location)
# DMRs.allgens.normalized.counts$location=row.names(DMRs.allgens.normalized.counts)
# DMRs.FC.F0=DMRs.allgens.normalized.counts[DMRs.allgens.normalized.counts$location %in% all.F0.FC.win,]
# F0.o=group.F0.obese$sample
# F0.c=group.F0.control$sample
# DMRs.f0.obese <- DMRs.FC.F0[, colnames(DMRs.FC.F0) %in% F0.o]
# DMRs.f0.control <- DMRs.FC.F0[, colnames(DMRs.FC.F0) %in% F0.c]
# DMRs.f0=cbind(DMRs.f0.control,DMRs.f0.obese)
# DMRs.f0=as.matrix(DMRs.f0)
# par(mar=c(1,1,1,1))
# DMRs.f0.1=log(DMRs.f0)
# DMRs.f0.1[DMRs.f0.1==-Inf] <- 0
# tiff("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/plots/heatmap.F0.DMR.FC.GBS.MeDIP.log.normalizedcounts.tiff", width = 100000, height = 10000, units = "px", res = 1000)
# heatmap.2(DMRs.f0.1,dendrogram='none', col = brewer.pal(n = 9, name = "Reds"),key.title ="",key=T ,key.ylab = "" ,key.xlab = "log(normalized counts)",margins=c(5,13),Rowv=FALSE, Colv=FALSE,trace='none',labCol = colnames(DMRs.f0),scale ="none")
# dev.off()
# brewer.pal(n = 3, name = "PuRd") Reds
# bluered(100)
# 
# #heatmap F1 ####
# f1.hipo
# f1.hipe
# all.F1.FC.win=c(f1.hipe$location,f1.hipo$location)
# DMRs.allgens.normalized.counts$location=row.names(DMRs.allgens.normalized.counts)
# DMRs.FC.F1=DMRs.allgens.normalized.counts[DMRs.allgens.normalized.counts$location %in% all.F1.FC.win,]
# F1.o=group.F1.obese$sample
# F1.c=group.F1.control$sample
# DMRs.F1.obese <- DMRs.FC.F1[, colnames(DMRs.FC.F1) %in% F1.o]
# DMRs.F1.control <- DMRs.FC.F1[, colnames(DMRs.FC.F1) %in% F1.c]
# DMRs.F1=cbind(DMRs.F1.control,DMRs.F1.obese)
# DMRs.F1=as.matrix(DMRs.F1)
# DMRs.F1.1=log(DMRs.F1)
# DMRs.F1.1[DMRs.F1.1==-Inf] <- 0
# tiff("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/plots/heatmap.F1.DMR.FC.GBS.MeDIP.log.normalizedcounts.tiff", width = 100000, height = 10000, units = "px", res = 1000)
# heatmap.2(DMRs.F1.1,dendrogram='none', col = brewer.pal(n = 9, name = "Reds"),key.title ="",key=T ,key.ylab = "" ,key.xlab = "log(normalized counts)",margins=c(5,13),Rowv=FALSE, Colv=FALSE,trace='none',labCol = colnames(DMRs.F1),scale ="none")
# dev.off()
# brewer.pal(n = 3, name = "PuRd") Reds
# bluered(100)
# 
# #heatmap F2 ####
# f2.hipo
# f2.hipe
# all.F2.FC.win=c(f2.hipe$location,f2.hipo$location)
# DMRs.allgens.normalized.counts$location=row.names(DMRs.allgens.normalized.counts)
# DMRs.FC.F2=DMRs.allgens.normalized.counts[DMRs.allgens.normalized.counts$location %in% all.F2.FC.win,]
# F2.o=group.F2.obese$sample
# F2.c=group.F2.control$sample
# DMRs.F2.obese <- DMRs.FC.F2[, colnames(DMRs.FC.F2) %in% F2.o]
# DMRs.F2.control <- DMRs.FC.F2[, colnames(DMRs.FC.F2) %in% F2.c]
# DMRs.F2=cbind(DMRs.F2.control,DMRs.F2.obese)
# DMRs.F2=as.matrix(DMRs.F2)
# DMRs.F2.1=log(DMRs.F2)
# DMRs.F2.1[DMRs.F2.1==-Inf] <- 0
# tiff("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/plots/heatmap.F2.DMR.FC.GBS.MeDIP.log.normalizedcounts.tiff", width = 100000, height = 10000, units = "px", res = 1000)
# heatmap.2(DMRs.F2.1,dendrogram='none', col = brewer.pal(n = 9, name = "Reds"),key.title ="",key=T ,key.ylab = "" ,key.xlab = "log(normalized counts)",margins=c(5,13),Rowv=FALSE, Colv=FALSE,trace='none',labCol = colnames(DMRs.F2),scale ="none")
# dev.off()
# brewer.pal(n = 3, name = "PuRd") Reds
# bluered(100)
# 
# #overlapping of DMRs between generations####
# #F2vsF1####
# overlapp.F2vsF1=inner_join(z.F2,z.F1,by="location")
# overlapp.F2vsF1$both[overlapp.F2vsF1$diffexpressed.x=="UP" & overlapp.F2vsF1$diffexpressed.y=="UP"] <- T
# overlapp.F2vsF1$both[overlapp.F2vsF1$diffexpressed.x=="DOWN" & overlapp.F2vsF1$diffexpressed.y=="DOWN"] <- T
# overlapp.F2vsF1=overlapp.F2vsF1[complete.cases(overlapp.F2vsF1),]
# write.csv(overlapp.F2vsF1,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/overlap.DMRs.F2vsF1.csv",quote = F,row.names = F)
# #F2vsF3####
# overlapp.F2vsF3=inner_join(z.F2,z.F3,by="location")
# overlapp.F2vsF3$both[overlapp.F2vsF3$diffexpressed.x=="UP" & overlapp.F2vsF3$diffexpressed.y=="UP"] <- T
# overlapp.F2vsF3$both[overlapp.F2vsF3$diffexpressed.x=="DOWN" & overlapp.F2vsF3$diffexpressed.y=="DOWN"] <- T
# overlapp.F2vsF3=overlapp.F2vsF3[complete.cases(overlapp.F2vsF3),]
# write.csv(overlapp.F2vsF3,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/overlap.DMRs.F2vsF3.csv",quote = F,row.names = F)
# 
# # let's do a vennagram plot to see if any of the methylated windows are also with SNPs####
# #we know that is 91 F1, 64 F2 and 142 F3, between each are 2 windows common:
# #install.packages("ggvenn")
# library(ggvenn)
# #install.packages("wesanderson")
# library(wesanderson)
# z.F1=z.F1[z.F1$diffexpressed=="UP"|z.F1$diffexpressed=="DOWN",]
# z.F2=z.F2[z.F2$diffexpressed=="UP"|z.F2$diffexpressed=="DOWN",]
# z.F3=z.F3[z.F3$diffexpressed=="UP"|z.F3$diffexpressed=="DOWN",]
# cololeh=wes_palette("GrandBudapest1", 3, type = c("discrete"))
# ggvenn(list(F1=z.F1$location,F2=z.F2$location,F3=z.F3$location),show_percentage=F,fill_color =cololeh,text_size = 10,digits=0)
# 
# 
# ########################### overlapp of SNPs and DMRs#####
# z.F3=z.F3[complete.cases(z.F3),]
# z.F1=z.F1[complete.cases(z.F1),]
# z.F2=z.F2[complete.cases(z.F2),]
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS_templeton/F1_to_F3/R_objects/comparision.c.F2.vs.F3.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS_templeton/F1_to_F3/R_objects/comparision.o.F2.vs.F3.vcf.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS_templeton/F1_to_F3/R_objects/comparision.F3.vcf.1.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS_templeton/F1_to_F3/R_objects/comparision.F2.vcf.1.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS_templeton/F1_to_F3/R_objects/comparision.F2controlvsF1.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS_templeton/F1_to_F3/R_objects/comparision.F2obesevsF1.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F3.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F2.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F1.rda")
# ###########################################################################################
# #comparision f2 methylation with SNPs inheritance from F1 and between novel SNPs in F2obese####
# comparision.F2obesevsF1.vcf$coordinates=rownames(comparision.F2obesevsF1.vcf)
# list.com=matrix(unlist(strsplit(comparision.F2obesevsF1.vcf$coordinates,split="_")),ncol=2,byrow=T)
# comparision.F2obesevsF1.vcf$chromosome=list.com[,1]
# comparision.F2obesevsF1.vcf$location=list.com[,2]
# 
# list.com=matrix(unlist(strsplit(z.F2$location,split=":")),ncol=2,byrow=T)
# z.F2$chromosome=list.com[,1]
# z.F2$coordinate=list.com[,2]
# list.com=matrix(unlist(strsplit(z.F2$coordinate,split="-")),ncol=2,byrow=T)
# z.F2$start=list.com[,1]
# z.F2$end=list.com[,2]
# row.names(z.F2)=z.F2$location
# 
# comparision.F2obesevsF1.vcf=subset(comparision.F2obesevsF1.vcf,chromosome%in%unique(z.F2$chromosome))
# 
# #now check if the SNP is in between the range of the windows####
# #between(x, left, right)
# overlapp.F2.DMRs.SNPs.F2o.vs.F1=as.data.frame(matrix(unlist(lapply(comparision.F2obesevsF1.vcf$location, function(x){between(x,lower = z.F2$start,upper = z.F2$end)})),ncol=64,byrow=T)) #THAT 64 IS BECAUSE IS THE NUMBER OF DMRS IF IS NOT THE NUMBER OF DMRS IT DOES NOT WORK!!!
# row.names(overlapp.F2.DMRs.SNPs.F2o.vs.F1)=comparision.F2obesevsF1.vcf$coordinates
# colnames(overlapp.F2.DMRs.SNPs.F2o.vs.F1)=z.F2$location
# logical_vector <- apply(overlapp.F2.DMRs.SNPs.F2o.vs.F1 == "TRUE", 1, any)
# overlapp.F2.DMRs.SNPs.F2o.vs.F1 <- subset(overlapp.F2.DMRs.SNPs.F2o.vs.F1, logical_vector)
# # save(overlapp.F2.DMRs.SNPs.F2o.vs.F1,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/overlapp.F2.DMRs.SNPs.F2o.vs.F1.rda")
# # write.csv(overlapp.F2.DMRs.SNPs.F2o.vs.F1,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/overlapp.F2.DMRs.SNPs.F2o.vs.F1.csv",quote = F,row.names = T)
# 
# ############################################################################################
# #comparision f2 methylation with SNPs inheritance from F1 and between novel SNPs in F2control####
# comparision.F2controlvsF1.vcf$coordinates=rownames(comparision.F2controlvsF1.vcf)
# list.com=matrix(unlist(strsplit(comparision.F2controlvsF1.vcf$coordinates,split="_")),ncol=2,byrow=T)
# comparision.F2controlvsF1.vcf$chromosome=list.com[,1]
# comparision.F2controlvsF1.vcf$location=list.com[,2]
# 
# list.com=matrix(unlist(strsplit(z.F2$location,split=":")),ncol=2,byrow=T)
# z.F2$chromosome=list.com[,1]
# z.F2$coordinate=list.com[,2]
# list.com=matrix(unlist(strsplit(z.F2$coordinate,split="-")),ncol=2,byrow=T)
# z.F2$start=list.com[,1]
# z.F2$end=list.com[,2]
# row.names(z.F2)=z.F2$location
# z.F2=z.F2[z.F2$diffexpressed=="UP"|z.F2$diffexpressed=="DOWN",]
# comparision.F2controlvsF1.vcf=subset(comparision.F2controlvsF1.vcf,chromosome%in%unique(z.F2$chromosome))
# 
# #now check if the SNP is in between the range of the windows####
# #between(x, left, right)
# overlapp.F2.DMRs.SNPs.F2c.vs.F1=as.data.frame(matrix(unlist(lapply(comparision.F2controlvsF1.vcf$location, function(x){between(x,lower = z.F2$start,upper = z.F2$end)})),ncol=64,byrow=T))
# row.names(overlapp.F2.DMRs.SNPs.F2c.vs.F1)=comparision.F2controlvsF1.vcf$coordinates
# colnames(overlapp.F2.DMRs.SNPs.F2c.vs.F1)=z.F2$location
# logical_vector <- apply(overlapp.F2.DMRs.SNPs.F2c.vs.F1 == "TRUE", 1, any)
# overlapp.F2.DMRs.SNPs.F2c.vs.F1 <- subset(overlapp.F2.DMRs.SNPs.F2c.vs.F1, logical_vector)
# save(overlapp.F2.DMRs.SNPs.F2c.vs.F1,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/overlapp.F2.DMRs.SNPs.F2c.vs.F1.rda")
# write.csv(overlapp.F2.DMRs.SNPs.F2c.vs.F1,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/overlapp.F2.DMRs.SNPs.F2c.vs.F1.csv",quote = F,row.names = T)
# 
# ############################################################################################
# #comparision f2 methylation with SNPs inheritance from F1 and between novel SNPs comparing F2 obese vs F2 control####
# comparision.F2.vcf.1$coordinates=rownames(comparision.F2.vcf.1)
# list.com=matrix(unlist(strsplit(comparision.F2.vcf.1$coordinates,split="_")),ncol=2,byrow=T)
# comparision.F2.vcf.1$chromosome=list.com[,1]
# comparision.F2.vcf.1$location=list.com[,2]
# 
# list.com=matrix(unlist(strsplit(z.F2$location,split=":")),ncol=2,byrow=T)
# z.F2$chromosome=list.com[,1]
# z.F2$coordinate=list.com[,2]
# list.com=matrix(unlist(strsplit(z.F2$coordinate,split="-")),ncol=2,byrow=T)
# z.F2$start=list.com[,1]
# z.F2$end=list.com[,2]
# row.names(z.F2)=z.F2$location
# z.F2=z.F2[z.F2$diffexpressed=="UP"|z.F2$diffexpressed=="DOWN",]
# comparision.F2.vcf.1=subset(comparision.F2.vcf.1,chromosome%in%unique(z.F2$chromosome))
# 
# #now check if the SNP is in between the range of the windows####
# #between(x, left, right)
# overlapp.F2.DMRs.SNPs.F2ovsF2c=as.data.frame(matrix(unlist(lapply(comparision.F2.vcf.1$location, function(x){between(x,lower = z.F2$start,upper = z.F2$end)})),ncol=64,byrow=T))
# row.names(overlapp.F2.DMRs.SNPs.F2ovsF2c)=comparision.F2.vcf.1$coordinates
# colnames(overlapp.F2.DMRs.SNPs.F2ovsF2c)=z.F2$location
# logical_vector <- apply(overlapp.F2.DMRs.SNPs.F2ovsF2c == "TRUE", 1, any)
# overlapp.F2.DMRs.SNPs.F2ovsF2c <- subset(overlapp.F2.DMRs.SNPs.F2ovsF2c, logical_vector)
# save(overlapp.F2.DMRs.SNPs.F2ovsF2c,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/overlapp.F2.DMRs.SNPs.F2ovsF2c.rda")
# write.csv(overlapp.F2.DMRs.SNPs.F2ovsF2c,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/overlapp.F2.DMRs.SNPs.F2ovsF2c.csv",quote = F,row.names = T)
# 
# ############################################################################################
# #comparision F3 methylation with SNPs inheritance from F1 and between novel SNPs comparing F3 obese vs F3 control####
# comparision.F3.vcf.1$coordinates=rownames(comparision.F3.vcf.1)
# list.com=matrix(unlist(strsplit(comparision.F3.vcf.1$coordinates,split="_")),ncol=2,byrow=T)
# comparision.F3.vcf.1$chromosome=list.com[,1]
# comparision.F3.vcf.1$location=list.com[,2]
# 
# list.com=matrix(unlist(strsplit(z.F3$location,split=":")),ncol=2,byrow=T)
# z.F3$chromosome=list.com[,1]
# z.F3$coordinate=list.com[,2]
# list.com=matrix(unlist(strsplit(z.F3$coordinate,split="-")),ncol=2,byrow=T)
# z.F3$start=list.com[,1]
# z.F3$end=list.com[,2]
# row.names(z.F3)=z.F3$location
# z.F3=z.F3[z.F3$diffexpressed=="UP"|z.F3$diffexpressed=="DOWN",]
# comparision.F3.vcf.1=subset(comparision.F3.vcf.1,chromosome%in%unique(z.F3$chromosome))
# 
# #now check if the SNP is in between the range of the windows####
# #between(x, left, right)
# overlapp.F3.DMRs.SNPs.F3ovsF3c=as.data.frame(matrix(unlist(lapply(comparision.F3.vcf.1$location, function(x){between(x,lower = z.F3$start,upper = z.F3$end)})),ncol=142,byrow=T))
# row.names(overlapp.F3.DMRs.SNPs.F3ovsF3c)=comparision.F3.vcf.1$coordinates
# colnames(overlapp.F3.DMRs.SNPs.F3ovsF3c)=z.F3$location
# logical_vector <- apply(overlapp.F3.DMRs.SNPs.F3ovsF3c == "TRUE", 1, any)
# overlapp.F3.DMRs.SNPs.F3ovsF3c <- subset(overlapp.F3.DMRs.SNPs.F3ovsF3c, logical_vector)
# # save(overlapp.F3.DMRs.SNPs.F3ovsF3c,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/overlapp.F3.DMRs.SNPs.F3ovsF3c.rda")
# # write.csv(overlapp.F3.DMRs.SNPs.F3ovsF3c,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/overlapp.F3.DMRs.SNPs.F3ovsF3c.csv",quote = F,row.names = T)
# 
# ############################################################################################
# #comparision F3 methylation with SNPs inheritance from F1 and between novel SNPs in F3control####
# comparision.c.F2.vs.F3.vcf$coordinates=rownames(comparision.c.F2.vs.F3.vcf)
# list.com=matrix(unlist(strsplit(comparision.c.F2.vs.F3.vcf$coordinates,split="_")),ncol=2,byrow=T)
# comparision.c.F2.vs.F3.vcf$chromosome=list.com[,1]
# comparision.c.F2.vs.F3.vcf$location=list.com[,2]
# 
# list.com=matrix(unlist(strsplit(z.F3$location,split=":")),ncol=2,byrow=T)
# z.F3$chromosome=list.com[,1]
# z.F3$coordinate=list.com[,2]
# list.com=matrix(unlist(strsplit(z.F3$coordinate,split="-")),ncol=2,byrow=T)
# z.F3$start=list.com[,1]
# z.F3$end=list.com[,2]
# row.names(z.F3)=z.F3$location
# z.F3=z.F3[z.F3$diffexpressed=="UP"|z.F3$diffexpressed=="DOWN",]
# 
# comparision.c.F2.vs.F3.vcf=subset(comparision.c.F2.vs.F3.vcf,chromosome%in%unique(z.F3$chromosome))
# 
# #now check if the SNP is in between the range of the windows####
# #between(x, left, right)
# overlapp.F3.DMRs.SNPs.F3c.vs.F2=as.data.frame(matrix(unlist(lapply(comparision.c.F2.vs.F3.vcf$location, function(x){between(x,lower = z.F3$start,upper = z.F3$end)})),ncol=142,byrow=T))
# row.names(overlapp.F3.DMRs.SNPs.F3c.vs.F2)=comparision.c.F2.vs.F3.vcf$coordinates
# colnames(overlapp.F3.DMRs.SNPs.F3c.vs.F2)=z.F3$location
# logical_vector <- apply(overlapp.F3.DMRs.SNPs.F3c.vs.F2 == "TRUE", 1, any)
# overlapp.F3.DMRs.SNPs.F3c.vs.F2 <- subset(overlapp.F3.DMRs.SNPs.F3c.vs.F2, logical_vector)
# save(overlapp.F3.DMRs.SNPs.F3c.vs.F2,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/overlapp.F3.DMRs.SNPs.F3c.vs.F2.rda")
# write.csv(overlapp.F3.DMRs.SNPs.F3c.vs.F2,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/overlapp.F3.DMRs.SNPs.F3c.vs.F2.csv",quote = F,row.names = T)
# 
# ############################################################################################
# #comparision F3 methylation with SNPs inheritance from F1 and between novel SNPs in F3obese####
# comparision.o.F2.vs.F3.vcf$coordinates=rownames(comparision.o.F2.vs.F3.vcf)
# list.com=matrix(unlist(strsplit(comparision.o.F2.vs.F3.vcf$coordinates,split="_")),ncol=2,byrow=T)
# comparision.o.F2.vs.F3.vcf$chromosome=list.com[,1]
# comparision.o.F2.vs.F3.vcf$location=list.com[,2]
# 
# list.com=matrix(unlist(strsplit(z.F3$location,split=":")),ncol=2,byrow=T)
# z.F3$chromosome=list.com[,1]
# z.F3$coordinate=list.com[,2]
# list.com=matrix(unlist(strsplit(z.F3$coordinate,split="-")),ncol=2,byrow=T)
# z.F3$start=list.com[,1]
# z.F3$end=list.com[,2]
# row.names(z.F3)=z.F3$location
# z.F3=z.F3[z.F3$diffexpressed=="UP"|z.F3$diffexpressed=="DOWN",]
# 
# comparision.o.F2.vs.F3.vcf=subset(comparision.o.F2.vs.F3.vcf,chromosome%in%unique(z.F3$chromosome))
# 
# #now check if the SNP is in between the range of the windows####
# #between(x, left, right)
# overlapp.F3.DMRs.SNPs.F3o.vs.F2=as.data.frame(matrix(unlist(lapply(comparision.o.F2.vs.F3.vcf$location, function(x){between(x,lower = z.F3$start,upper = z.F3$end)})),ncol=142,byrow=T))
# row.names(overlapp.F3.DMRs.SNPs.F3o.vs.F2)=comparision.o.F2.vs.F3.vcf$coordinates
# colnames(overlapp.F3.DMRs.SNPs.F3o.vs.F2)=z.F3$location
# logical_vector <- apply(overlapp.F3.DMRs.SNPs.F3o.vs.F2 == "TRUE", 1, any)
# overlapp.F3.DMRs.SNPs.F3o.vs.F2 <- subset(overlapp.F3.DMRs.SNPs.F3o.vs.F2, logical_vector)
# save(overlapp.F3.DMRs.SNPs.F3o.vs.F2,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/overlapp.F3.DMRs.SNPs.F3o.vs.F2.rda")
# write.csv(overlapp.F3.DMRs.SNPs.F3o.vs.F2,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/overlapp.F3.DMRs.SNPs.F3o.vs.F2.csv",quote = F,row.names = T)
# 
# #let's do a PCA to see if the DMRs separate in the same way as the SNP PCA####
# #first get all the DMR, raw counts for F1
# z.F1=z.F1%>%subset(z.F1$diffexpressed=="UP"|z.F1$diffexpressed=="DOWN")
# z.F2=z.F2%>%subset(z.F2$diffexpressed=="UP"|z.F2$diffexpressed=="DOWN")
# z.F3=z.F3%>%subset(z.F3$diffexpressed=="UP"|z.F3$diffexpressed=="DOWN")
# stats.all.gen=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/meth_countmatrix.txt")
# rownames(stats.all.gen) <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# stats.all.gen = stats.all.gen[rowSums(stats.all.gen[,4:59]) >=1,]
# location.stats.all.gen <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# stats.all.gen$location <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# stats.F1=left_join(z.F1,stats.all.gen,by="location")
# stats.F2=left_join(z.F2,stats.all.gen,by="location")
# stats.F3=left_join(z.F3,stats.all.gen,by="location")
# stats.F1=stats.F1[,c(2,11:66)]
# stats.F2=stats.F2[,c(2,15:70)]
# stats.F3=stats.F3[,c(2,15:70)]
# DMRs.allgens.raw.counts=bind_rows(stats.F1,stats.F2)
# DMRs.allgens.raw.counts=bind_rows(DMRs.allgens.raw.counts,stats.F3)
# DMRs.allgens.raw.counts=unique(DMRs.allgens.raw.counts)
# #save(DMRs.allgens.raw.counts,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/DMRs.allgens.raw.counts.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/DMRs.allgens.raw.counts.rda")
# row.names(DMRs.allgens.raw.counts)=DMRs.allgens.raw.counts$location
# DMRs.allgens.raw.counts.1=DMRs.allgens.raw.counts[,2:57]
# DMRs.allgens.raw.counts=DMRs.allgens.raw.counts[,2:57]
# #DMRs.allgens.raw.counts=t(DMRs.allgens.raw.counts)
# #DMRs.allgens.raw.counts=matrix(DMRs.allgens.raw.counts,ncol=288,nrow = 56)
# DMRs.allgens.raw.counts <- DGEList(counts=DMRs.allgens.raw.counts[,], genes=row.names(DMRs.allgens.raw.counts))
# rownames(DMRs.allgens.raw.counts$counts) <- rownames(DMRs.allgens.raw.counts$genes) <- row.names(DMRs.allgens.raw.counts)
# DMRs.allgens.raw.counts <- calcNormFactors(DMRs.allgens.raw.counts)
# efective.DMRs= as.data.frame(DMRs.allgens.raw.counts$samples$lib.size*DMRs.allgens.raw.counts$samples$norm.factors)
# DMRs.allgens.raw.counts=data.frame(mapply(`*`,DMRs.allgens.raw.counts.1[,],t(efective.DMRs)))
# DMRs.allgens.normalized.counts=DMRs.allgens.raw.counts
# #save(DMRs.allgens.normalized.counts,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/DMRs.allgens.normalized.counts.rda")
# row.names(DMRs.allgens.normalized.counts)=row.names(DMRs.allgens.raw.counts.1)
# pca <- prcomp(t(DMRs.allgens.raw.counts), scale. = TRUE)
# biplot(pca,var.axes = F,cex=c(1,0.1))
# 
# pca <- prcomp(t(DMRs.allgens.normalized.counts), scale. = TRUE)
# groups=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)
# groups.DMRs=as.data.frame(row.names(pca$x))
# colnames(groups.DMRs)=c("sample")
# groups=left_join(groups.DMRs,groups,by="sample")
# pca=as.data.frame(pca$x)
# ggplot(pca,aes(x=PC1,y=PC2,color=interaction(groups$patient,groups$sex)))+geom_point(size=15)+
#   geom_hline(yintercept = 0,linetype="dotted")+labs(y= "PC2", x = "PC1")+
#   geom_vline(xintercept = 0, linetype="dotted")+#ggtitle("Principal component analysis \n with DMRs")+
#   theme(legend.title = element_blank(),axis.text.x=element_text(size=50),axis.text.y =element_text(size=50),legend.text = element_text(size=50),axis.title = element_text(size=50))
# 
# #we see that there are a couple of outliers, F1.8, F2.21 and F1.2, let's exclude them, plotting the 1 and 2 component#####
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/DMRs.allgens.normalized.counts.rda")
# groups=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)
# groups$sample <- gsub("-", ".", groups$sample)
# pca <- prcomp(t(DMRs.allgens.normalized.counts[,c(-4,-10,-23)]), scale. = TRUE,rank. = 3)
# pca.ind=as.data.frame(pca$x)
# groups.DMRs=as.data.frame(row.names(pca$x))
# colnames(groups.DMRs)=c("sample")
# groups=left_join(groups.DMRs,groups,by="sample")
# ggplot(pca.ind,aes(x=PC1,y=PC2,color=interaction(groups$patient,groups$sex)))+geom_point(size=15)+
#   geom_hline(yintercept = 0,linetype="dotted")+labs(y= "PC2", x = "PC1")+
#   geom_vline(xintercept = 0, linetype="dotted")+scale_color_manual(values=c("#d2a9b3", "#e74269", "#9b1746","#c2dd9b","#85bc37","#204c26"))+
#   theme(legend.title = element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=50),axis.text.y =element_text(size=50),legend.text = element_text(size=50),axis.title = element_text(size=50))
# 
# #let's do the 3D PCA ####
# intera.group=interaction(groups$patient,groups$sex)
# #plot3D::points3D(x=pca.ind$PC1, y=pca.ind$PC2, z=pca.ind$PC3, xlab ="PC1", ylab ="PC2", zlab ="PC3",col=as.integer(intera.group))
# plot3d(x=pca.ind$PC1, y=pca.ind$PC2, z=pca.ind$PC3,size = 8,type = "p", xlab ="PC1", ylab ="PC2", zlab ="PC3",col=as.integer(intera.group))
# text3d(pca.ind[,1:3],
#        texts=c(intera.group), 
#        cex= 0.7, pos=3)
# legend3d("topright", legend = paste(unique(intera.group)), pch = 16, col = unique(as.integer(intera.group)), cex=1, inset=c(0.02))
# rgl.snapshot('C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/plots/3dplot_GBS_MEDIP_ICR.png', fmt = 'png')
# 
# #with specific colors ####
# #want we need to do is to do a vector of the same lenght of the pca matrix 
# groups$color="NA"
# groups$color[groups$patient=="F0" & groups$sex=="Control"] <- "#d2a9b3"
# groups$color[groups$patient=="F0" & groups$sex=="Obese"] <- "#c2dd9b"
# groups$color[groups$patient=="F1" & groups$sex=="Control"] <- "#e74269"
# groups$color[groups$patient=="F1" & groups$sex=="Obese"] <- "#85bc37"
# groups$color[groups$patient=="F2" & groups$sex=="Control"] <- "#9b1746"
# groups$color[groups$patient=="F2" & groups$sex=="Obese"] <- "#204c26"
# intera.group=interaction(groups$patient,groups$sex)
# #plot3D::points3D(x=pca.ind$PC1, y=pca.ind$PC2, z=pca.ind$PC3, xlab ="PC1", ylab ="PC2", zlab ="PC3",col=as.integer(intera.group))
# plot3d(x=pca.ind$PC1, y=pca.ind$PC2, z=pca.ind$PC3,size = 8,type = "p", xlab ="PC1", ylab ="PC2", zlab ="PC3",col = groups$color)
# text3d(pca.ind[,1:3],
#        texts=c(intera.group), 
#        cex= 0.7, pos=3)
# legend3d("topright", legend = paste(unique(intera.group)), pch = 16, col =c("#d2a9b3","#c2dd9b", "#e74269","#85bc37", "#9b1746","#204c26"), cex=1, inset=c(0.02))
# rgl.snapshot('C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/plots/3dplot_GBS_MEDIP_ICR.png', fmt = 'png')
# 
# #let's do the 3D PCA but by adding gens
# pca.ind.f0=pca.ind[row.names(pca.ind) %in% groups.f0$sample,]
# #only F0
# plot3d(x=pca.ind.f0$PC1, y=pca.ind.f0$PC2, z=pca.ind.f0$PC3,size = 8,type = "p", xlab ="PC1", ylab ="PC2", zlab ="PC3",col=as.integer(interaction(groups.f0$patient,groups.f0$sex)))
# text3d(pca.ind.f0[,1:3],
#        texts=c(interaction(groups.f0$patient,groups.f0$sex)), 
#        cex= 0.7, pos=3)
# legend3d("topright", legend = paste(unique(interaction(groups.f0$patient,groups.f0$sex))), pch = 16, col = unique(as.integer(interaction(groups.f0$patient,groups.f0$sex))), cex=1, inset=c(0.02))
# 
# #let's do the 3D PCA but by adding gens
# group.f0.f1=rbind(groups.f0,groups.F1)
# pca.ind.f0=pca.ind[row.names(pca.ind) %in% group.f0.f1$sample,]
# #only F0
# plot3d(x=pca.ind.f0$PC1, y=pca.ind.f0$PC2, z=pca.ind.f0$PC3,size = 8,type = "p", xlab ="PC1", ylab ="PC2", zlab ="PC3",col=as.integer(interaction(group.f0.f1$patient,group.f0.f1$sex)))
# text3d(pca.ind.f0[,1:3],
#        texts=c(interaction(group.f0.f1$patient,group.f0.f1$sex)), 
#        cex= 0.7, pos=3)
# legend3d("topright", legend = paste(unique(interaction(group.f0.f1$patient,group.f0.f1$sex))), pch = 16, col = unique(as.integer(interaction(group.f0.f1$patient,group.f0.f1$sex))), cex=1, inset=c(0.02))
# 
# summary(pca)
# biplot(pca,var.axes = F,cex=c(1,0.1))
# pca$x
# 
# #we see that there are a couple of outliers, F1.8, F2.21 and F1.2, let's exclude them, plotting the 2 and 3 component#####
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/DMRs.allgens.normalized.counts.rda")
# groups=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)
# groups$sample <- gsub("-", ".", groups$sample)
# pca <- prcomp(t(DMRs.allgens.normalized.counts[,c(-4,-10,-23)]), scale. = TRUE,rank. = 3)
# pca.ind=as.data.frame(pca$x)
# groups.DMRs=as.data.frame(row.names(pca$x))
# colnames(groups.DMRs)=c("sample")
# groups=left_join(groups.DMRs,groups,by="sample")
# ggplot(pca.ind,aes(x=PC2,y=PC3,color=interaction(groups$patient,groups$sex)))+geom_point(size=15)+
#   geom_hline(yintercept = 0,linetype="dotted")+labs(y= "PC3", x = "PC2")+
#   geom_vline(xintercept = 0, linetype="dotted")+
#   theme(legend.title = element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=50),axis.text.y =element_text(size=50),legend.text = element_text(size=50),axis.title = element_text(size=50))
# 
# summary(pca)
# biplot(pca,var.axes = F,cex=c(1,0.1))
# pca$x
# 
# ##### PCA+BCA approach #####
# #PCA to have the eigthteen vectors and to see how your data looks like
# pca_TMM_CPM.DMrs<-dudi.pca(t(DMRs.allgens.raw.counts),scannf=F,nf=5) #keep the first 5 axis
# fviz_pca_biplot(pca_TMM_CPM.DMrs,col.var='contrib',title = "PCA - Biplot, DMR in all generations")
# 
# # let's annotate the coordinates####
# # BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
# # BiocManager::install("GenomicRanges")
# library(GenomicRanges)
# library(BSgenome.Mmusculus.UCSC.mm10)
# # BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# # BiocManager::install("org.Mm.eg.db")
# library(org.Mm.eg.db)
# # BiocManager::install("annotate")
# library(annotate)
# # if (!require("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# # install.packages("stringr")
# library(stringr)
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/DMRs.allgens.raw.counts.rda")
# row.names(DMRs.allgens.normalized.counts)=DMRs.allgens.raw.counts[,1]
# #first we need to put from the location column a dataframe with chromosome, start and end position
# list.com=matrix(unlist(strsplit(DMRs.allgens.raw.counts$location,split=":")),ncol=2,byrow=T)
# DMRs.allgens.coordinates=data.frame(Chromosome=list.com[,1])
# list.com=matrix(unlist(strsplit(list.com[,2],split="-")),ncol=2,byrow=T)
# DMRs.allgens.coordinates$Start=list.com[,1]
# DMRs.allgens.coordinates$End=list.com[,2]
# 
# # The function "makeGRangesFromDataFrame" from the library 
# # GenomicRanges makes an object GRanges called "intervals" from "myIntervals"
# DMRs.allgens.intervals = GenomicRanges::makeGRangesFromDataFrame(DMRs.allgens.coordinates)
# txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
# # extract the list of all gene_Ids from the txdb object
# genes = genes(txdb)
# 
# 
# # Make the annotating function. It will annotate the intervals with gene_Ids
# annotateIntervals <-  function(intervals, txdb)
#   {
#     stopifnot(is(intervals, "GRanges"), is(txdb, "TxDb"))
#     anno = genes(txdb)
#     olaps = findOverlaps(intervals, anno)
#     mcols(olaps)$gene_id = genes$gene_id[subjectHits(olaps)]
#     intervals_factor = factor(queryHits(olaps), levels=seq_len(queryLength(olaps)))
#     intervals$gene_id = splitAsList(mcols(olaps)$gene_id, intervals_factor)
#     intervals
#   }        
# 
# # Use the "annotateIntervals" funtion in order to annotate 
# #the intervals with gene_Ids and produce "myAnnotation" data.frame 
# myAnnotation.DMRs <- as.data.frame(annotateIntervals(DMRs.allgens.intervals, txdb))
# # Make an empty data.frame for append all the annotations, 
# #(we can call it "the master")
# myDf_master <- data.frame()
# 
# 
# # Now we want Hugo gene names in our annotations! 
# #So, for each annotated interval get hugo gene names...
# for (i in 1:length(myAnnotation.DMRs$gene_id)) {
#   # if the gene list is not empty...
#   if(length(c(na.omit(myAnnotation.DMRs$gene_id[i])[[1]])) != 0) {
#     # annotate the interval and copy into a myDf data.frame
#     myDf <- data.frame(myAnnotation.DMRs$seqnames[i], myAnnotation.DMRs$start[i], 
#                        myAnnotation.DMRs$end[i], toString(unname(getSYMBOL(c(na.omit(myAnnotation.DMRs$gene_id[i])[[1]]), data='org.Mm.eg.db'))))
#     # append tge myDF annotations with rbind into the myDf_master
#     myDf_master <- rbind(myDf_master, myDf)
#   }
# }
# 
# # Make a new header for the master dataframe: 
# # it will be the header of our output file also...
# myDf_header <- c("Chromosome", "Start", "End", "Genes")
# names(myDf_master) <- myDf_header
# myDf_master$location <- paste0(myDf_master$Chromosome, ":", myDf_master$Start, "-", myDf_master$End)
# 
# # write myDf_master content in a comma separated output file.
# write.csv(myDf_master, file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/annotated.DMRs.allgen.csv", row.names = FALSE,quote = F)
# 
# #now let's put together all the info for each of the generations
# DMRs.annotated.F1=left_join(z.F1,myDf_master,by="location")
# DMRs.annotated.F2=left_join(z.F2,myDf_master,by="location")
# DMRs.annotated.F3=left_join(z.F3,myDf_master,by="location")
# DMRs.annotated.F1=DMRs.annotated.F1[,c(-7,-8,-9,-10)]
# DMRs.annotated.F2=DMRs.annotated.F2[,c(-7,-8,-9,-10)]
# DMRs.annotated.F3=DMRs.annotated.F3[,c(-7,-8,-9,-10)]
# 
# #put the NCBI name + a little info con the genes####
# relation.table=fread(file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/b543f503-9da9-415e-9034-ad8b13c44873",col.names = T)
# relation.table$Genes=relation.table$`Source Id`
# DMRs.annotated.F1=left_join(DMRs.annotated.F1,relation.table,by="Genes")
# DMRs.annotated.F1=DMRs.annotated.F1[,c(-8,-10)]
# DMRs.annotated.F2=left_join(DMRs.annotated.F2,relation.table,by="Genes")
# DMRs.annotated.F2=DMRs.annotated.F2[,c(-8,-10)]
# DMRs.annotated.F3=left_join(DMRs.annotated.F3,relation.table,by="Genes")
# DMRs.annotated.F3=DMRs.annotated.F3[,c(-8,-10)]
# save(DMRs.annotated.F1,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/DMRs.annotated.F1.rda")
# write.csv(DMRs.annotated.F1, file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/DMRs.annotated.F1.csv", row.names = FALSE,quote = F)
# save(DMRs.annotated.F2,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/DMRs.annotated.F2.rda")
# write.csv(DMRs.annotated.F2, file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/DMRs.annotated.F2.csv", row.names = FALSE,quote = F)
# save(DMRs.annotated.F3,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/DMRs.annotated.F3.rda")
# write.csv(DMRs.annotated.F3, file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/DMRs.annotated.F3.csv", row.names = FALSE,quote = F)
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/DMRs.annotated.F1.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/DMRs.annotated.F2.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/DMRs.annotated.F3.rda")
# 
# ###for the pathway enrichment####
# #first create their Text file format for expression dataset
# #first column with gene names, second with description, the rest with the counts, normalized!
# stats.all.gen=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/meth_countmatrix.txt")
# groups.mus.musculus=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)
# rownames(stats.all.gen) <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# stats.all.gen = stats.all.gen[rowSums(stats.all.gen[,4:59]) >=1,]
# location.stats.all.gen <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# stats.all.gen$location <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# DMRs.annotated.F1=left_join(DMRs.annotated.F1,stats.all.gen,by="location")
# DMRs.counts.F1=DMRs.annotated.F1[DMRs.annotated.F1$diffexpressed == "UP" | DMRs.annotated.F1$diffexpressed == "DOWN",]
# DMRs.counts.F1=DMRs.counts.F1[complete.cases(DMRs.counts.F1),]
# DMRs.counts.F1=DMRs.counts.F1[,c(7,9,13:68)]
# names(DMRs.counts.F1)[names(DMRs.counts.F1) == 'Genes'] <- 'Name'
# names(DMRs.counts.F1)[names(DMRs.counts.F1) == 'Gene Description'] <- 'Description'
# DMRs.counts.F1=DMRs.counts.F1[,1:13]
# write.table(DMRs.counts.F1,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/DMRs.counts.F1.txt",row.names = F, quote = F,sep = "\t")
# 
# DMRs.annotated.F2=left_join(DMRs.annotated.F2,stats.all.gen,by="location")
# DMRs.counts.F2=DMRs.annotated.F2[DMRs.annotated.F2$diffexpressed == "UP" | DMRs.annotated.F2$diffexpressed == "DOWN",]
# DMRs.counts.F2=DMRs.counts.F2[complete.cases(DMRs.counts.F2),]
# DMRs.counts.F2=DMRs.counts.F2[,c(7,9,13:68)]
# names(DMRs.counts.F2)[names(DMRs.counts.F2) == 'Genes'] <- 'Name'
# names(DMRs.counts.F2)[names(DMRs.counts.F2) == 'Gene Description'] <- 'Description'
# DMRs.counts.F2=DMRs.counts.F2[,c(1:2,14:34)]
# write.table(DMRs.counts.F2,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/DMRs.counts.F2.txt",row.names = F, quote = F,sep = "\t")
# 
# DMRs.annotated.F3=left_join(DMRs.annotated.F3,stats.all.gen,by="location")
# DMRs.counts.F3=DMRs.annotated.F3[DMRs.annotated.F3$diffexpressed == "UP" | DMRs.annotated.F3$diffexpressed == "DOWN",]
# DMRs.counts.F3=DMRs.counts.F3[complete.cases(DMRs.counts.F3),]
# DMRs.counts.F3=DMRs.counts.F3[,c(7,9,13:68)]
# names(DMRs.counts.F3)[names(DMRs.counts.F3) == 'Genes'] <- 'Name'
# names(DMRs.counts.F3)[names(DMRs.counts.F3) == 'Gene Description'] <- 'Description'
# DMRs.counts.F3=DMRs.counts.F3[,c(1:2,35:58)]
# write.table(DMRs.counts.F3,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/DMRs.counts.F3.txt",row.names = F, quote = F,sep = "\t")
# 
# #now let's create the Categorical class file format
# groups=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)
# groups.DMRs=as.data.frame(colnames(DMRs.counts.F1[,c(-1,-2)]))
# colnames(groups.DMRs)=c("sample")
# groups=left_join(groups.DMRs,groups,by="sample")
# groups.F1=groups[1:11,]
# groups.F1=groups.F1$sex
# groups.F2=groups[12:32,]
# groups.F2=groups.F2$sex
# groups.F3=groups[33:56,]
# groups.F3=groups.F3$sex
# write.csv(groups.F1,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/groups.f1.csv",quote = F,row.names = F)
# write.csv(groups.F2,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/groups.f2.csv",quote = F,row.names = F)
# write.csv(groups.F3,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/groups.f3.csv",quote = F,row.names = F)
# 
# #so let's do again the files for gsea, the thing is now we are going to annotate everything that we have, not only what is significant
# # let's annotate ALL the coordinates####
# library(GenomicRanges)
# library(BSgenome.Mmusculus.UCSC.mm10)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# library(org.Mm.eg.db)
# library(annotate)
# library(stringr)
# library(data.table)
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F1.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F2.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F3.rda")
# stats.all.gen=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/meth_countmatrix.txt")
# rownames(stats.all.gen) <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# stats.all.gen = stats.all.gen[rowSums(stats.all.gen[,4:59]) >=1,]
# location.stats.all.gen <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# stats.all.gen$location <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# edgeR.F1.MW$location=z.F1$location
# edgeR.F2.MW$location=z.F2$location
# edgeR.F3.MW$location=z.F3$location
# stats.F1=left_join(z.F1,edgeR.F1.MW,by="location")
# stats.F2=left_join(z.F2,edgeR.F2.MW,by="location")
# stats.F3=left_join(z.F3,edgeR.F3.MW,by="location")
# stats.F1=stats.F1[,c(2,8:18)]
# stats.F2=stats.F2[,c(2,8:28)]
# stats.F3=stats.F3[,c(2,8:31)]
# # save(stats.F1,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/normalized.counts.F1.rda")
# # save(stats.F2,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/normalized.counts.F2.rda")
# # save(stats.F3,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/normalized.counts.F3.rda")
# # F1 annotation all locations ####
# row.names(stats.F1)=stats.F1[,1]
# #first we need to put from the location column a dataframe with chromosome, start and end position
# list.com=matrix(unlist(strsplit(stats.F1$location,split=":")),ncol=2,byrow=T)
# DMRs.F1.coordinates=data.frame(Chromosome=list.com[,1])
# list.com=matrix(unlist(strsplit(list.com[,2],split="-")),ncol=2,byrow=T)
# DMRs.F1.coordinates$Start=list.com[,1]
# DMRs.F1.coordinates$End=list.com[,2]
# 
# # The function "makeGRangesFromDataFrame" from the library 
# # GenomicRanges makes an object GRanges called "intervals" from "myIntervals"
# DMRs.F1.intervals = GenomicRanges::makeGRangesFromDataFrame(DMRs.F1.coordinates)
# txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
# # extract the list of all gene_Ids from the txdb object
# genes = genes(txdb)
# 
# 
# # Make the annotating function. It will annotate the intervals with gene_Ids
# annotateIntervals <-  function(intervals, txdb)
# {
#   stopifnot(is(intervals, "GRanges"), is(txdb, "TxDb"))
#   anno = genes(txdb)
#   olaps = findOverlaps(intervals, anno)
#   mcols(olaps)$gene_id = genes$gene_id[subjectHits(olaps)]
#   intervals_factor = factor(queryHits(olaps), levels=seq_len(queryLength(olaps)))
#   intervals$gene_id = splitAsList(mcols(olaps)$gene_id, intervals_factor)
#   intervals
# }        
# 
# # Use the "annotateIntervals" funtion in order to annotate 
# #the intervals with gene_Ids and produce "myAnnotation" data.frame 
# myAnnotation.F1.DMRs <- as.data.frame(annotateIntervals(DMRs.F1.intervals, txdb))
# # Make an empty data.frame for append all the annotations, 
# #(we can call it "the master")
# myDf.F1_master <- data.frame()
# 
# 
# # Now we want Hugo gene names in our annotations! 
# #So, for each annotated interval get hugo gene names...
# for (i in 1:length(myAnnotation.F1.DMRs$gene_id)) {
#   # if the gene list is not empty...
#   if(length(c(na.omit(myAnnotation.F1.DMRs$gene_id[i])[[1]])) != 0) {
#     # annotate the interval and copy into a myDf data.frame
#     myDf <- data.frame(myAnnotation.F1.DMRs$seqnames[i], myAnnotation.F1.DMRs$start[i], 
#                        myAnnotation.F1.DMRs$end[i], toString(unname(getSYMBOL(c(na.omit(myAnnotation.F1.DMRs$gene_id[i])[[1]]), data='org.Mm.eg.db'))))
#     # append tge myDF annotations with rbind into the myDf_master
#     myDf.F1_master <- rbind(myDf.F1_master, myDf)
#   }
# }
# 
# # Make a new header for the master dataframe: 
# # it will be the header of our output file also...
# myDf_header <- c("Chromosome", "Start", "End", "Genes")
# names(myDf.F1_master) <- myDf_header
# myDf.F1_master$location <- paste0(myDf.F1_master$Chromosome, ":", myDf.F1_master$Start, "-", myDf.F1_master$End)
# 
# # write myDf_master content in a comma separated output file.
# write.csv(myDf.F1_master, file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/annotated.windows.F1.csv", row.names = FALSE,quote = F)
# 
# # F2 annotation all locations ####
# row.names(stats.F2)=stats.F2[,1]
# #first we need to put from the location column a dataframe with chromosome, start and end position
# list.com=matrix(unlist(strsplit(stats.F2$location,split=":")),ncol=2,byrow=T)
# DMRs.F2.coordinates=data.frame(Chromosome=list.com[,1])
# list.com=matrix(unlist(strsplit(list.com[,2],split="-")),ncol=2,byrow=T)
# DMRs.F2.coordinates$Start=list.com[,1]
# DMRs.F2.coordinates$End=list.com[,2]
# 
# # The function "makeGRangesFromDataFrame" from the library 
# # GenomicRanges makes an object GRanges called "intervals" from "myIntervals"
# DMRs.F2.intervals = GenomicRanges::makeGRangesFromDataFrame(DMRs.F2.coordinates)
# txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
# # extract the list of all gene_Ids from the txdb object
# genes = genes(txdb)
# 
# 
# # Make the annotating function. It will annotate the intervals with gene_Ids
# annotateIntervals <-  function(intervals, txdb)
# {
#   stopifnot(is(intervals, "GRanges"), is(txdb, "TxDb"))
#   anno = genes(txdb)
#   olaps = findOverlaps(intervals, anno)
#   mcols(olaps)$gene_id = genes$gene_id[subjectHits(olaps)]
#   intervals_factor = factor(queryHits(olaps), levels=seq_len(queryLength(olaps)))
#   intervals$gene_id = splitAsList(mcols(olaps)$gene_id, intervals_factor)
#   intervals
# }        
# 
# # Use the "annotateIntervals" funtion in order to annotate 
# #the intervals with gene_Ids and produce "myAnnotation" data.frame 
# myAnnotation.F2.DMRs <- as.data.frame(annotateIntervals(DMRs.F2.intervals, txdb))
# # Make an empty data.frame for append all the annotations, 
# #(we can call it "the master")
# myDf.F2_master <- data.frame()
# 
# 
# # Now we want Hugo gene names in our annotations! 
# #So, for each annotated interval get hugo gene names...
# for (i in 1:length(myAnnotation.F2.DMRs$gene_id)) {
#   # if the gene list is not empty...
#   if(length(c(na.omit(myAnnotation.F2.DMRs$gene_id[i])[[1]])) != 0) {
#     # annotate the interval and copy into a myDf data.frame
#     myDf <- data.frame(myAnnotation.F2.DMRs$seqnames[i], myAnnotation.F2.DMRs$start[i], 
#                        myAnnotation.F2.DMRs$end[i], toString(unname(getSYMBOL(c(na.omit(myAnnotation.F2.DMRs$gene_id[i])[[1]]), data='org.Mm.eg.db'))))
#     # append tge myDF annotations with rbind into the myDf_master
#     myDf.F2_master <- rbind(myDf.F2_master, myDf)
#   }
# }
# 
# # Make a new header for the master dataframe: 
# # it will be the header of our output file also...
# myDf_header <- c("Chromosome", "Start", "End", "Genes")
# names(myDf.F2_master) <- myDf_header
# myDf.F2_master$location <- paste0(myDf.F2_master$Chromosome, ":", myDf.F2_master$Start, "-", myDf.F2_master$End)
# 
# # write myDf_master content in a comma separated output file.
# write.csv(myDf.F2_master, file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/annotated.windows.F2.csv", row.names = FALSE,quote = F)
# 
# # F3 annotation all locations ####
# row.names(stats.F3)=stats.F3[,1]
# #first we need to put from the location column a dataframe with chromosome, start and end position
# list.com=matrix(unlist(strsplit(stats.F3$location,split=":")),ncol=2,byrow=T)
# DMRs.F3.coordinates=data.frame(Chromosome=list.com[,1])
# list.com=matrix(unlist(strsplit(list.com[,2],split="-")),ncol=2,byrow=T)
# DMRs.F3.coordinates$Start=list.com[,1]
# DMRs.F3.coordinates$End=list.com[,2]
# 
# # The function "makeGRangesFromDataFrame" from the library 
# # GenomicRanges makes an object GRanges called "intervals" from "myIntervals"
# DMRs.F3.intervals = GenomicRanges::makeGRangesFromDataFrame(DMRs.F3.coordinates)
# txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
# # extract the list of all gene_Ids from the txdb object
# genes = genes(txdb)
# 
# 
# # Make the annotating function. It will annotate the intervals with gene_Ids
# annotateIntervals <-  function(intervals, txdb)
# {
#   stopifnot(is(intervals, "GRanges"), is(txdb, "TxDb"))
#   anno = genes(txdb)
#   olaps = findOverlaps(intervals, anno)
#   mcols(olaps)$gene_id = genes$gene_id[subjectHits(olaps)]
#   intervals_factor = factor(queryHits(olaps), levels=seq_len(queryLength(olaps)))
#   intervals$gene_id = splitAsList(mcols(olaps)$gene_id, intervals_factor)
#   intervals
# }        
# 
# # Use the "annotateIntervals" funtion in order to annotate 
# #the intervals with gene_Ids and produce "myAnnotation" data.frame 
# myAnnotation.F3.DMRs <- as.data.frame(annotateIntervals(DMRs.F3.intervals, txdb))
# # Make an empty data.frame for append all the annotations, 
# #(we can call it "the master")
# myDf.F3_master <- data.frame()
# 
# 
# # Now we want Hugo gene names in our annotations! 
# #So, for each annotated interval get hugo gene names...
# for (i in 1:length(myAnnotation.F3.DMRs$gene_id)) {
#   # if the gene list is not empty...
#   if(length(c(na.omit(myAnnotation.F3.DMRs$gene_id[i])[[1]])) != 0) {
#     # annotate the interval and copy into a myDf data.frame
#     myDf <- data.frame(myAnnotation.F3.DMRs$seqnames[i], myAnnotation.F3.DMRs$start[i], 
#                        myAnnotation.F3.DMRs$end[i], toString(unname(getSYMBOL(c(na.omit(myAnnotation.F3.DMRs$gene_id[i])[[1]]), data='org.Mm.eg.db'))))
#     # append tge myDF annotations with rbind into the myDf_master
#     myDf.F3_master <- rbind(myDf.F3_master, myDf)
#   }
# }
# 
# # Make a new header for the master dataframe: 
# # it will be the header of our output file also...
# myDf_header <- c("Chromosome", "Start", "End", "Genes")
# names(myDf.F3_master) <- myDf_header
# myDf.F3_master$location <- paste0(myDf.F3_master$Chromosome, ":", myDf.F3_master$Start, "-", myDf.F3_master$End)
# 
# # write myDf_master content in a comma separated output file.
# write.csv(myDf.F3_master, file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/annotated.windows.F3.csv", row.names = FALSE,quote = F)
# 
# ###for the pathway enrichment####
# #first create their Text file format for expression dataset
# #first column with gene names, second with description, the rest with the counts, normalized!
# stats.all.gen=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/meth_countmatrix.txt")
# groups.mus.musculus=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)
# rownames(stats.all.gen) <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# stats.all.gen = stats.all.gen[rowSums(stats.all.gen[,4:59]) >=1,]
# location.stats.all.gen <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# stats.all.gen$location <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$stop)
# #F1####
# windows.annotated.F1=left_join(myDf.F1_master,stats.F1,by="location")
# windows.annotated.F1$Description="NA"
# windows.counts.F1=windows.annotated.F1[,c(4,17,6:16)]
# names(windows.counts.F1)[names(windows.counts.F1) == 'Genes'] <- 'Name'
# write.table(windows.counts.F1,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/windows.counts.F1.txt",row.names = F, quote = F,sep = "\t")
# #F2####
# windows.annotated.F2=left_join(myDf.F2_master,stats.F2,by="location")
# windows.annotated.F2$Description="NA"
# windows.counts.F2=windows.annotated.F2[,c(4,27,6:26)]
# names(windows.counts.F2)[names(windows.counts.F2) == 'Genes'] <- 'Name'
# write.table(windows.counts.F2,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/windows.counts.F2.txt",row.names = F, quote = F,sep = "\t")
# #F3####
# windows.annotated.F3=left_join(myDf.F3_master,stats.F3,by="location")
# windows.annotated.F3$Description="NA"
# windows.counts.F3=windows.annotated.F3[,c(4,30,6:29)]
# names(windows.counts.F3)[names(windows.counts.F3) == 'Genes'] <- 'Name'
# write.table(windows.counts.F3,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/windows.counts.F3.txt",row.names = F, quote = F,sep = "\t")
# 
# #tinny thing, histogram of chr8:20647650-20647763 for F2 as it says that is downregulated and there is a G to C SNP in it in the obese group
# plot(DMRs.allgens.raw.counts[141,13:32])
# F2=data.frame(DMRs.allgens.raw.counts%>%dplyr::select(contains(c("F2"))))
# groups=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)
# groups$sample <- gsub("-", ".", groups$sample)
# # Filter dataframe B for samples with treatment "obese"
# obese_samples <- groups %>% filter(sex == "Obese") %>%  select(sample)
# # Select the relevant columns from dataframe A
# F2.obese <- F2 %>% select(one_of(obese_samples$sample))
# # Filter dataframe B for samples with treatment "Control"
# Control_samples <- groups %>% filter(sex == "Control") %>%  select(sample)
# # Select the relevant columns from dataframe A
# F2.control <- F2 %>% select(one_of(Control_samples$sample))
# F2.obese=F2.obese[141,]
# F2.control=F2.control[141,]
# ggplot()
# #with marta's input let's expand the annotation ####
# 
# 
# #so these is the gen set GSEA, which is not exactly what I want, we want pathways:
# # install.packages("remotes")
# library(remotes)
# # install.packages("devtools")
# library(devtools)
# # install_bitbucket("sonnhammergroup/anubix")
# library(ANUBIX)
# 
# ## Verificacion con MACS3 windows ####
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F3.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F2.rda")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/z.F1.rda")
# z.F1$gen = "F0" 
# z.F2$gen = "F1"
# z.F3$gen = "F2"
# f0.hiper=z.F1[ z.F1$diffexpressed=="Hypermethylated",]
# f1.hipe=z.F2[ z.F2$diffexpressed=="UP",]
# f2.hipe=z.F3[ z.F3$diffexpressed=="UP",]
# f0.hipo=z.F1[ z.F1$diffexpressed=="Hypomethylated",]
# f1.hipo=z.F2[ z.F2$diffexpressed=="DOWN",]
# f2.hipo=z.F3[ z.F3$diffexpressed=="DOWN",]
# all.win.DMR.logFC=rbind(f0.hiper,f1.hipe,f2.hipe,f0.hipo,f1.hipo,f2.hipo)
# write.table(all.win.DMR.logFC,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/merged_DMRs.txt",row.names = F,quote =F,sep = "\t")
# MACS3_win=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/merged.peaks.GBSandGBS-MEDIP.ICR.strain_summits.bed",header=F)
# MACS3_win=data.frame(chr=MACS3_win$V1,start=MACS3_win$V2,end=MACS3_win$V3,name=MACS3_win$V4)
# MACS3_win$location <- paste0(MACS3_win$chr, ":", MACS3_win$start, "-", MACS3_win$end)
# MACS3_win$row=row.names(MACS3_win)
# #verified_wind=inner_join(MACS3_win,all.win.DMR.logFC,by="location")
# #let's try and see if the peak is within the windows####
# library("dplyr")
# all.win.DMR.logFC[c('chrandstart', 'stop')] <- str_split_fixed(all.win.DMR.logFC$location, '-', 2)
# all.win.DMR.logFC[c('chr', 'start')] <- str_split_fixed(all.win.DMR.logFC$chrandstart, ':', 2)
# #verified_wind=apply(as.numeric(MACS3_win$V2),1,function(x) {dplyr::between(x,as.numeric(all.win.DMR.logFC$start),as.numeric(all.win.DMR.logFC$stop))})
# is_in_range <- function(number, df) {
#   #any(as.numeric(number) >= as.numeric(df$start) & number <= as.numeric(df$stop))
#   between(as.numeric(number),df$start,df$stop)
# }
# results1=vector()
# for (i in 1:19) {
#   # Subset both dataframes for the current "chr" value
#   subsetA <- subset(all.win.DMR.logFC, chr == paste0("chr", i))
#   subsetB <- subset(MACS3_win, chr == paste0("chr", i))
#   subsetB$row=row.names(subsetB)
#   subsetB.1=subsetB$start
#   results <- (lapply(subsetB.1, is_in_range, df = subsetA))
#   names(results) <- subsetB$row
#   results1=c(results1,results)
# }
# #taking out all the FALSE values 
# contains_true <- function(x) {
#   any(x == TRUE)
# }
# contains_true_vector <- sapply(results1, contains_true)
# results_TRUE <- results1[contains_true_vector]
# win_verified=data.frame(row=names(results_TRUE))
# verified_windows_DMRs_ICR=inner_join(win_verified,MACS3_win,by="row")
# 
# # Doing a join to see if the peaks from MQ10 GBS-MEDIP are the same as GBS-MEDIP uniquely mapped ####
# MACS3_win_MQ10=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/merged.peaks.GBSandGBS-MEDIP.MQ10.ICR.strain_summits.bed",header=F)
#   MACS3_win_MQ10$location=paste0(MACS3_win_MQ10$V1, ":", MACS3_win_MQ10$V2, "-", MACS3_win_MQ10$V3)
# 
# verified_wind_MQ10=inner_join(MACS3_win_MQ10, verified_windows_DMRs_ICR,by="location")
# verified_wind_MQ10_location=data.frame(location=verified_wind_MQ10$location,name.MQ10=verified_wind_MQ10$V4)
# save(verified_wind_MQ10_location,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/verified_windows_MACS3_MQ10_GBS_MEDIP.rda")
# write.table(verified_wind_MQ10_location,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/verified_windows_GBS-MEDIP_by_MACS3_MQ10.txt",row.names = F,quote =F,sep = "\t")
# load("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/Rresults/verified_windows_MACS3_MQ10_GBS_MEDIP.rda")
# #write the windows in bed format####
# verified_wind_MQ10_location[c('chrandstart', 'stop')] <- str_split_fixed(verified_wind_MQ10_location$location, '-', 2)
# verified_wind_MQ10_location[c('chr', 'start')] <- str_split_fixed(verified_wind_MQ10_location$chrandstart, ':', 2)  
# verified_wind_MQ10_location=verified_wind_MQ10_location[,c(5,6,4,2)]
# write.table(verified_wind_MQ10_location,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/verified_windows_GBS-MEDIP_by_MACS3_MQ10.bed",row.names = F,quote =F,sep = "\t")
# 
# #let's see if the width of the peaks from MACS3 is the same as in the DMRs ####
# narrow_peaks_MACS3_MQ10=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/merged.peaks.GBSandGBS-MEDIP.MQ10.ICR.strain_peaks.narrowPeak",header = F)
# narrow_peaks_MACS3_MQ10$location=paste0(narrow_peaks_MACS3_MQ10$V1, ":", narrow_peaks_MACS3_MQ10$V2, "-", narrow_peaks_MACS3_MQ10$V3)
# narrow_peaks_MACS3_MQ10=data.frame(location=narrow_peaks_MACS3_MQ10$location,name.MQ10=narrow_peaks_MACS3_MQ10$V4)
# location_whole_windows_MACS3=inner_join(narrow_peaks_MACS3_MQ10,verified_wind_MQ10_location,by="name.MQ10")
# colnames(location_whole_windows_MACS3)[colnames(location_whole_windows_MACS3) == "location.x"]="location"
# inner_join(all.win.DMR.logFC,location_whole_windows_MACS3,by="location")
# 
# #now that we have the verified windows, let's do an histogram to check their lenght####
# win.parentals=fread("C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/verification_MACS3/win.size.parentals",header = T)
# freq.win.f2.1=table(win.parentals$F2_1)
# freq.win.f2.2=table(win.parentals$F2_2)
# freq.win.f2.14=table(win.parentals$F2_14)
# freq.win.f2.15=table(win.parentals$F2_15)
# freq.win.f2.16=table(win.parentals$F2_16)
# 
# hist(win.parentals$F2_1, col='red',breaks = 200,xlim = c(0,1800))
# hist(win.parentals$F2_2, col='green',breaks = 200, add=TRUE,xlim = c(0,1800))
# hist(win.parentals$F2_14, col='yellow',breaks = 200, add=TRUE,xlim = c(0,1800))
# hist(win.parentals$F2_15, col='pink',breaks = 200, add=TRUE,xlim = c(0,1800))
# hist(win.parentals$F2_16, col='blue',breaks = 200, add=TRUE,xlim = c(0,1800))
# legend(x= c("F2_1","F2_2","F2_14","F2_15","F2_16"),col = c('red','green','yellow','pink','blue'))
# 
# # let's do the PCA with the verified DMRs ####
# "C:\Users\User\Box\Templeton_mus_musuculus\GBS-MEDIP\counts.verified.wind.MACS3.bed"
# #we see that there are a couple of outliers, F1.8, F2.21 and F1.2, let's exclude them, plotting the 1 and 2 component#####
# load("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/Rresults/normalized.counts.F2.rda")
# load("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/Rresults/normalized.counts.F3.rda")
# load("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/Rresults/normalized.counts.F2.rda")
# load("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/Rresults/normalized.counts.F3.rda")
# #"C:\Users\viode560\Box\Templeton_mus_musuculus\GBS-MEDIP\counts.verified.wind.MACS3.bed"
# win=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/counts.verified.wind.MACS3.txt", header = T)
# win=fread("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/counts.verified.wind.MACS3.txt", header = T)
# win$location <- paste0(win$chr, ":", win$start, "-", win$stop)
# win=as.data.frame(win$location)
# win=as.data.frame(win)
# colnames(win)=c("location")
# groups=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt", header = T)
# groups=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt", header = T)
# colnames(win) <- sub("_", ".", colnames(win))
# verified.windows.MACS3=inner_join(win,stats.F2,by="location")
# #verified.windows.MACS3.1=inner_join(win,stats.F3,by="location")
# colnames(verified.windows.MACS3) <- sub("F2", "F1", colnames(verified.windows.MACS3))
# colnames(verified.windows.MACS3) <- sub("F3", "F2", colnames(verified.windows.MACS3))
# pca <- prcomp(t(verified.windows.MACS3[,-1]), scale. = F,rank. = 3)
# pca.ind=as.data.frame(pca$x)
# groups.DMRs=as.data.frame(row.names(pca$x))
# colnames(groups.DMRs)=c("sample")
# groups=left_join(groups.DMRs,groups,by="sample")
# ggplot(pca.ind,aes(x=PC1,y=PC2,color=interaction(groups$patient,groups$sex)))+geom_point(size=15)+
#   geom_hline(yintercept = 0,linetype="dotted")+labs(y= "PC2", x = "PC1")+
#   geom_vline(xintercept = 0, linetype="dotted")+scale_color_manual(values=c( "#e74269", "#9b1746","#85bc37","#204c26"))+
#   theme(legend.title = element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=50),axis.text.y =element_text(size=50),legend.text = element_text(size=50),axis.title = element_text(size=50))
# 
# ggplot(eigenvectors,aes(x=PC1,y=PC2,color=interaction(groups$patient,groups$sex),label=groups$sample,size=(groups$coverage)))+geom_point()+
#   geom_hline(yintercept = 0)+labs(y= "PC2 (1.93%)", x = "PC1 (2.15%)")+
#   geom_label_repel(show.legend = F,max.overlaps=53)+
#   geom_vline(xintercept = 0)+scale_color_manual(values=c("#d2a9b3", "#e74269", "#9b1746","#85bc37","#204c26"))+
#   theme(legend.title = element_blank(),legend.key=element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=20),axis.text.y =element_text(size=20),legend.text = element_text(size=20),axis.title = element_text(size=20))+
#   guides(color = guide_legend(override.aes = list(size = 10,fill=NA))) 
# 
# 

##########################################################################################################
#DMRs with the new counts from windows from MACS3 and counts from FeatureCounts #####
##########################################################################################################
#stats.all.gen=fread("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/MACS3.final.count.matrix.featureCounts.txt",header = T,fill = T)
stats.all.gen=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/MACS3.final.count.matrix.featureCounts.txt",header = T)
#groups.mus.musculus=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T,fill = T)
groups.mus.musculus=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T,fill = T)
# #index to filter####
# non_zero_counts <- apply(stats.all.gen[,6:60], 1, function(x) sum(x != 0))
# #filter by having at least 5 ind with counts####
# stats.all.gen <- stats.all.gen[non_zero_counts >= 5, ]
#format####
stats.all.gen <- stats.all.gen[!grepl("random", stats.all.gen$end), ]
stats.all.gen <- stats.all.gen %>% select(where(~ all(!is.na(.))))
location.stats.all.gen <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$end)
stats.all.gen$location=paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$end)
individuals.stats.all.gen=as.data.frame(colnames(stats.all.gen[,4:58]))
colnames(individuals.stats.all.gen)=c("sample")
groups.mus.musculus=left_join(individuals.stats.all.gen,groups.mus.musculus,by="sample")

#No stats for the F0 as we only have controls ####
##### mann whitney F1 #####
gen <- "F1"
ind_gen <- groups.mus.musculus$sample[groups.mus.musculus$patient == gen]
group_gen <- groups.mus.musculus$sex[groups.mus.musculus$patient == gen]
F1.stats <- stats.all.gen[, c(..ind_gen,"location")]
row.names(F1.stats)=F1.stats$location
#index to filter####
non_zero_counts <- apply(F1.stats[,1:21], 1, function(x) sum(x != 0))
#filter by having at least 5 ind with counts####
F1.stats <- F1.stats[non_zero_counts >= 2, ]
edgeR.F1 <- DGEList(counts=F1.stats[,1:21], genes=F1.stats$location)
# rownames(edgeR.F1$counts) <- rownames(edgeR.F1$genes) <- F1.stats[,12]
edgeR.F1 <- calcNormFactors(edgeR.F1)
efective.F1= as.data.frame(edgeR.F1$samples$lib.size*edgeR.F1$samples$norm.factors)
edgeR.F1.MW=data.frame(mapply(`*`,F1.stats[,1:21],t(efective.F1)))
p.values.F1.mann.whitney=apply(edgeR.F1.MW, 1, function(x){wilcox.test(x~group_gen,exact=F)$p.value})
p.values.F1.mann.whitney=as.data.frame(p.values.F1.mann.whitney)
colnames(p.values.F1.mann.whitney)=c("p_value")
p.values.F1.mann.whitney$location=F1.stats$location
p.values.F1.mann.whitney$corrected_p_value.BH = p.adjust(((p.values.F1.mann.whitney$p_value)), method = "BH")
p.values.F1.mann.whitney$corrected_p_value = p.adjust(((p.values.F1.mann.whitney$p_value)), method = "bonferroni")
#p.values.F1.mann.whitney=p.values.F1.mann.whitney[p.values.F1.mann.whitney$p_value<0.05,]
#write.csv(p.values.F1.mann.whitney,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/pvalues_testis_F1_MW_less_0.05.csv",quote = F,row.names = F)
# getting the windows that are marked as significant
MW.bonfe.sig.F1=p.values.F1.mann.whitney[p.values.F1.mann.whitney$p_value<0.05,]

#No stats for the F0 as we only have controls ####
##### mann whitney F2 #####
gen <- "F2"
ind_gen <- groups.mus.musculus$sample[groups.mus.musculus$patient == gen]
group_gen <- groups.mus.musculus$sex[groups.mus.musculus$patient == gen]
F2.stats <- stats.all.gen[, c(..ind_gen,"location")]
row.names(F2.stats)=F2.stats$location
#index to filter####
non_zero_counts <- apply(F2.stats[,1:24], 1, function(x) sum(x != 0))
#filter by having at least 5 ind with counts####
F2.stats <- F2.stats[non_zero_counts >= 2, ]
edgeR.F2 <- DGEList(counts=F2.stats[,1:24], genes=F2.stats$location)
# rownames(edgeR.F2$counts) <- rownames(edgeR.F2$genes) <- F2.stats[,12]
edgeR.F2 <- calcNormFactors(edgeR.F2)
efective.F2= as.data.frame(edgeR.F2$samples$lib.size*edgeR.F2$samples$norm.factors)
edgeR.F2.MW=data.frame(mapply(`*`,F2.stats[,1:24],t(efective.F2)))
p.values.F2.mann.whitney=apply(edgeR.F2.MW, 1, function(x){wilcox.test(x~group_gen, exact=F)$p.value})
p.values.F2.mann.whitney=as.data.frame(p.values.F2.mann.whitney)
colnames(p.values.F2.mann.whitney)=c("p_value")
p.values.F2.mann.whitney$location=F2.stats$location
p.values.F2.mann.whitney$corrected_p_value.BH = p.adjust(((p.values.F2.mann.whitney$p_value)), method = "BH")
p.values.F2.mann.whitney$corrected_p_value = p.adjust(((p.values.F2.mann.whitney$p_value)), method = "bonferroni")
#p.values.F2.mann.whitney=p.values.F2.mann.whitney[p.values.F2.mann.whitney$p_value<0.05,]
#write.csv(p.values.F2.mann.whitney,file = "C:/Users/viode560/Documents/templeton_2022_2023_mouse/GBS-MEDIP/pvalues_testis_F2_MW_less_0.05.csv",quote = F,row.names = F)
# getting the windows that are marked as significant
MW.bonfe.sig.F2=p.values.F2.mann.whitney[p.values.F2.mann.whitney$p_value<0.05,]

design.F2=model.matrix(~group_gen)
edgeR.F2 <- estimateDisp(edgeR.F2, design.F2)
x.F2=glmFit(edgeR.F2, design.F2)
x.F2 = glmLRT(x.F2)
y.F2=topTags(x.F2,n=nrow(edgeR.F2))
logFC.F2=data.frame(logFC=y.F2$table$logFC,location=y.F2$table$genes)
z.F2=left_join(p.values.F2.mann.whitney,logFC.F2,by="location")
z.F2$logFC <- -1 * z.F2$logFC
z.F2$diffexpressed <- "NoDiff"
z.F2$diffexpressed[z.F2$logFC > 0.6 & z.F2$p_value < 0.05] <- "Hypermethylated"
z.F2$diffexpressed[z.F2$logFC < -0.6 & z.F2$p_value < 0.05] <- "Hypomethylated"
z.F2$delabel <- NA
z.F2$delabel[z.F2$diffexpressed != "NoDiff"] <- z.F2$location[z.F2$diffexpressed != "NoDiff"]

#let's get the normalized counts for the two gens of the DMRs ####
edgeR.F2.MW$location=F2.stats$location
edgeR.F1.MW$location=F1.stats$location
DMRs.norm.counts=inner_join(edgeR.F2.MW,edgeR.F1.MW,by="location")
DMRs.norm.counts=inner_join(z.F2,DMRs.norm.counts, by="location")

#let's do a dot plot with the meth levels in the Y axis and the different classes on the X axis
individuals.f1.f2.gen=as.data.frame(colnames(DMRs.norm.counts[,8:52]))
colnames(individuals.f1.f2.gen)=c("sample")
groups.mus.musculus.f1.f2=left_join(individuals.f1.f2.gen,groups.mus.musculus,by="sample")
DMRs.co=as.data.frame(t(DMRs.norm.counts[DMRs.norm.counts$diffexpressed!= "NoDiff",]))
DMRs.co=DMRs.co[-c(1:6),]
colnames(DMRs.co)=DMRs.co[1,]
DMRs.co$sample=DMRs.co$delabel
DMRs.co=DMRs.co[-c(1),]
DMRs.co$sample=row.names(DMRs.co)
DMRs.co=inner_join(DMRs.co,groups.mus.musculus.f1.f2,by="sample")
DMRs.co=pivot_longer(DMRs.co, cols = 1:5, names_to = c("window"), values_to = "methylation")

ggplot(DMRs.co , aes(x = window, y = log(as.numeric(methylation)), color = sex, shape= patient)) +
  geom_point(size = 3,position = position_dodge(width = 0.5)) + 
  labs(title = "Dot Plot of Generations by Treatment",
       x = "Generation",
       y = "Methylation level") +
  theme_minimal()+labs(color ="Treatment", shape="Generation")+scale_color_manual(values = c("Control" = "#e74269", "Obese" = "#85bc37"))

#do colors by interaction####
ggplot(DMRs.co , aes(x = window, y = log(as.numeric(methylation)), color = interaction(sex, patient))) +
  geom_point(size = 7,position = position_dodge(width = 0.5)) + 
  labs(x = "Differentially Methylated Region (DMR)",
       y = "Methylation level") +
  theme(legend.title = element_blank(),legend.key=element_blank(),panel.background=element_blank(),
        plot.background=element_blank(),axis.text.x=element_text(size=20),axis.text.y =element_text(size=20),
        legend.text = element_text(size=20),axis.title = element_text(size=30))+
  scale_color_manual(values=c( "#e74269","#85bc37", "#9b1746","#204c26"))

write.table(DMRs.norm.counts, file="C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/DMRs.F1andF2.txt",quote=F,sep="\t",row.names=F)

#### heatmap DMRs ####
#save(DMRs.norm.counts,file = "C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/DMRs.F1andF2.rda")
load("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/DMRs.F1andF2.rda")
groups.mus.musculus=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)
load("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/DMRs.F1andF2.rda")
groups.mus.musculus=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)

DMRs.norm.counts=DMRs.norm.counts[DMRs.norm.counts$diffexpressed!= "NoDiff",]
#write.table(DMRs.norm.counts,file = "C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/DMRs_final_intragenerational.txt",sep = "\t",quote = F, row.names = F)
locations_wind=DMRs.norm.counts$location
DMRs.norm.counts=DMRs.norm.counts[,8:52]
indv=data.frame(sample=names(DMRs.norm.counts))
groups=left_join(indv,groups.mus.musculus,by="sample")
groups <- groups %>%
  arrange(patient, sex, family, sample)
DMRs.norm.counts <- DMRs.norm.counts[,groups$sample]
row.names(DMRs.norm.counts)=locations_wind
DMRs.norm.counts=log(DMRs.norm.counts)
DMRs.norm.counts[DMRs.norm.counts==-Inf] <- 0
DMRs.norm.counts=as.matrix(DMRs.norm.counts)

# groups.F2 = groups.mus.musculus[groups.mus.musculus$patient == "F2", ]
# group.F2.obese = groups.F2[groups.F2$sex == "Obese", ]
# group.F2.control = groups.F2[groups.F2$sex == "Control", ]
# F2.o=group.F2.obese$sample
# F2.c=group.F2.control$sample
# DMRs.F2.obese.scaled <- DMRs.norm.counts.1[, colnames(DMRs.norm.counts.1) %in% F2.o]
# DMRs.F2.control.scaled <- DMRs.norm.counts.1[, colnames(DMRs.norm.counts.1) %in% F2.c]
# groups.F1 = groups.mus.musculus[groups.mus.musculus$patient == "F1", ]
# group.F1.obese = groups.F1[groups.F1$sex == "Obese", ]
# group.F1.control = groups.F1[groups.F1$sex == "Control", ]
# F1.o=group.F1.obese$sample
# F1.c=group.F1.control$sample
# DMRs.F1.obese.scaled <- DMRs.norm.counts.1[, colnames(DMRs.norm.counts.1) %in% F1.o]
# DMRs.F1.control.scaled <- DMRs.norm.counts.1[, colnames(DMRs.norm.counts.1) %in% F1.c]

# DMRs.F2.scaled=cbind(DMRs.F1.control.scaled,DMRs.F1.obese.scaled,DMRs.F2.control.scaled,DMRs.F2.obese.scaled)
# row.names(DMRs.F2.scaled)=DMRs.norm.counts.1$location
# DMRs.f2=as.matrix(DMRs.F2.scaled)
# par(mar=c(1,1,1,1))
# DMRs.f0.1=log(DMRs.f2)
# DMRs.f0.1[DMRs.f0.1==-Inf] <- 0
# par(cex.axis=0.8)
# tiff("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/heatmap.F2.DMR.GBS.MeDIP.log.normalizedcounts.tiff", width = 100000, height = 10000, units = "px", res = 1000)
# heatmap.2(DMRs.f0.1,dendrogram='none', col = brewer.pal(n = 9, name = "Reds"),key.title ="",key=T ,key.ylab = "" ,key.xlab = "log(normalized counts)",margins=c(5,20),Rowv=FALSE, Colv=F,trace='none',labCol = colnames(DMRs.f2),scale ="none")
# dev.off()

######## Multilayered heatmap to plot methylation levels and depth ####
# depth_DMRs=read.table("C:/Users/viode560/Box/Templeton_mus_musuculus/info_templeton_proyect/DMRs_depth.txt",header = T,fill = T)
# row.names(depth_DMRs)=depth_DMRs$location
# pa="F0"
# columns_to_keep <- !grepl(pa, names(depth_DMRs))
# depth_DMRs <- depth_DMRs[, columns_to_keep]
# depth_DMRs=depth_DMRs[,-1]
# family_map <- groups %>%
#   select(sample, patient, family) %>%
#   arrange(patient, family, sample) %>%
#   mutate(col_order = row_number()) %>%
#   with(setNames(col_order, sample))
# ordered_cols <- names(depth_DMRs)[order(family_map[names(depth_DMRs)])]
# depth_DMRs <- depth_DMRs[ordered_cols]

# #do heatmap for meth levels####
# col_fun = colorRamp2(c(0, 15.43), c("#fff3f3", "#990000"))
# column_ha = HeatmapAnnotation(Generation=anno_block(labels = c("F1","F1","F2","F2"),labels_gp = gpar(col = "black", fontsize = 20)),Group=anno_block(labels = c("Control","Obese","Control","Obese"),labels_gp = gpar(col = "black", fontsize = 20)),annotation_legend_param = list(title_size = 20, labels_size = 20))
# Heatmap_methlevels=Heatmap(DMRs.f0.1, name = "log(Methylation levels)",
#                            bottom_annotation  = column_ha,border = TRUE,column_gap = unit(1, "mm"), 
#                            col = col_fun,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),
#                            column_title ="Methylation level",show_column_names = FALSE,row_title = "Genomic locations",
#                            cluster_rows = FALSE,cluster_columns = FALSE,
#                            heatmap_legend_param = list(title_size = 20, labels_size = 20,legend_position = "topleft"),
#                            row_title_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20),
#                            row_names_gp = gpar(fontsize = 20),cell_fun = function(j, i, x, y, width, height, fill) {
#                              grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "grey", lwd = 0.5))})
# #column_ha = HeatmapAnnotation(Generation=c(rep("F1",21),rep("F2",24)),Group=c(rep("Control",10),rep("Obese",11),rep("Control",12),rep("Obese",12)),col = list(Group = c("Obese"="#85bc37","Control"="#e74269"),Generation= c("F1"="grey","F2"="black")))                                                    
# 
# #do heatmap for depth level GBS-MEDIP ####
# depth_DMRs=as.matrix(depth_DMRs[,-1])
# depth_DMRs <- depth_DMRs[, colnames(DMRs.f0.1)]
# # col_fun = colorRamp2(c(0,7.28,14.56,29.125,58.25,111.5, 223), c("#FBFDFB","#BFE1B0","#74C67A","#39A96B","#1D9A6C", "#137177","#0A2F51"))
# # heatmap_depth=Heatmap(depth_DMRs, name = "Average depth per site",bottom_annotation  = column_ha,border = TRUE,column_gap = unit(1, "mm"), col = col_fun,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),column_title ="Individuals",row_title = "Genomic locations",cluster_rows = FALSE,cluster_columns = FALSE)
# # 
# # Heatmap_methlevels+heatmap_depth ####
# depth_DMRs=log(depth_DMRs)
# col_fun1 = colorRamp2(c(0,5.41), c("#FBFDFB","#0A2F51"))
# heatmap_depth=Heatmap(depth_DMRs, name = "log(Average depth per site)",bottom_annotation  = column_ha,border = TRUE,
#                       column_gap = unit(1, "mm"), col = col_fun1,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),
#                       column_title ="Average depth per site",show_column_names = FALSE,row_title = "Genomic locations",
#                       cluster_rows = FALSE,cluster_columns = FALSE,
#                       heatmap_legend_param = list(title_size = 20, labels_size = 20,legend_position = "topleft"),
#                       row_title_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20),
#                       row_names_gp = gpar(fontsize = 20))
# ht_list=Heatmap_methlevels+heatmap_depth
# par(mar=c(1,1,1,1))
# tiff("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/DMRs_and_depth.tiff", width = 30000, height = 10000, units = "px", res = 1000)
# draw(ht_list, heatmap_legend_side = "left", padding = unit(c(4, 3, 3, 29), "mm"))
# dev.off()

#do heatmap for depth level GBS ####
depth_DMRs=read.table("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/Quality_control/DMRs_depth_GBS.txt",header = T,fill = T)
depth_DMRs=read.table("C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/Quality_control/DMRs_depth_GBS.txt",header = T,fill = T)
row.names(depth_DMRs)=depth_DMRs$location
pa="F0"
columns_to_keep <- !grepl(pa, names(depth_DMRs))
depth_DMRs <- depth_DMRs[, columns_to_keep]
groups <- groups %>%
  arrange(patient, sex, family, sample)
depth_DMRs <- depth_DMRs[,groups$sample]
depth_DMRs=as.matrix(depth_DMRs[,-46])
# depth_DMRs=log(depth_DMRs)
# depth_DMRs[depth_DMRs==-Inf] <- 0
# depth_DMRs <- depth_DMRs[, colnames(DMRs.f0.1)]
# col_fun = colorRamp2(c(0,7.28,14.56,29.125,58.25,111.5, 223), c("#FBFDFB","#BFE1B0","#74C67A","#39A96B","#1D9A6C", "#137177","#0A2F51"))
# heatmap_depth=Heatmap(depth_DMRs, name = "Average depth per site",bottom_annotation  = column_ha,border = TRUE,column_gap = unit(1, "mm"), col = col_fun,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),column_title ="Individuals",row_title = "Genomic locations",cluster_rows = FALSE,cluster_columns = FALSE)
# 
# Heatmap_methlevels+heatmap_depth ####
#depth_DMRs=log(depth_DMRs)
# col_fun1 = colorRamp2(c(0,10.8), c("#FBFDFB","#0A2F51"))
# heatmap_depth=Heatmap(depth_DMRs, name = "log(Average depth per site GBS)",bottom_annotation  = column_ha,border = TRUE,
#                       column_gap = unit(1, "mm"), col = col_fun1,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),
#                       column_title ="Average depth per site in GBS",show_column_names = FALSE,row_title = "Genomic locations",
#                       cluster_rows = FALSE,cluster_columns = FALSE,
#                       heatmap_legend_param = list(title_size = 20, labels_size = 20,legend_position = "topleft"),
#                       row_title_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20),
#                       row_names_gp = gpar(fontsize = 20))
# ht_list=Heatmap_methlevels+heatmap_depth
# save(ht_list,file = "C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/heatmap_GBSdepth_and_DMRs.rda")
# par(mar=c(1,1,1,1))
# tiff("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/DMRs_and_depthGBS.tiff", width = 30000, height = 10000, units = "px", res = 1000)
# draw(ht_list, heatmap_legend_side = "left", padding = unit(c(4, 3, 3, 29), "mm"))
# dev.off()
############################################################################
# Heatmap with the family info and above GBS-MEDIP y abajo GBS depth #####
# col_fun = colorRamp2(c(0, 15.43), c("#fff3f3", "#990000"))
# families_ICR=c(" "," "," ","2.3.3"," ","2.3.3"," ","3.6.4","2.3.3","4.8.4"," "," "," "," "," "," "," ","5.4.3","5.4.1","8.1.3","10.1.5","4.8.4","4.8.4","4.8.4","3.6.4","3.6.4","3.6.4","3.6.4"," ","4.8.4","4.8.4","4.8.4","3.6.4","5.4.3","8.1.3","8.1.3","5.4.3","5.4.1","5.4.1","5.4.1","5.4.3","8.1.3","10.1.5","10.1.5","")
# align_vector <- seq_along(families_ICR) 
# column_ha = HeatmapAnnotation(Family= anno_text(families_ICR))
# # column_ha = HeatmapAnnotation(
# #   Family = anno_textbox(
# #     text = families_ICR,
# #     align_to = align_vector,background_gp = gpar(fill = "white", col = "white"),
# #     gp = gpar(fontsize = 12, fontface = "bold", lwd = 1), 
# #     just = "center")
# #   )
# Heatmap_methlevels.1=Heatmap(DMRs.f0.1, name = "log(Methylation levels)",height = unit(7, "mm")*5,width = unit(10, "mm")*56,
#                            top_annotation  = column_ha,border = TRUE,column_gap = unit(1, "mm"), 
#                            col = col_fun,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),
#                            column_title ="Methylation level",show_column_names = FALSE,row_title = "",
#                            cluster_rows = FALSE,cluster_columns = FALSE,
#                            heatmap_legend_param = list(title_size = 20, labels_size = 20,legend_position = "topleft"),
#                            row_title_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20),
#                            row_names_gp = gpar(fontsize = 20),cell_fun = function(j, i, x, y, width, height, fill) {
#                              grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "grey", lwd = 0.5))})
# col_fun1 = colorRamp2(c(0,10.8), c("#FBFDFB","#0A2F51"))
# column_ha.1 = HeatmapAnnotation(Generation=anno_block(labels = c("F1","F1","F2","F2"),labels_gp = gpar(col = "black", fontsize = 20)),Group=anno_block(labels = c("Control","Obese","Control","Obese"),labels_gp = gpar(col = "black", fontsize = 20)),annotation_legend_param = list(title_size = 20, labels_size = 20))
# heatmap_depth.1=Heatmap(depth_DMRs, name = "log(Average depth per site GBS)",bottom_annotation  = column_ha.1,border = TRUE,height = unit(7, "mm")*5,width = unit(10, "mm")*56,
#                       column_gap = unit(1, "mm"), col = col_fun1,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),
#                       column_title ="Average depth per site in GBS",show_column_names = FALSE,row_title = "Genomic locations",
#                       cluster_rows = FALSE,cluster_columns = FALSE,
#                       heatmap_legend_param = list(title_size = 20, labels_size = 20,legend_position = "topleft"),
#                       row_title_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20),
#                       row_names_gp = gpar(fontsize = 20))
# ht_list=Heatmap_methlevels.1%v%heatmap_depth.1
# draw(ht_list)
# par(mar=c(1,1,1,1))
# tiff("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/heatmap_families_with_numbers.tiff", width = 30000, height = 10000, units = "px", res = 1000)
# draw(ht_list, heatmap_legend_side = "left", padding = unit(c(4, 3, 3, 29), "mm"))
# dev.off()

#################################################################################
# Now create the legend with colors, which will be better ######
col_fun = colorRamp2(c(0, 15.43), c("#fff3f3", "#990000"))
families_ICR=as.numeric(factor(groups$family))
colors <- rainbow(21)
col_fun.4 = colorRamp2(1:21, colors)
column_ha.3 = HeatmapAnnotation(Family= families_ICR, col = list(Family = col_fun.4),show_legend = FALSE )
Heatmap_methlevels.1=Heatmap(DMRs.norm.counts, name = "log(Methylation levels)",height = unit(10, "mm")*5,width = unit(10, "mm")*45,
                             top_annotation  = column_ha.3,border = TRUE,column_gap = unit(1, "mm"), 
                             col = col_fun,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),
                             row_title  ="",show_column_names = FALSE,column_title = "Methylation level",
                             cluster_rows = FALSE,cluster_columns = FALSE,
                             heatmap_legend_param = list(title_size = 20, labels_size = 20),
                             row_title_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20),
                             row_names_gp = gpar(fontsize = 20),cell_fun = function(j, i, x, y, width, height, fill) {
                               grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "grey", lwd = 0.5))})
tiff("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/heatmap_families_with_color.tiff", width = 5000, height = 1000, units = "px", res = 200)
Heatmap_methlevels.4=Heatmap(DMRs.norm.counts, name = "log(Methylation levels)",height = unit(10, "mm")*5,width = unit(10, "mm")*45,
                             top_annotation  = column_ha.3,bottom_annotation  = column_ha.1,border = TRUE,column_gap = unit(1, "mm"), 
                             col = col_fun,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),
                             row_title  ="",show_column_names = FALSE,column_title = "Methylation level",
                             cluster_rows = FALSE,cluster_columns = FALSE,
                             heatmap_legend_param = list(title_size = 20, labels_size = 20),
                             row_title_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20),
                             row_names_gp = gpar(fontsize = 20),cell_fun = function(j, i, x, y, width, height, fill) {
                               grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "grey", lwd = 0.5))})

draw(Heatmap_methlevels.4,heatmap_legend_side = "left",padding = unit(c(4, 3, 3, 29), "mm"))
dev.off()
col_fun1 = colorRamp2(c(0,10.8), c("#FBFDFB","#0A2F51"))
column_ha.1 = HeatmapAnnotation(Generation=anno_block(labels = c("F1","F1","F2","F2"),labels_gp = gpar(col = "black", fontsize = 20)),Group=anno_block(labels = c("Control","Overnutrition","Control","Overnutrition"),labels_gp = gpar(col = "black", fontsize = 20)),annotation_legend_param = list(title_size = 20, labels_size = 20))
heatmap_depth.1=Heatmap(depth_DMRs, name = "Average depth per site GBS",bottom_annotation  = column_ha.1,border = TRUE,height = unit(10, "mm")*5,width = unit(10, "mm")*45,
                        column_gap = unit(1, "mm"), col = col_fun1,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),
                        row_title  ="Average depth\nper site in GBS",show_column_names = FALSE,column_title = "",
                        cluster_rows = FALSE,cluster_columns = FALSE,
                        heatmap_legend_param = list(title_size = 20, labels_size = 20),
                        row_title_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20),
                        row_names_gp = gpar(fontsize = 20))
ht_list=Heatmap_methlevels.1%v%heatmap_depth.1

draw(ht_list)
par(mar=c(1,1,1,1))
tiff("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/heatmap_families_with_color.tiff", width = 30000, height = 10000, units = "px", res = 1000)
tiff("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/heatmap_families_with_color_and_depth.tiff", width = 30000, height = 10000, units = "px", res = 1000)
draw(ht_list, heatmap_legend_side = "left", padding = unit(c(4, 3, 3, 29), "mm"))
dev.off()

# Annotation of genomic features GBS-MEDIP ####
gbs_MEDIP_annot=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/info_templeton_proyect/annotation_GBS_medip.txt")
gbs_MEDIP_annot=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/annotation_GBS_medip.txt")
gbs_annot=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/info_templeton_proyect/annotation_GBS.txt")
gbs_annot=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/annotation_GBS.txt")

gbs_MEDIP_annot_genomicF=gbs_MEDIP_annot[1:6,-c(2)]
gbs_MEDIP_annot_TE=gbs_MEDIP_annot[7:13,-c(2)]
gbs_annot_genomicF=gbs_annot[1:6,-c(2)]
gbs_annot_TE=gbs_annot[7:13,-c(2)]


gbs_MEDIP_annot_genomicF=pivot_longer(gbs_MEDIP_annot_genomicF,cols = 2:5,names_to = "Generation", values_to = "Percentage")
gbs_MEDIP_annot_TE=pivot_longer(gbs_MEDIP_annot_TE,cols = 2:5,names_to = "Generation", values_to = "Percentage")

gbs_annot_genomicF=pivot_longer(gbs_annot_genomicF,cols = 2:5,names_to = "Generation", values_to = "Percentage")
gbs_annot_TE=pivot_longer(gbs_annot_TE,cols = 2:5,names_to = "Generation", values_to = "Percentage")
gbs_annot_TE$log=log(gbs_annot_TE$Percentage)

gbs_GF=ggplot(gbs_annot_genomicF, aes(x = Generation, y = Percentage, fill =Genomic_feature )) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Percentage", title = "GBS genomic feature\ndistribution") +
  scale_fill_manual(values = c("#a559aa", "#59a89c", "#f0c571", "#cecece", "#e02b35", "#082a54")) +
  theme_minimal()+theme(legend.position = "none")+
  scale_y_continuous(labels = scales::percent_format())

gbs_MEDIP_GF=ggplot(gbs_MEDIP_annot_genomicF, aes(x = Generation, y = Percentage, fill =Genomic_feature )) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "", title = "GBS-MeDIP genomic\nfeature distribution") +
  scale_fill_manual(values = c("#a559aa", "#59a89c", "#f0c571", "#cecece", "#e02b35", "#082a54")) +
  theme_minimal()+
  scale_y_continuous(labels = scales::percent_format())
combined_plot <- gbs_GF + gbs_MEDIP_GF + plot_layout(ncol = 2)
library(patchwork)
tiff("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/Annotation/GBS_and_GBS_MEDIP_GF.tiff",width = 2000, height = 1000, units = "px",res =200)
tiff("C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/Annotation/GBS_and_GBS_MEDIP_GF.tiff",width = 2500, height = 1000, units = "px",res =200)
combined_plot
dev.off()

custom_labels <- function(x) {
  paste0(round(exp(x), 1), "%")
}

gbs_annot_TE$Percentage <- ifelse(gbs_annot_TE$Percentage == 0, 0.1, gbs_annot_TE$Percentage)
gbs_TE=ggplot(gbs_annot_TE, aes(x = Generation, y = Percentage, fill = Genomic_feature )) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Percentage with the axis in log scale", title = "GBS transposable element distribution") +
  scale_fill_manual(values = c("#008B8B", "#FF7F50", "#5F9EA0", "#DAA520", "#6B8E23", "#708090", "#B22222")) +
  theme_minimal()+
  scale_y_log10()

library(scales) 
gbs_TE=ggplot(gbs_annot_TE, aes(x = Generation, y = Percentage, fill = Genomic_feature )) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Percentage with the axis in log scale", title = "GBS transposable element distribution") +
  scale_fill_manual(values = c("#008B8B", "#FF7F50", "#5F9EA0", "#DAA520", "#6B8E23", "#708090", "#B22222")) +
  theme_minimal()+
  scale_y_continuous(
    trans = 'log10',
    labels = percent_format(accuracy = 1),
    breaks = 10^(-3:2)  # Breaks at 0.1%, 1%, 10%, 100%
  ) +
  coord_cartesian(ylim = c(0.1, 100))

gbs_TE=ggplot(gbs_annot_TE, aes(x = Generation, y = Percentage, fill = Genomic_feature )) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Percentage", title = "GBS transposable element distribution") +
  scale_fill_manual(values = c("#008B8B", "#FF7F50", "#5F9EA0", "#DAA520", "#6B8E23", "#708090", "#B22222")) +
  theme_minimal()+theme(legend.position = "none")+
  scale_y_continuous(labels = scales::percent_format())

gbs_MEDIP_TE=ggplot(gbs_MEDIP_annot_TE, aes(x = Generation, y = Percentage, fill = Genomic_feature )) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Percentage", title = "GBS-MEDIP transposable element distribution") +
  scale_fill_manual(values = c("#008B8B", "#FF7F50", "#5F9EA0", "#DAA520", "#6B8E23", "#708090", "#B22222")) +
  theme_minimal()+
  scale_y_continuous(labels = scales::percent_format())
combined_plot.1 <- gbs_TE + gbs_MEDIP_TE + plot_layout(ncol = 2)
tiff("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/Annotation/GBS_and_GBS_MEDIP_RE.tiff",width = 2500, height = 1000, units = "px",res =200)
tiff("C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/Annotation/GBS_and_GBS_MEDIP_RE.tiff",width = 2500, height = 1000, units = "px",res =200)
combined_plot.1
dev.off()

# PCA with the new counts and with group and family as variables ####
# Necesito una matrix con las ventanas x columnas, y el ID de cada familia en una columna (aqui no importan individuos)
stats.all.gen=fread("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/MACS3.final.count.matrix.featureCounts.txt",header = T,fill = T)
stats.all.gen=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/MACS3.final.count.matrix.featureCounts.txt",header = T)
groups.mus.musculus=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T,fill = T)
groups.mus.musculus=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T,fill = T)
stats.all.gen <- stats.all.gen[!grepl("random", stats.all.gen$end), ]
stats.all.gen <- stats.all.gen %>% select(where(~ all(!is.na(.))))
location.stats.all.gen <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$end)
stats.all.gen$location=paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$end)
individuals.stats.all.gen=as.data.frame(colnames(stats.all.gen[,4:58]))
colnames(individuals.stats.all.gen)=c("sample")
groups.mus.musculus=left_join(individuals.stats.all.gen,groups.mus.musculus,by="sample")

normalized_win <- DGEList(counts=stats.all.gen[,4:48], genes=location.stats.all.gen)
# rownames(edgeR.F1$counts) <- rownames(edgeR.F1$genes) <- F1.stats[,12]
normalized_win <- calcNormFactors(normalized_win)
efective.all= as.data.frame(normalized_win$samples$lib.size*normalized_win$samples$norm.factors)
normalized_count_all=data.frame(mapply(`*`,stats.all.gen[,4:48],t(efective.all)))
normalized_count_all=as.data.frame(t(normalized_count_all))
colnames(normalized_count_all)=location.stats.all.gen
normalized_count_all=rownames_to_column(normalized_count_all,"sample")
normalized_count_all=left_join(normalized_count_all,groups.mus.musculus,by="sample")
rownames(normalized_count_all)=normalized_count_all$sample
normalized_count_all$patient=as.factor(normalized_count_all$patient)
normalized_count_all$sex=as.factor(normalized_count_all$sex)
normalized_count_all$family=as.factor(normalized_count_all$family)
#normalized_count_all[1:45,2:1600]=as.numeric(normalized_count_all[1:45,2:1600])
#save(normalized_count_all,file = "C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/normalized_counts_and_metadata.rda")
load("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/normalized_counts_and_metadata.rda")
load("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/normalized_counts_and_metadata.rda")

#scale the data #### F1.6 and F2.19 are clear outliers in PC1, PC2 and PC3
# normalized_count_all <- normalized_count_all %>%  filter(!str_detect(sample, "^F1\\.2$"))
# normalized_count_all <- normalized_count_all %>%  filter(!str_detect(sample, "^F2\\.16$"))
# normalized_count_all <- normalized_count_all %>%  filter(!str_detect(sample, "^F1\\.5$"))
normalized_count_all <- normalized_count_all %>%  filter(!str_detect(sample, "^F1\\.6$"))
normalized_count_all <- normalized_count_all %>%  filter(!str_detect(sample, "^F2\\.19$"))
# normalized_count_all <- normalized_count_all %>%  filter(!str_detect(sample, "^F2\\.1$"))
normalized_count_all.1=normalized_count_all[,c(2:1600)]
normalized_count_all.2=normalized_count_all[,c(2:1600)]
normalized_count_all.3=normalized_count_all[,c(2:1600)]
normalized_count_all.2$family=normalized_count_all$family
normalized_count_all.3$group=normalized_count_all$sex
normalized_count_all.3$family=normalized_count_all$family

#gbs_medip_pca <- dudi.pca((normalized_count_all[,2:1600]),scannf = FALSE, nf = 5,scale = TRUE)
#gbs_medip_pca.1=PCA(normalized_count_all,ncp = 5, graph = F,quali.sup =c(1,1601,1602,1603,1604))
gbs_medip_pca.1=PCA(normalized_count_all.1,ncp = 5, graph = F,quali.sup = 1600)
gbs_medip_pca.2=PCA(normalized_count_all.2,ncp = 5, graph = F,quali.sup = 1600)
gbs_medip_pca.3=PCA(normalized_count_all.3,ncp = 5, graph = F,quali.sup =c(1600,1601))
normalized_count_all.1$group=normalized_count_all$sex
normalized_count_all.1$group <- gsub("Obese", "Overnutrition", normalized_count_all.1$group)
normalized_count_all.1$gen=normalized_count_all$patient
normalized_count_all.1$family=normalized_count_all$family
dim.fam.group.gbsmedip=dimdesc(gbs_medip_pca.3, axes=c(1,2),proba = 1)
#test for normality ####
coord.ind=as.data.frame(gbs_medip_pca.3[["ind"]][["coord"]])
coord.ind$sample=row.names(coord.ind)
coord.ind=left_join(coord.ind,groups.mus.musculus,by="sample")
ks.test(coord.ind$Dim.1,"pnorm")
ks.test(coord.ind$Dim.2,"pnorm")
kruskal.test(Dim.1 ~ family, data = coord.ind)
kruskal.test(Dim.2 ~ family, data = coord.ind)
kruskal.test(Dim.1 ~ sex, data = coord.ind.grou)
kruskal.test(Dim.2 ~ sex, data = coord.ind.grou)

#see the percentage of explained variance by each component ####
fviz_eig(gbs_medip_pca)
fviz_eig(gbs_medip_pca.1,addlabels = TRUE)

get_eigenvalue(gbs_medip_pca.1)

# Do PCA for the windows ####
win_pca <- get_pca_var(gbs_medip_pca.1)
fviz_pca_var(gbs_medip_pca.1, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
# only family ####
p=fviz_pca_var(gbs_medip_pca.2, col.var = "contrib",geom = "arrow",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),select.var = list(contrib = 150))
p <- fviz_add(p, gbs_medip_pca.2$quali.sup$coord, color = "black")
# p = p + ylim(c(-5,2.5))
# p = p + xlim(c(-3,2.5))
p = p + labs(x="PC1 (6%)",y="PC2 (5.2%)")+guides(color = "none")+ggtitle("Display of variable family with GBS-MEDIP data")
p
######
p=fviz_pca_var(gbs_medip_pca.2, col.var = "contrib",geom = "arrow",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),select.var = list(contrib = 150))
#p <- fviz_add(p, gbs_medip_pca.2$quali.sup$coord, color = "black")
#group=c("C","C","O","O","O","O","C","C","C","C","C","C","C","C","O","O","O","O","O","O","O")
group=c("red","red","blue","blue","blue","blue","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue")
p <- fviz_add(p, gbs_medip_pca.2$quali.sup$coord, color = factor(group))
p = p + ylim(c(-5,2.5))
p = p + xlim(c(-3,2.5))
p = p + labs(x="PC1 (6%)",y="PC2 (5.2%)")+guides(color = "none")+ggtitle("Display of variable family with GBS-MEDIP data")
p
# only groups ####
d=fviz_pca_var(gbs_medip_pca.1, col.var = "contrib",geom = "arrow",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),select.var = list(contrib = 150))
d <- fviz_add(d, gbs_medip_pca.1$quali.sup$coord, color = "black")
d=d+labs(x="PC1 (5.3%)",y="PC2 (4.7%)",color = "Contribution\nof the windows\nto the PCs")+ ggtitle("Display of variable treatment")
#d = d + xlim(c(-1.5,1.5))
d

# Have those both plots one next to the other ####
combined_plot = p + d 
tiff(filename = "C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/PCA_windows_contribution_family_treatment.tiff",
     height = 1200,width = 2000,res = 150)
tiff(filename = "C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/PCA_windows_contribution_family_treatment.tiff",
     height = 1200,width = 3500,res = 150)
print(combined_plot)
dev.off()
#only groups and the windows that contribute the most ####
d=fviz_pca_var(gbs_medip_pca.1, col.var = "contrib",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),select.var = list(contrib = 150))
d <- fviz_add(d, gbs_medip_pca.1$quali.sup$coord, color = "black")
d
select.ind = list(contrib = 5)

#p=fviz_pca_var(gbs_medip_pca.2)

# fviz_add(p, gbs_medip_pca.2$quali.sup$coord, color = c("red"))
# fviz_add(p, gbs_medip_pca.2$quanli.sup$coord, 
#          geom = c("arrow", "text"), 
#          color = "red")
#fviz_pca_var(gbs_medip_pca.1, alpha.var = "cos2")

#Controbution of the windows to the PC ####
library("corrplot")

corrplot(win_pca$contrib, is.corr=FALSE) 

# Do the PCA for individuals ####
fviz_pca_ind(gbs_medip_pca.1,axes = c(1,2))
fviz_pca_ind(gbs_medip_pca.1,axes = c(2,3))
# e <- fviz_pca_ind(gbs_medip_pca.1, 
#                   axes = c(2, 3), 
#                   habillage = interaction(normalized_count_all.1$group, normalized_count_all.1$gen),
#                   palette = c("#e74269", "#85bc37", "#9b1746", "#204c26"),
  #                 pointsize = 5,  # Set point size
  #                 labelsize = 10,  
  #                 mean.point = FALSE) +  
  # ylim(c(-5, 2.5)) +  # Set y-axis limits
  # xlim(c(-3, 1)) +  # Set x-axis limits
  # labs(x = "PC2 (4.7%)", y = "PC3 (4.6%)") +  # Set axis labels
  # theme(legend.position = "right",  # Customize legend position
  #       legend.title = element_blank(),  # Remove legend title
  #       legend.text = element_text(size = 20),  # Adjust legend text size
  #       axis.text = element_text(size = 20),  # Adjust axis text size
  #       axis.title = element_text(size = 20)) +  # Adjust axis title size
  # guides(color = guide_legend(override.aes = list(size = 5))) +  # Set point size in legend
  # ggtitle(NULL)  # Remove title

########################
#put the family labels
gbs_medip_eigenvectors=as.data.frame(gbs_medip_pca.1$ind$coord)
gbs_medip_eigenvectors$group=normalized_count_all.1$group
gbs_medip_eigenvectors$gen=normalized_count_all.1$gen
gbs_medip_eigenvectors$family=normalized_count_all.1$family
PCA_noF0_gbs_medip=ggplot(gbs_medip_eigenvectors,aes(x=Dim.2,y=Dim.3,color=interaction(group, gen),
                                                     label=gbs_medip_eigenvectors$family))+
  geom_point(size=7)+ylim(c(-5, 2.5))+xlim(c(-3, 1))

PCA_noF0_gbs_medip=PCA_noF0_gbs_medip+geom_hline(yintercept = 0)+
  labs(x = "PC1 (4.7%)", y = "PC2 (4.6%)")+
  geom_vline(xintercept = 0)+
  scale_color_manual(values=c("#e74269", "#85bc37","#9b1746","#204c26"))

PCA_noF0_gbs_medip=PCA_noF0_gbs_medip+theme(legend.title = element_blank(),legend.key=element_blank(),
                                            panel.background=element_blank(),
                                            plot.background=element_blank(),
                                            axis.text.x=element_text(size=20),
                                            axis.text.y =element_text(size=20),
                                            legend.text = element_text(size=20),
                                            axis.title = element_text(size=20))+
  guides(color = guide_legend(override.aes = list(size = 10,fill=NA)))

PCA_noF0_gbs_medip=PCA_noF0_gbs_medip+geom_label_repel(show.legend = F,max.overlaps=43)

#put the plot screen to the max or it will give error!!!!!
tiff(filename = "C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/individuals_PCA_all_meth_win.tiff",width = 3000, height = 2000, units = "px",res =200)
tiff(filename = "C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/individuals_PCA_all_meth_win.tiff",width = 3000, height = 2000, units = "px",res =200)
print(PCA_noF0_gbs_medip)
dev.off()

a=fviz_pca_ind(gbs_medip_pca.1,axes = c(2,3),
             repel = TRUE,geom.ind = "point",
             fill.ind = normalized_count_all.1$group,pointshape = 21,pointsize=4
)
a=a+xlim(c(-4,1.5))
a=a+ylim(c(-10,2.5))
a <- fviz_add(a, gbs_medip_pca.1$quali.sup$coord, color = "red")
a
#now do it but put the group ####
b=fviz_pca_biplot(gbs_medip_pca.1, 
                # Individuals
                geom.ind = "point",
                fill.ind = normalized_count_all.1$group, col.ind = "black",pointsize = "contrib",
                pointshape = 21,
                palette = "jco",
                # Variables
                col.var = "contrib",
                gradient.cols = "RdYlBu",select.var = list(contrib = 150),
                legend.title = list(fill = "Group", color = "Contribution of the\nwindows to the PC")
)
b=b+xlim(c(-4,5))
b=b+ylim(c(-5,5))
b
individuals=c("F2.9","F2.8","F2.7","F2.5","F2.4","F2.3","F2.27","F2.25","F2.24","F2.22","F2.21","F2.20","F2.2","F2.18",
  "F2.17","F2.15","F2.14","F2.13","F2.12","F2.11","F2.10","F1.9","F1.8","F1.7","F1.4",
  "F1.3","F1.24","F1.23","F1.21","F1.20","F1.19","F1.18","F1.17","F1.16","F1.15","F1.14","F1.13","F1.12","F1.1")

fviz_pca_ind(gbs_medip_pca.1, pointsize = "contrib", 
             pointshape = 21, fill = "#E7B800",
             repel=TRUE,select.ind=list(name=individuals)
)

fviz_pca_biplot(gbs_medip_pca.1, select.ind = list(contrib = 5), 
                select.var = list(contrib = 5),
                ggtheme = theme_minimal())

## Intergenerational STATS,let's do F1vsF2 in control and obese group #####
stats.all.gen=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/MACS3.final.count.matrix.featureCounts.txt",header = T)
stats.all.gen=fread("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/MACS3.final.count.matrix.featureCounts.txt",header = T)
groups.mus.musculus=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T,fill = T)
groups.mus.musculus=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T,fill = T)
stats.all.gen <- stats.all.gen[!grepl("random", stats.all.gen$end), ]
stats.all.gen <- stats.all.gen %>% select(where(~ all(!is.na(.))))
location.stats.all.gen <- paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$end)
stats.all.gen$location=paste0(stats.all.gen$chr, ":", stats.all.gen$start, "-", stats.all.gen$end)
individuals.stats.all.gen=as.data.frame(colnames(stats.all.gen[,4:58]))
colnames(individuals.stats.all.gen)=c("sample")
groups.mus.musculus=left_join(individuals.stats.all.gen,groups.mus.musculus,by="sample")

#let's do first the control ####
gen <- "Control"
ind_gen <- data.frame(v1=groups.mus.musculus$sample[groups.mus.musculus$sex == gen])
group_gen <- data.frame(v1=groups.mus.musculus$patient[groups.mus.musculus$sex == gen])
pattern <- "F0\\."
cells<- !grepl(pattern = pattern, ind_gen$v1, ignore.case = TRUE)
ind_gen=ind_gen[cells,]
pattern <- "F0"
cells<- !grepl(pattern = pattern, group_gen$v1, ignore.case = TRUE)
group_gen=group_gen[cells,]
control.stats <- stats.all.gen[, c(..ind_gen,"location")]
row.names(control.stats)=control.stats$location
#index to filter####
non_zero_counts <- apply(control.stats[,1:22], 1, function(x) sum(x != 0))
#filter by having at least 2 ind with counts####
control.stats <- control.stats[non_zero_counts >= 2, ]
edgeR.control <- DGEList(counts=control.stats[,1:22], genes=control.stats$location)
# rownames(edgeR.control$counts) <- rownames(edgeR.control$genes) <- control.stats[,12]
edgeR.control <- calcNormFactors(edgeR.control)
efective.control= as.data.frame(edgeR.control$samples$lib.size*edgeR.control$samples$norm.factors)
edgeR.control.MW=data.frame(mapply(`*`,control.stats[,1:22],t(efective.control)))
p.values.control.mann.whitney=apply(edgeR.control.MW, 1, function(x){wilcox.test(x~group_gen, exact=F)$p.value})
p.values.control.mann.whitney=data.frame(p_value=p.values.control.mann.whitney)
p.values.control.mann.whitney$location=control.stats$location
p.values.control.mann.whitney$corrected_p_value.BH = p.adjust(((p.values.control.mann.whitney$p_value)), method = "BH")
p.values.control.mann.whitney$corrected_p_value = p.adjust(((p.values.control.mann.whitney$p_value)), method = "bonferroni")

## let's do now the obese #####
gen <- "Overnutrition"
ind_gen <- data.frame(v1=groups.mus.musculus$sample[groups.mus.musculus$sex == gen])
group_gen <- data.frame(v1=groups.mus.musculus$patient[groups.mus.musculus$sex == gen])
pattern <- "F0\\."
cells<- !grepl(pattern = pattern, ind_gen$v1, ignore.case = TRUE)
ind_gen=ind_gen[cells,]
pattern <- "F0"
cells<- !grepl(pattern = pattern, group_gen$v1, ignore.case = TRUE)
group_gen=group_gen[cells,]
obese.stats <- stats.all.gen[, c(..ind_gen,"location")]
row.names(obese.stats)=obese.stats$location
#index to filter####
non_zero_counts <- apply(obese.stats[,1:23], 1, function(x) sum(x != 0))
#filter by having at least 2 ind with counts####
obese.stats <- obese.stats[non_zero_counts >= 2, ]
edgeR.obese <- DGEList(counts=obese.stats[,1:23], genes=obese.stats$location)
# rownames(edgeR.obese$counts) <- rownames(edgeR.obese$genes) <- obese.stats[,12]
edgeR.obese <- calcNormFactors(edgeR.obese)
efective.obese= as.data.frame(edgeR.obese$samples$lib.size*edgeR.obese$samples$norm.factors)
edgeR.obese.MW=data.frame(mapply(`*`,obese.stats[,1:23],t(efective.obese)))
p.values.obese.mann.whitney=apply(edgeR.obese.MW, 1, function(x){wilcox.test(x~group_gen, exact=F)$p.value})
p.values.obese.mann.whitney=data.frame(p_value=p.values.obese.mann.whitney)
p.values.obese.mann.whitney$location=obese.stats$location
p.values.obese.mann.whitney$corrected_p_value.BH = p.adjust(((p.values.obese.mann.whitney$p_value)), method = "BH")
p.values.obese.mann.whitney$corrected_p_value = p.adjust(((p.values.obese.mann.whitney$p_value)), method = "bonferroni")
write.table(p.values.control.mann.whitney,file = "C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/p_values_intergenerational_DMRs_control.txt",sep = "\t",quote = F,row.names = F)
write.table(p.values.obese.mann.whitney,file = "C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/p_values_intergenerational_DMRs_obese.txt",sep = "\t",quote = F,row.names = F)

#save(edgeR.obese.MW,file = "C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/Rresults/normalized_counts_obese_intergenerational_DMRs.rda")

#as only in the obese we have DMRs, let's plot those ####
p.values.obese.mann.whitney=p.values.obese.mann.whitney[p.values.obese.mann.whitney$p_value<0.05,]
edgeR.obese.MW$location=obese.stats$location
edgeR.obese.MW=inner_join(edgeR.obese.MW,p.values.obese.mann.whitney,by="location")
edgeR.obese.MW=edgeR.obese.MW[,-c(24:30)]
pattern <- "F0\\."
cells<- !grepl(pattern = pattern, groups.mus.musculus$sample, ignore.case = TRUE)
groups.mus.musculus=groups.mus.musculus[cells,]
# first order by families #
family_map <- groups.mus.musculus %>%
  select(sample, patient, family) %>%
  arrange(patient, family, sample) %>%
  mutate(col_order = row_number()) %>%
  with(setNames(col_order, sample))
ordered_cols <- names(edgeR.obese.MW)[order(family_map[names(edgeR.obese.MW)])]
edgeR.obese.MW <- edgeR.obese.MW[ordered_cols]

row.names(edgeR.obese.MW)=p.values.obese.mann.whitney$location
obese.DMRs=log(edgeR.obese.MW)
obese.DMRs=as.matrix(obese.DMRs)
obese.DMRs[obese.DMRs == -Inf] <- 0
ind=data.frame(sample=row.names(edgeR.obese.MW))
all=inner_join(ind, groups.mus.musculus,by="sample")

col_fun = colorRamp2(c(0, 15.12), c("#fff3f3", "#990000"))
families_ICR=as.numeric(factor(all$family))
col_fun.4 = colorRamp2(1:11, c("black", "#FFC107", "#2196F3", "#FF5722", "#9C27B0", "#E91E63", "#009688", "#FFEB3B", "#6A5ACD", "#A0522D", "#808000"))
column_ha.3 = HeatmapAnnotation(Family= (families_ICR), col = list(Family = col_fun.4),show_legend = FALSE )
column_ha.1 = HeatmapAnnotation(Generation=anno_block(labels = c("F1","F2"),labels_gp = gpar(col = "black", fontsize = 20)),annotation_legend_param = list(title_size = 20, labels_size = 20))
Heatmap_methlevels.1=Heatmap(obese.DMRs, name = "log(Methylation levels)",height = unit(10, "mm")*5,width = unit(10, "mm")*23,
                             top_annotation  = column_ha.3,bottom_annotation  = column_ha.1,border = TRUE,column_gap = unit(1, "mm"), 
                             col = col_fun,column_split = c(rep("F1_O",11),rep("F2_O",12)),
                             column_title ="Differential Methylated Regions across generations in the treated group",show_column_names = FALSE,row_title = "",
                             cluster_rows = FALSE,cluster_columns = FALSE,
                             heatmap_legend_param = list(title_size = 20, labels_size = 20),
                             row_title_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20),
                             row_names_gp = gpar(fontsize = 20),cell_fun = function(j, i, x, y, width, height, fill) {
                               grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "grey", lwd = 0.5))})
par(mar=c(1,1,1,1))
tiff("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/heatmap_DMRs_intergenerational_obese_families_with_color.tiff", width = 3000, height = 1000, units = "px", res = 200)
draw(Heatmap_methlevels.1, heatmap_legend_side = "left", padding = unit(c(4, 3, 3, 29), "mm"))
dev.off()

### density plots of the methylated windows ####
load("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/normalized_counts_and_metadata.rda")
normalized_count_all=normalized_count_all[,-c(10601:10604)]
groups.mus.musculus=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T,fill = T)
#groups.mus.musculus=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T,fill = T)
individuals.stats.all.gen=data.frame(sample=normalized_count_all$sample)
groups.mus.musculus=left_join(individuals.stats.all.gen,groups.mus.musculus,by="sample")
normalized_count_all=left_join(normalized_count_all,groups.mus.musculus,by="sample")
#all together in one plot ####
all.in.one=pivot_longer(normalized_count_all, cols = 2:10600, names_to = "windows", values_to = "methylation_levels")
tiff(filename = "C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/density_plot_all_CG_VS_OG.tiff",width = 2500, height = 2000, units = "px",res =200)
ggplot(all.in.one, aes(x = log(methylation_levels), color = interaction(sex,patient))) +
  geom_density(alpha=0.5, size=1.5) +
  labs(title = "",
       x = "log(Normalized methylation Level)",
       y = "Density", color="Group") +scale_color_manual(values=c("#e74269", "#85bc37","#9b1746","#204c26")) +
  theme_minimal()+theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
dev.off()
## in the F1 generation ####
gen <- "F1"
ind_gen <- groups.mus.musculus$sample[groups.mus.musculus$patient == gen]
F1 <- normalized_count_all %>% filter(sample %in% ind_gen)
F1=pivot_longer(F1, cols = 2:10600, names_to = "windows", values_to = "methylation_levels")
F1=F1%>%arrange(windows)
# do the density plot ####
tiff(filename = "C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/density_plot_f1_CG_VS_OG.tiff",width = 2500, height = 2000, units = "px",res =200)
ggplot(F1, aes(x = log(methylation_levels), color = factor(sex))) +
  geom_density(alpha=0.5, size=1.5) +
  labs(title = "Density Plot of Methylation Levels in F1",
       x = "Methylation Level",
       y = "Density", color="Group") +scale_color_manual(values = c("Control" = "#e74269", "Overnutrition" = "#85bc37")) +
  theme_minimal()+theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
dev.off()

## in the F2 generation ####
gen <- "F2"
ind_gen <- groups.mus.musculus$sample[groups.mus.musculus$patient == gen]
F2 <- normalized_count_all %>% filter(sample %in% ind_gen)
F2=pivot_longer(F2, cols = 2:10600, names_to = "windows", values_to = "methylation_levels")
F2=F2%>%arrange(windows)
# do the density plot ####
tiff(filename = "C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/density_plot_F2_CG_VS_OG.tiff",width = 2500, height = 2000, units = "px",res =200)
ggplot(F2, aes(x = log(methylation_levels), color = factor(sex))) +
  geom_density(alpha=0.5, size=1.5) +
  labs(title = "Density Plot of Methylation Levels in F2",
       x = "Methylation Level",
       y = "Density", color="Group") +scale_color_manual(values = c("Control" = "#e74269", "Overnutrition" = "#85bc37")) +
  theme_minimal()+theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
dev.off()

# calculate the AUC of these plots #####
groups <- unique(interaction(all.in.one$sex, all.in.one$patient))
auc_results <- data.frame(Group = character(), AUC = numeric())

for (group in groups) {
  # Subset the data for the current group
  subset_data <- subset(all.in.one, interaction(sex, patient) == group)
  
  # Calculate the density
  density_estimate <- density(log(subset_data$methylation_levels))
  
  # Calculate the AUC using the trapezoidal rule
  auc <- sum(diff(density_estimate$x) * (density_estimate$y[-1] + density_estimate$y[-length(density_estimate$y)])) / 2
  
  # Store the results
  auc_results <- rbind(auc_results, data.frame(Group = group, AUC = auc))
}

######################################################################################
# see methylation levels  paired with types of RE to see if there is any 
#influence
#maybe also do a correlation?
library(data.table)
file_paths <- list.files(path="C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/repeated_elements/" , pattern="\\.RE.UCSC.win.MACS3.MQ10.bed",full.names = T)
list_bed_meth=lapply(file_paths,function(x){fread(x,header = F) })
file_paths <- gsub("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/repeated_elements/", "", file_paths)
file_paths <- gsub("C13", "", file_paths)
file_paths <- gsub(".RE.UCSC.win.MACS3.MQ10.bed", "", file_paths)
file_paths <- gsub("_", ".", file_paths)
names(list_bed_meth) <- file_paths
list_bed_meth <- lapply(list_bed_meth, function(x) {
  x[, -c(5, 6, 8, 9, 10, 11, 12)]
})
my_colnames <- c("chr", "start", "end", "methylation", "type_RE", "RE")
list_bed_meth <- lapply(list_bed_meth, function(df) {
  colnames(df) <- my_colnames
  df
})

#we need to normalize the counts individually ####
list_bed_meth_norm=list()
list_bed_meth_norm <- lapply(list_bed_meth, function(df) {
  dge <- DGEList(counts = df$methylation, genes = df[, c("chr", "start", "end")])
  dge <- calcNormFactors(dge)
  
  # Effective scaling factor
  scaling_factor <- dge$samples$lib.size * dge$samples$norm.factors
  
  # Multiply methylation counts by the factor
  df$meth_norm <- df$methylation * scaling_factor
  
  return(df)
})
#now that the counts are normalized, let's merge by type of RE ####
merged_dt <- rbindlist(list_bed_meth_norm, idcol = "sample")
merged_dt[RE == ".", RE := "not_RE"]

# now let's do the heatmap ####
library(ComplexHeatmap)
library(data.table)
library(tidyr)
library(dplyr)
merged_dt[, sample := gsub("F1", "F0", sample)]
merged_dt[, sample := gsub("F2", "F1", sample)]
merged_dt[, sample := gsub("F3", "F2", sample)]
# need to have the info on the individuals to be able to separate them ####
groups.mus.musculus=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T,fill = T)
merged_dt=left_join(merged_dt,groups.mus.musculus,by="sample")
# save(merged_dt,file = "C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/repeated_elements/merged_norm_count_meth_and_RE.rda")
load("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/repeated_elements/merged_norm_count_meth_and_RE.rda")
merged_dt <- merged_dt[!grepl("F0", merged_dt$sample),] #in here is a dataframe
re_levels <- sort(unique(merged_dt$type_RE[merged_dt$type_RE != "."]))
re_levels <- c(re_levels, ".")
merged_dt=setDT(merged_dt) #convert to datatable
merged_dt[, type_RE := factor(type_RE, levels = re_levels)]
setorder(merged_dt, type_RE)
merged_wider=pivot_wider(merged_dt,names_from = "sample",values_from = c("meth_norm"),values_fn = mean)

# tenemos lineas duplicadas, hay solo ~400 RE
merged_wider_nocoord=merged_wider[,c(5:6,11:ncol(merged_wider))]
meta_cols <- setdiff(names(merged_wider_nocoord)[1:2], "RE")
sample_cols <- names(merged_wider_nocoord)[3:ncol(merged_wider_nocoord)]

merged_compacted <- merged_wider %>%
  group_by(RE) %>%
  summarise(# For metadata, take the first non-NA (or just the first if consistent)
    across(all_of(meta_cols), ~ first(na.omit(.x))),
    # For sample columns, combine across rows keeping the non-NA (like coalesce across rows)
    across(all_of(sample_cols), ~ if (all(is.na(.x))) NA_real_ else .x[!is.na(.x)][1]),
    .groups = "drop")
# save(merged_compacted,file = "C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/repeated_elements/merged_compacted.rda")
load("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/repeated_elements/merged_compacted.rda")
#vamos a hacer maximos y minimos para ver que rango tenemos q poner los colores ####
library(circlize)
#install.packages("magick")
library(magick)

max(merged_dt$meth_norm)
#294064190
min(merged_dt$meth_norm)
#7670
#number of types of RE
merged_compacted %>%
  count(type_RE)
#is not ordered by type of RE! ####
desired_levels <- c("DNA", "LINE", "LTR", "Retroposon", "Satellite", "Simple_repeat", "SINE", ".")
merged_compacted$type_RE <- factor(merged_compacted$type_RE, levels = desired_levels)
merged_compacted_ordered <- merged_compacted %>%
  arrange(type_RE)
merged_compacted_filtered <- merged_compacted_ordered %>%
  filter(type_RE != ".")
merged_compacted_filtered$type_RE <- droplevels(merged_compacted_filtered$type_RE)
unique(merged_compacted_filtered$type_RE)

# matrix_re=as.matrix(merged_wider[,c(11:55)])
matrix_re=as.matrix(merged_compacted_filtered[,3:ncol(merged_compacted_filtered)])
matrix_re_log2 <- log2(matrix_re)
# row.names(matrix_re)=merged_wider$RE
row.names(matrix_re)=merged_compacted_filtered$type_RE

col_fun = colorRamp2(c(min(matrix_re_log2,na.rm = T), max(matrix_re_log2,na.rm = T)), c("#fff3f3", "#990000"))
column_ha = HeatmapAnnotation(Generation=anno_block(labels = c("F1","F1","F2","F2"),
                                                    labels_gp = gpar(col = "black", fontsize = 70)),
                              Group=anno_block(labels = c("Control","Obese","Control","Obese"),
                                               labels_gp = gpar(col = "black", fontsize = 70)),
                              annotation_legend_param = list(title_size = 70, labels_size = 70))

row_ha = rowAnnotation(Type_RE = anno_block(
  panel_fun = function(index, levels) {
    # Ensure levels is a character string
    label_text <- as.character(levels)
    
    # Replace long labels with line-broken versions
    label_text <- sub("Simple_repeat", "Simple repeat", label_text)
    # Adjust font size for smaller labels
    label_fontsize <- ifelse(label_text == "Retroposon", 50,
                             ifelse(label_text == "Satellite", 45, 70))
    # Draw box and centered label
    grid.rect(gp = gpar(fill = "white", col = "black",lwd=2))
    grid.text(label_text, x = 0.5, y = 0.5, just = "center",
              gp = gpar(col = "black", fontsize = label_fontsize))
  },
  width = unit(18, "cm")
))
######################################################################################################
#test trial code
# row_ha = rowAnnotation(Type_RE = anno_block(
#   panel_fun = function(index, levels) {
#     # Ensure levels is a character string
#     label_text <- as.character(levels)
#     
#     # Replace long labels with line-broken versions
#     label_text <- sub("Simple_repeat", "Simple repeat", label_text)
#     # Adjust font size for smaller labels
#     label_fontsize <- ifelse(label_text %in% c("Retroposon", "Satellite"), 50, 70)
#     
#     # Draw box and centered label
#     grid.rect(gp = gpar(fill = "white", col = "black"))
#     grid.text(label_text, x = 0.5, y = 0.5, just = "center",
#               gp = gpar(col = "black", fontsize = label_fontsize))
#   },
#   width = unit(18, "cm")
# ))



#######################################################################################################
Heatmap_methlevels=Heatmap(matrix_re_log2, name = "Methylation levels in\nRepeated elements",
                           height = unit(5, "mm")*356,width = unit(10, "mm")*45,
                           bottom_annotation  = column_ha,border = TRUE,column_gap = unit(3, "mm"),
                           col = col_fun,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),
                           column_title ="Methylation level",show_column_names = FALSE,row_title = "Types of Repeated element",
                           left_annotation = row_ha,show_row_names = F,
                           row_split = merged_compacted_filtered$type_RE,
                           cluster_rows = FALSE,cluster_columns = FALSE,
                           na_col = "grey",use_raster = T,show_heatmap_legend = FALSE,
                           row_title_gp = gpar(fontsize = 70), column_title_gp = gpar(fontsize = 70),
                           row_names_gp = gpar(fontsize = 70),layer_fun = function(j, i, x, y, width, height, fill) {
                             grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "grey", lwd = 0.5))})
# print the legend in another tiff image 


###################################
#test trial code
# Heatmap_methlevels=Heatmap(matrix_re, name = "Methylation levels in\nRepeated elements",
#                            height = unit(5, "mm")*356,width = unit(10, "mm")*45,
#                            bottom_annotation  = column_ha,border = TRUE,column_gap = unit(1, "mm"),
#                            col = col_fun,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),
#                            column_title ="Methylation level",show_column_names = FALSE,row_title = "Types of Repeated element",
#                            left_annotation = row_ha,show_row_names = F,
#                            row_split = merged_compacted_ordered$type_RE,
#                            cluster_rows = FALSE,cluster_columns = FALSE,
#                            na_col = "grey",use_raster = T,show_heatmap_legend = FALSE,
#                            row_title_gp = gpar(fontsize = 70), column_title_gp = gpar(fontsize = 70),
#                            row_names_gp = gpar(fontsize = 70),layer_fun = function(j, i, x, y, width, height, fill) {
#                              grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "grey", lwd = 0.5))})
###################################
tiff("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/heatmap_RE_log2meth_norm_levels_no_legend.tiff", 
     width = 5000, height = 11200, units = "px", res = 150)
# pdf("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/heatmap_RE_meth_norm_levels.pdf",
#     width = 30, height = 100, useDingbats = FALSE)
draw(Heatmap_methlevels,heatmap_legend_side = "right",padding = unit(c(4, 3, 10, 60), "mm"))
dev.off()

#############################################################################

#Now let's do this but with hierarchical trees on each of the types of RE
#we need to convert to 0 all NA otherwise it cannot do the clustering
merged_compacted_filtered <- merged_compacted_filtered %>%
  mutate(across(3:ncol(.), ~ replace_na(., 0)))

matrix_re=as.matrix(merged_compacted_filtered[,3:ncol(merged_compacted_filtered)])
matrix_re_log2 <- log2(matrix_re)
# row.names(matrix_re)=merged_wider$RE
row.names(matrix_re)=merged_compacted_filtered$type_RE

col_fun = colorRamp2(c(min(matrix_re_log2,na.rm = T), max(matrix_re_log2,na.rm = T)), c("#fff3f3", "#990000"))
column_ha = HeatmapAnnotation(Generation=anno_block(labels = c("F1","F1","F2","F2"),
                                                    labels_gp = gpar(col = "black", fontsize = 70)),
                              Group=anno_block(labels = c("Control","Obese","Control","Obese"),
                                               labels_gp = gpar(col = "black", fontsize = 70)),
                              annotation_legend_param = list(title_size = 70, labels_size = 70))

row_ha = rowAnnotation(Type_RE = anno_block(
  panel_fun = function(index, levels) {
    # Ensure levels is a character string
    label_text <- as.character(levels)
    
    # Replace long labels with line-broken versions
    label_text <- sub("Simple_repeat", "Simple repeat", label_text)
    # Adjust font size for smaller labels
    label_fontsize <- ifelse(label_text == "Retroposon", 50,
                             ifelse(label_text == "Satellite", 45, 70))
    # Draw box and centered label
    grid.rect(gp = gpar(fill = "white", col = "black",lwd=2))
    grid.text(label_text, x = 0.5, y = 0.5, just = "center",
              gp = gpar(col = "black", fontsize = label_fontsize))
  },
  width = unit(18, "cm")
))
Heatmap_methlevels=Heatmap(matrix_re_log2, name = "Methylation levels in\nRepeated elements",
                           height = unit(5, "mm")*356,width = unit(10, "mm")*45,
                           bottom_annotation  = column_ha,border = TRUE,column_gap = unit(3, "mm"),
                           col = col_fun,column_split = c(rep("F1_C",10),rep("F1_O",11),rep("F2_C",12),rep("F2_O",12)),
                           column_title ="Methylation level",show_column_names = FALSE,row_title = "Types of Repeated element",
                           left_annotation = row_ha,show_row_names = F,
                           row_split = merged_compacted_filtered$type_RE,
                           cluster_rows = T, clustering_distance_rows="euclidean",cluster_columns = FALSE,
                           na_col = "grey",row_dend_side="right",use_raster = T,show_heatmap_legend = T,
                           row_title_gp = gpar(fontsize = 70), column_title_gp = gpar(fontsize = 70),
                           row_names_gp = gpar(fontsize = 70),layer_fun = function(j, i, x, y, width, height, fill) {
                             grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "grey", lwd = 0.5))})
tiff("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/heatmap_RE_log2meth_norm_levels_clustered_no_legend.tiff", 
     width = 5000, height = 11200, units = "px", res = 150)
# pdf("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/heatmap_RE_meth_norm_levels.pdf",
#     width = 30, height = 100, useDingbats = FALSE)
draw(Heatmap_methlevels,heatmap_legend_side = "right",padding = unit(c(4, 3, 10, 60), "mm"))
dev.off()

tiff("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/heatmap_RE_log2meth_norm_levels_clustered_legend.tiff", 
     width = 5000, height = 11200, units = "px", res = 150)
# pdf("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/heatmap_RE_meth_norm_levels.pdf",
#     width = 30, height = 100, useDingbats = FALSE)
draw(Heatmap_methlevels,heatmap_legend_side = "right",padding = unit(c(4, 3, 10, 60), "mm"))
dev.off()

##################################################################################
# now let's do the stats test to see if we have any of the types of RE significant
load("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/repeated_elements/merged_compacted.rda")
load("C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/repeated_elements/merged_compacted.rda")

# first test normality: Shapiro–Wilk, as we have less than 50 samples
mat <- as.matrix(merged_compacted[, 3:47])
sw <- t(apply(mat, 1, function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]          # drop NA/NaN/Inf
  n <- length(x)
  if (n < 3 || n > 5000) {      # shapiro.test requires 3..5000 obs
    return(c(W = NA_real_, p.value = NA_real_, n = n))
  }
  st <- shapiro.test(x)
  c(W = unname(st$statistic), p.value = st$p.value, n = n)
}))

sw_df <- cbind(merged_compacted[c("RE", "type_RE")], as.data.frame(sw))
non_normal_dist_meth_re=subset(sw_df, p.value < 0.05)
nrow(non_normal_dist_meth_re)/nrow(merged_compacted)
normal_dist_meth_re=subset(sw_df, p.value > 0.05)
nrow(normal_dist_meth_re)/nrow(merged_compacted)

#there is only 95 from 400 that are actually normal

#now let's do the non-normal as all the p-value significant means is not normal, Mann Whitney
non_normal_dist_meth_re=left_join(non_normal_dist_meth_re,merged_compacted,by = c("RE","type_RE"))
normal_dist_meth_re=left_join(normal_dist_meth_re,merged_compacted,by = c("RE","type_RE"))
significant_re=rbind(non_normal_dist_meth_re,normal_dist_meth_re)
f1_over_samples <- groups.mus.musculus %>%
  filter(patient == "F1", sex == "Overnutrition") %>%
  pull(sample)
dist_meth_re_F1_OG= intersect(f1_over_samples, names(significant_re))

f1_cont_samples <- groups.mus.musculus %>%
  filter(patient == "F1", sex == "Control") %>%
  pull(sample)
dist_meth_re_F1_CG= intersect(f1_cont_samples, names(significant_re))

f1_over_df <- significant_re %>%
  select(
    RE, type_RE, any_of(dist_meth_re_F1_OG)
  )

f1_cont_df <- significant_re %>%
  select(
    RE, type_RE, any_of(dist_meth_re_F1_CG)
  )

mw.RE <- vector("list", length(significant_re))

for (i in seq_len(nrow(significant_re))){
  print(i)
  og=f1_over_df[i,]
  cg=f1_cont_df[i,]
  if(og$RE==cg$RE){
    x=as.vector(as.numeric(og[,3:13]))
    x=x[is.finite(x)]
    y=as.vector(as.numeric(cg[,3:12]))
    y=y[is.finite(y)]
    if(length(x)==0L || length(y)==0L){
      mw.RE[[i]]=NA
    }else{
      c.m=wilcox.test(x,y)
      mw.RE[[i]]=c.m$p.value
    }
  }
}

mw.RE=data.frame(p_value=unlist(mw.RE),RE=significant_re$RE,type_RE=significant_re$type_RE)
write.table(mw.RE,file = "C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/repeated_elements/p_values_RE_mann_whitney.txt",row.names=F,quote =F)

#now let's do F2
f2_over_samples <- groups.mus.musculus %>%
  filter(patient == "F2", sex == "Overnutrition") %>%
  pull(sample)
dist_meth_re_F2_OG= intersect(f2_over_samples, names(significant_re))

f2_cont_samples <- groups.mus.musculus %>%
  filter(patient == "F2", sex == "Control") %>%
  pull(sample)
dist_meth_re_F2_CG= intersect(f2_cont_samples, names(significant_re))

f2_over_df <- significant_re %>%
  select(
    RE, type_RE, any_of(dist_meth_re_F2_OG)
  )

f2_cont_df <- significant_re %>%
  select(
    RE, type_RE, any_of(dist_meth_re_F2_CG)
  )

mw.RE.F2 <- vector("list", length(significant_re))

for (i in seq_len(nrow(significant_re))){
  print(i)
  og=f2_over_df[i,]
  cg=f2_cont_df[i,]
  if(og$RE==cg$RE){
    x=as.vector(as.numeric(og[,3:13]))
    x=x[is.finite(x)]
    y=as.vector(as.numeric(cg[,3:12]))
    y=y[is.finite(y)]
    if(length(x)==0L || length(y)==0L){
      mw.RE.F2[[i]]=NA
    }else{
      c.m=wilcox.test(x,y)
      mw.RE.F2[[i]]=c.m$p.value
    }
  }
}

mw.RE.F2=data.frame(p_value=unlist(mw.RE.F2),RE=significant_re$RE,type_RE=significant_re$type_RE)
write.table(mw.RE.F2,file = "C:/Users/User/Box/Templeton_mus_musuculus/GBS-MEDIP/repeated_elements/p_valuesF2_RE_mann_whitney.txt",row.names=F,quote =F)















