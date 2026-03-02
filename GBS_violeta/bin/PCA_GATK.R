library(tidyr)
library(tidyselect)
library(dplyr)
library(ggplot2)
library("viridis") 
library(gplots)
library(qqman)
library(data.table)
library(diffdf)
library(useful)
library(stringr)
library(rrapply)
library(ggrepel)
library(patchwork)
library(plotrix)
library(ggplot2)
library(gridExtra)
#install.packages("ggbreak")
library(ggbreak) 
#install.packages("plotrix")
library(FactoMineR)
library(factoextra)

## PCA for Variant calling from GATK ####
gatk.pca=fread("C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/PCA_ICR_GBS_gatk.eigenvec")
groups=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)
#row.names(gatk.pca)=gatk.pca$V1
#row.names(gatk.pca) <- gsub(x = row.names(gatk.pca), pattern = "C13", replacement = "")
groups.eigenvec=gatk.pca[,1]
colnames(groups.eigenvec)=c("sample")
groups.cpg=left_join(groups.eigenvec,groups,by="sample")
eigenvectors=gatk.pca[,3:4]
colnames(eigenvectors)=c("PC1","PC2")
row.names(eigenvectors)=groups.cpg$sample
# ggplot(eigenvectors,aes(x=PC1,y=PC2,color=factor(groups$patient),shape=factor(groups$sex)))+geom_point(size=15)+
#   geom_hline(yintercept = 0)+labs(y= "PC2", x = "PC1")+
#   geom_vline(xintercept = 0)+
#   theme(legend.title = element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=50),axis.text.y =element_text(size=50),legend.text = element_text(size=50),axis.title = element_text(size=50))

#Save it as jpeg with 5000 of width####
ggplot(eigenvectors,aes(x=PC1,y=PC2,color=interaction(groups.cpg$patient,groups.cpg$sex)))+geom_point(size=15)+
  geom_hline(yintercept = 0)+labs(y= "PC2 (1.93%)", x = "PC1 (2.15%)")+
  geom_vline(xintercept = 0)+scale_color_manual(values=c("#d2a9b3", "#e74269", "#9b1746","#c2dd9b","#85bc37","#204c26"))+
  theme(legend.title = element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=50),axis.text.y =element_text(size=50),legend.text = element_text(size=50),axis.title = element_text(size=50))

#####There are two outliers, and the whole F0 is control, so let's change that ####
#outlier line 48 and 21####
eigenvectors=eigenvectors[-c(21,48),]
groups=groups.cpg[-c(21,48),]
row.names(eigenvectors)=groups$sample
ggplot(eigenvectors,aes(x=PC1,y=PC2,color=interaction(groups.cpg$patient,groups.cpg$sex)))+geom_point(size=15)+
  geom_hline(yintercept = 0)+labs(y= "PC2 (1.93%)", x = "PC1 (2.15%)")+
  geom_vline(xintercept = 0)+scale_color_manual(values=c("#d2a9b3", "#e74269", "#9b1746","#c2dd9b","#85bc37","#204c26"))+
  theme(legend.title = element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=50),axis.text.y =element_text(size=50),legend.text = element_text(size=50),axis.title = element_text(size=50))

#Let's plot by coverage ####
coverage=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/coverage_indv.txt")
colnames(coverage)=c("sample","coverage")
groups=left_join(groups,coverage,by="sample")

ggplot(eigenvectors,aes(x=PC1,y=PC2,color=interaction(groups$patient,groups$sex),size=(groups$coverage)))+geom_point()+
  geom_hline(yintercept = 0)+labs(y= "PC2 (1.93%)", x = "PC1 (2.15%)")+
  geom_label(label=groups$sample,show.legend=F,position=position_jitter())+
  geom_vline(xintercept = 0)+scale_color_manual(values=c("#d2a9b3", "#e74269", "#9b1746","#c2dd9b","#85bc37","#204c26"))+
  theme(legend.title = element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=20),axis.text.y =element_text(size=20),legend.text = element_text(size=20),axis.title = element_text(size=20))

ggplot(eigenvectors,aes(x=PC1,y=PC2,color=interaction(groups$patient,groups$sex),size=(groups$coverage)))+geom_point()+
  geom_hline(yintercept = 0)+labs(y= "PC2 (1.93%)", x = "PC1 (2.15%)")+
  geom_text(label=groups$sample,check_overlap=T,show.legend=F,hjust="left")+
  geom_vline(xintercept = 0)+scale_color_manual(values=c("#d2a9b3", "#e74269", "#9b1746","#c2dd9b","#85bc37","#204c26"))+
  theme(legend.title = element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=20),axis.text.y =element_text(size=20),legend.text = element_text(size=20),axis.title = element_text(size=20))

ggplot(eigenvectors,aes(x=PC1,y=PC2,color=interaction(groups$patient,groups$sex),size=(groups$coverage)))+geom_point()+
  geom_hline(yintercept = 0)+labs(y= "PC2 (1.93%)", x = "PC1 (2.15%)")+
  geom_label_repel(data=eigenvectors,label=groups$sample,show.legend=F)+
  geom_vline(xintercept = 0)+scale_color_manual(values=c("#d2a9b3", "#e74269", "#9b1746","#c2dd9b","#85bc37","#204c26"))+
  theme(legend.title = element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=20),axis.text.y =element_text(size=20),legend.text = element_text(size=20),axis.title = element_text(size=20))

ggplot(eigenvectors,aes(x=PC1,y=PC2,color=interaction(groups$patient,groups$sex),label=groups$sample,size=(groups$coverage)))+geom_point()+
  geom_hline(yintercept = 0)+labs(y= "PC2 (1.93%)", x = "PC1 (2.15%)")+
  geom_label_repel(show.legend = F,max.overlaps=53)+
  geom_vline(xintercept = 0)+scale_color_manual(values=c("#d2a9b3", "#e74269", "#9b1746","#85bc37","#204c26"))+
  theme(legend.title = element_blank(),legend.key=element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=20),axis.text.y =element_text(size=20),legend.text = element_text(size=20),axis.title = element_text(size=20))+
  guides(color = guide_legend(override.aes = list(size = 10,fill=NA))) 

#Let's do this but without the F0 ####
gatk.pca=fread("C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/PCA_ICR_GBS_gatk.eigenvec")
gatk.pca=gatk.pca[!grep("F0",gatk.pca$V1),]
groups=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T)
groups=groups[!grep("F0",groups$patient),]
coverage=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/coverage_indv.txt")
colnames(coverage)=c("sample","coverage")
groups=left_join(groups,coverage,by="sample")
groups.eigenvec=gatk.pca[,1]
colnames(groups.eigenvec)=c("sample")
groups.cpg=left_join(groups.eigenvec,groups,by="sample")
eigenvectors=gatk.pca[,3:4]
colnames(eigenvectors)=c("PC1","PC2")
row.names(eigenvectors)=groups.cpg$sample
eigenvectors=eigenvectors[-c(10,37),]
groups=groups.cpg[-c(10,37),]
row.names(eigenvectors)=groups$sample
PCA_noF0=ggplot(eigenvectors,aes(x=PC1,y=PC2,color=interaction(groups$patient,groups$sex)))+geom_point(size=15)+
  geom_hline(yintercept = 0)+labs(y= "PC2 (1.93%)", x = "PC1 (2.15%)")+
  geom_vline(xintercept = 0)+scale_color_manual(values=c("#e74269", "#9b1746","#85bc37","#204c26"))+
  theme(legend.title = element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=50),axis.text.y =element_text(size=50),legend.text = element_text(size=50),axis.title = element_text(size=50))
tiff(filename = "C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/PCA_noF0_no_depth.tiff",width = 6000, height = 4000, units = "px",res =200)
PCA_noF0
dev.off()

PCA_noF0.1=ggplot(eigenvectors,aes(x=PC1,y=PC2,color=interaction(groups$patient,groups$sex),label=groups$sample,size=(groups$coverage)))+geom_point()+
  geom_hline(yintercept = 0)+labs(y= "PC2 (1.93%)", x = "PC1 (2.15%)")+
  geom_label_repel(show.legend = F,max.overlaps=53)+
  geom_vline(xintercept = 0)+scale_color_manual(values=c("#e74269", "#9b1746","#85bc37","#204c26"))+
  theme(legend.title = element_blank(),legend.key=element_blank(),panel.background=element_blank(),plot.background=element_blank(),axis.text.x=element_text(size=20),axis.text.y =element_text(size=20),legend.text = element_text(size=20),axis.title = element_text(size=20))+
  guides(color = guide_legend(override.aes = list(size = 10,fill=NA))) 
tiff(filename = "C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/PCA_noF0_depth.tiff",width = 3000, height = 2000, units = "px",res =200)
PCA_noF0.1
dev.off()


### PCA similar to the GBS-MEDIP one to see how the treatment and the family posicionates with the SNP data #####
#install.packages("vcfR")
library(vcfR)

# Read the VCF file
vcf_data <- read.vcfR("C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/ICR_genotypes_filtered_gatk.vcf")
#save(vcf_data,file = "C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/ICR_genotypes_filtered_gatk.rda")
load("C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/ICR_genotypes_filtered_gatk.rda")

genotype_matrix <- extract.gt(vcf_data, element = "GT", as.numeric = TRUE)
colnames(genotype_matrix) <- gsub("C13", "", colnames(genotype_matrix))
colnames(genotype_matrix) <- gsub("F1", "F0", colnames(genotype_matrix))
colnames(genotype_matrix) <- gsub("F2", "F1", colnames(genotype_matrix))
colnames(genotype_matrix) <- gsub("F3", "F2", colnames(genotype_matrix))
genotype_matrix=as.data.frame(t(genotype_matrix))
groups.mus.musculus=fread("C:/Users/User/OneDrive - Uppsala universitet/PhD_projects/Templeton/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T,fill = T)
individuals_GBS=data.frame(sample=rownames(genotype_matrix))
groups.mus.musculus=left_join(individuals_GBS,groups.mus.musculus,by="sample")
#save(genotype_matrix,file = "C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/genotype_matrix_ICR_GBS.rda")
load("C:/Users/User/OneDrive - Uppsala universitet/PhD_projects/Templeton/GBS_ICR_Sperm/gatk_best_practices/genotype_matrix_ICR_GBS.rda")
load("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/genotype_matrix_ICR_GBS.rda")
pattern <- "F0\\."
genotype_matrix <- genotype_matrix[!grepl(pattern, rownames(genotype_matrix)), ]

# try by putting everything that is NA as a random number ####
#genotype_matrix[is.na(genotype_matrix)] = 3

#load("C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/genotype_matrix_ICR_GBS.rda")
library(FactoMineR)
library(factoextra)

#let's do the PCA of the family with the 3 as NA ####
# genotype_matrix$family=groups.mus.musculus$family
# gbs_pca_fam=PCA(genotype_matrix,ncp = 5, graph = F,quali.sup = 393002)
# fviz_eig(gbs_pca_fam,addlabels = TRUE)
# 
# p=fviz_pca_var(gbs_pca_fam, col.var = "contrib",
#                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) #,select.var = list(contrib = 150)
# p <- fviz_add(p, gbs_pca_fam$quali.sup$coord, color = "black")
# p = p + ylim(c(-5,2.5))
# p = p + xlim(c(-3,2.5))
# p = p + guides(color = "none")+ggtitle("Display of variable family with GBS data")
# p

#####################################
# delete every row that has a NA ####
#####################################
genotype_matrix.1 <- genotype_matrix[, colSums(is.na(genotype_matrix)) == 0]
individuals_GBS=data.frame(sample=rownames(genotype_matrix.1))
groups.mus.musculus=left_join(individuals_GBS,groups.mus.musculus,by="sample")

#calculate SNPs for each group ####
snp_control=genotype_matrix.1[rownames(genotype_matrix.1) %in% groups.mus.musculus$sample[groups.mus.musculus$sex == "Control"], ]
sum(snp_control == 1)
snp_obese=genotype_matrix.1[rownames(genotype_matrix.1) %in% groups.mus.musculus$sample[groups.mus.musculus$sex == "Overnutrition"], ]
sum(snp_obese == 1)

#let's do the PCA of the family without any NAs ####
genotype_matrix.1$family=groups.mus.musculus$family
gbs_pca_fam.1=PCA(genotype_matrix.1,ncp = 5, graph = F,quali.sup = 33907)
fviz_eig(gbs_pca_fam.1,addlabels = TRUE)

p=fviz_pca_var(gbs_pca_fam.1,axes = c(2, 3), col.var = "contrib",geom = "arrow",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),select.var = list(contrib = 500))
p <- fviz_add(p, gbs_pca_fam.1$quali.sup$coord, color = "black")
p = p + guides(color = "none")+ggtitle("Display of variable family with GBS data")
p

#let's try to do the same but with the axis break ####
genotype_matrix.1$family=groups.mus.musculus$family
gbs_pca_fam.1=PCA(genotype_matrix.1,ncp = 5, graph = F,quali.sup = 33907)
fviz_eig(gbs_pca_fam.1,addlabels = TRUE)

p=fviz_pca_var(gbs_pca_fam.1,axes = c(2, 3), col.var = "contrib",geom = "arrow",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),select.var = list(contrib = 500))
p <- fviz_add(p, gbs_pca_fam.1$quali.sup$coord, color = "black")
p = p + guides(color = "none")+ggtitle("Display of variable family with GBS data")
p=p+labs(x="PC2 (5.1%)",y="PC3 (4.6%)")+ xlim(c(-56.5,1))
p=p+ylim(c(-31,31))
p=p+ scale_x_break(c(-55.5,-1))
# p=p+ scale_y_break(c(1,17))
p
#to make the association between the coordinates and the family ####
dim.family.gbs=dimdesc(gbs_pca_fam.1, axes=c(2,3),proba = 1)
#test for normality ####
coord.ind=as.data.frame(gbs_pca_fam.1[["ind"]][["coord"]])
coord.ind$sample=row.names(coord.ind)
coord.ind=left_join(coord.ind,groups.mus.musculus,by="sample")
ks.test(coord.ind$Dim.2,"pnorm")
ks.test(coord.ind$Dim.3,"pnorm")
kruskal.test(Dim.2 ~ family, data = coord.ind)
kruskal.test(Dim.3 ~ family, data = coord.ind)


factor(coord.ind$family)
# only groups without any NAs ####
genotype_matrix_group.1=genotype_matrix[, colSums(is.na(genotype_matrix)) == 0]
genotype_matrix_group.1$group=groups.mus.musculus$sex
gbs_pca_group.1=PCA(genotype_matrix_group.1,ncp = 5, graph = F,quali.sup = 33907)
d=fviz_pca_var(gbs_pca_group.1,axes = c(2, 3), col.var = "contrib",geom = "arrow",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),select.var = list(contrib = 400))
d <- fviz_add(d, gbs_pca_group.1$quali.sup$coord, color = "black")
d=d+labs(color = "Contribution\nof SNP\nto the PCs")+ ggtitle("Display of variable treatment")
d=d+labs(x="PC2 (5.1%)",y="PC3 (4.6%)")+ xlim(c(-113,1))
d=d+ scale_x_break(c(-97,-1))
d=d+scale_x_break(c(-111,-100))
d
dim.group.gbs=dimdesc(gbs_pca_group.1, axes=c(2,3),proba = 1)
dim.group.gbs=dimdesc(gbs_pca_group.1, axes=c(4,5),proba = 1)
#test for normality ####
coord.ind.grou=as.data.frame(gbs_pca_group.1[["ind"]][["coord"]])
coord.ind.grou$sample=row.names(coord.ind.grou)
coord.ind.grou=left_join(coord.ind.grou,groups.mus.musculus,by="sample")
ks.test(coord.ind.grou$Dim.2,"pnorm")
ks.test(coord.ind.grou$Dim.3,"pnorm")
kruskal.test(Dim.2 ~ sex, data = coord.ind.grou)
kruskal.test(Dim.3 ~ sex, data = coord.ind.grou)

#combine the two plots ####
group.PCA=as.data.frame(rbind(gbs_pca_group.1[["quali.sup"]][["coord"]],gbs_pca_group.1[["var"]][["coord"]]))
combined_plot = p + d 
tiff(filename = "C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/PCA_SNPs_contribution_family_treatment.tiff",
     height = 1200,width = 2500,res = 150)
print(combined_plot)
dev.off()
# individuals without NAs ####
genotype_matrix_group.1$generation=groups.mus.musculus$patient
gbs_pca_group.gen=PCA(genotype_matrix_group.1,ncp = 5, graph = F,quali.sup = c(33907:33908))
# individuals=c("F1.12","F1.13","F1.14","F1.15",
#               "F1.16","F1.17","F1.18","F1.19","F1.1","F1.20","F1.21","F1.23","F1.24","F1.2","F1.3","F1.4","F1.5","F1.6","F1.7",
#               "F1.8","F1.9","F2.10","F2.11","F2.12","F2.13","F2.14","F2.15","F2.16","F2.17","F2.18","F2.19","F2.1","F2.20","F2.21",
#               "F2.22","F2.24","F2.25","F2.27","F2.2","F2.3","F2.4","F2.5","F2.7","F2.8","F2.9")

a=fviz_pca_ind(gbs_pca_group.gen,axes = c(2, 3),
               repel = TRUE,select.ind=list(name=individuals),
               fill.ind = genotype_matrix_group.1$group,geom.ind = c("point", "text"),palette = "jco",pointshape = 21 
)
a=fviz_add(a, gbs_pca_group.gen$quali.sup$coord, color = "red")
a

### do it like in the GBS-MEDIP PCA with factoextra
gbs_eigenvectors=as.data.frame(gbs_pca_group.gen$ind$coord)
gbs_eigenvectors$group=genotype_matrix_group.1$group
gbs_eigenvectors$gen=genotype_matrix_group.1$generation
genotype_matrix_group.1$family=groups.mus.musculus$family
gbs_eigenvectors$family=genotype_matrix_group.1$family
PCA_noF0_gbs=ggplot(gbs_eigenvectors,aes(x=Dim.2,y=Dim.3,color=interaction(group, gen),
                                                     label=gbs_eigenvectors$family))+
  geom_point(size=7)+xlim(-40,45)+ylim(-40,30)

PCA_noF0_gbs=PCA_noF0_gbs+geom_hline(yintercept = 0)+
  labs(x = "PC2 (5.1%)", y = "PC3 (4.6%)")+
  geom_vline(xintercept = 0)+
  scale_color_manual(values=c("#e74269", "#85bc37","#9b1746","#204c26"))

PCA_noF0_gbs=PCA_noF0_gbs+theme(legend.title = element_blank(),legend.key=element_blank(),
                                            panel.background=element_blank(),
                                            plot.background=element_blank(),
                                            axis.text.x=element_text(size=20),
                                            axis.text.y =element_text(size=20),
                                            legend.text = element_text(size=20),
                                            axis.title = element_text(size=20))+
  guides(color = guide_legend(override.aes = list(size = 10,fill=NA)))

PCA_noF0_gbs=PCA_noF0_gbs+geom_label_repel(show.legend = F,max.overlaps=45)
tiff(filename = "C:/Users/User/OneDrive - Uppsala universitet/PhD_projects/Templeton/GBS_ICR_Sperm/gatk_best_practices/plots/GBS_PCA_with_family_names.tiff",width = 3000, height = 2000, units = "px",res =200)
print(PCA_noF0_gbs)
dev.off()
#Do this with the PC1 and PC2 to say that PC1 did not make sense ####
PCA_noF0_gbs=ggplot(gbs_eigenvectors,aes(x=Dim.1,y=Dim.2,color=interaction(group, gen)))+
  geom_point(size=7)
PCA_noF0_gbs=PCA_noF0_gbs+geom_hline(yintercept = 0)+
  labs(x = "PC1 (26.8%)", y = "PC2 (5.1%)")+
  geom_vline(xintercept = 0)+
  scale_color_manual(values=c("#e74269", "#85bc37","#9b1746","#204c26"))
PCA_noF0_gbs=PCA_noF0_gbs+theme(legend.title = element_blank(),legend.key=element_blank(),
                                panel.background=element_blank(),
                                plot.background=element_blank(),
                                axis.text.x=element_text(size=20),
                                axis.text.y =element_text(size=20),
                                legend.text = element_text(size=20),
                                axis.title = element_text(size=20))+
  guides(color = guide_legend(override.aes = list(size = 10,fill=NA)))
PCA_noF0_gbs=PCA_noF0_gbs+geom_label_repel(aes(label = family),show.legend = F,max.overlaps=45)

#put the plot screen to the max or it will give error!!!!!
tiff(filename = "C:/Users/viode560/Box/Templeton_mus_musuculus/GBS-MEDIP/plots/individuals_PCA_all_meth_win.tiff",width = 3000, height = 2000, units = "px",res =200)
tiff(filename = "C:/Users/User/OneDrive - Uppsala universitet/PhD_projects/Templeton/GBS_ICR_Sperm/gatk_best_practices/plots/PC1_PC2_individuals_PCA_factoextra.tiff",width = 3000, height = 2000, units = "px",res =200)
print(PCA_noF0_gbs)
dev.off()

#################################################################################
#let's do a hieralquical tree ####
library(data.table)
library(dplyr)
#install.packages("ggplot2")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("SNPRelate",force = TRUE)
#install.packages("data.tree")
#install.packages("dartR")
library(data.tree)
library(SNPRelate)
library(dartR)
load("C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/genotype_matrix_ICR_GBS.rda")
load("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/genotype_matrix_ICR_GBS.rda")
groups.mus.musculus=fread("C:/Users/User/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T,fill = T)
groups.mus.musculus=fread("C:/Users/viode560/Box/Templeton_mus_musuculus/info_templeton_proyect/groups_and_generations_mus_musculus.txt",header = T,fill = T)
gen <- "F0"
ind_gen <- groups.mus.musculus$sample[!groups.mus.musculus$patient == gen]
group_gen <- groups.mus.musculus$sex[!groups.mus.musculus$patient == gen]
genotype_matrix <- genotype_matrix[rownames(genotype_matrix) %in% ind_gen, ]
ind_gen=data.frame(sample=ind_gen)
groups.mus.musculus=left_join(ind_gen,groups.mus.musculus,by="sample")

# #let's do the distance matrix with hclust ####
# scaled_matrix=scale(normalized_count_all[,-c(1,10601:10604)])
# distance_matrix <- (dist(normalized_count_all[,-c(1,10601:10604)], method = "euclidean"))
# distance_matrix=as.matrix(distance_matrix)
# hclust(d = distance_matrix)
# so we need to do the dissimilarity matrix better cause eucledian distances is not good ####
#install.packages("fastcluster")
library(fastcluster)
library(SNPRelate)
genotype_matrix=as.matrix(genotype_matrix)
snp_data_gds <- snpgdsCreateGeno("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/snp_data.gds", genmat = (genotype_matrix), sample.id = rownames(genotype_matrix), snp.id = colnames(genotype_matrix), snpfirstdim = F)
genofile <- snpgdsOpen("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/snp_data.gds")
# i am doing identity by state and not identity by descend to see if assuming these individuals are unrelated will get me the families or the treatment ####
#can be used to determine similarity between any two individuals####
ibs <- snpgdsIBS(genofile, num.thread = 2)
#transforming the IBS matrix into format distance matrix ####
dissimilarity_matrix <- 1 - ibs$ibs
dist_object <- as.dist(dissimilarity_matrix)
# do the hierarchical clustering####
hc_fast <- fastcluster::hclust(dist_object)
# Plot the dendrogram as only groups
tiff(filename = "C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/plots/hierarchical_tree_identity_by_state.tiff",width = 3000, height = 2000, units = "px",res =200)
tiff(filename = "C:/Users/viode560/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/plots/hierarchical_tree_identity_by_state.tiff",width = 3000, height = 2000, units = "px",res =200)
plot(hc_fast, labels = group_gen,main="", xlab="")
dev.off()

# Plot the dendrogram as gen + groups
tiff(filename = "C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/plots/hierarchical_tree_identity_by_state.tiff",width = 3000, height = 2000, units = "px",res =200)
tiff(filename = "C:/Users/viode560/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/plots/hierarchical_tree_identity_by_stategensandgroup.tiff",width = 3000, height = 2000, units = "px",res =200)
plot(hc_fast, labels =interaction(groups.mus.musculus$patient,groups.mus.musculus$sex),main="", xlab="")
dev.off()

# Plot the dendrogram as only families
tiff(filename = "C:/Users/User/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/plots/hierarchical_tree_identity_by_state.tiff",width = 3000, height = 2000, units = "px",res =200)
tiff(filename = "C:/Users/viode560/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/gatk_best_practices/plots/hierarchical_tree_identity_by_state_families.tiff",width = 3000, height = 2000, units = "px",res =200)
plot(hc_fast, labels = groups.mus.musculus$family,main="", xlab="")
dev.off()

#install.packages("dendextend")
library(dendextend)
dend <- as.dendrogram(hc_fast)
patient_colors <- groups.mus.musculus$patient
unique_groups <- unique(patient_colors)
group_colors <- rainbow(length(unique_groups))
names(group_colors) <- unique_groups
colored_dend <- dend %>%
  set("labels_col", value = group_colors[patient_colors]) %>%
  set("branches_k_color", value = group_colors[patient_colors], k = length(unique_groups))
plot(colored_dend)
