#circos plot ####
setwd("C:/Users/viode560/Box/Templeton_mus_musuculus/GBS_ICR_Sperm/cnv_detection/circos_plot_files")

library(circlize)
library(dplyr)
# install.packages("ComplexHeatmap")
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
set.seed(123)
load("CNV_annotated_curated.rda")
CNV_annotated$general_type=CNV_annotated$family.x
table(CNV_annotated$general_type)
CNV_annotated$general_type[CNV_annotated$family.x%in%c("rRNA","scRNA","snRNA","srpRNA","tRNA")]="RE_RNA"
CNV_annotated$general_type[CNV_annotated$family.x%in%c("Satellite")]="Simple_repeat"
CNV_annotated=CNV_annotated[!CNV_annotated$general_type=="RC",]
table(CNV_annotated$general_type)
load("cnv_f1_f2_for_graph.rda")
filtered_calls=filtered_calls[,2:6]
names(filtered_calls)[names(filtered_calls)=="ind"]="sample"
names(filtered_calls)[names(filtered_calls)=="chr"]="chrom"
CNV_annotated=merge(CNV_annotated,filtered_calls,by=c("sample","chrom","start","end"),all.x=T)

# circos plot per generation ####
# The structure of the dataframe should be chrom, start, end, whatever else!!! ####
# F1_CNV=CNV_annotated%>%filter(patient == "F1")
# F1_CNV=F1_CNV[,c(2,3,4,16)]
# circos.initializeWithIdeogram(species = "mm39")
# color_CNV_log=colorRamp2(c(-1,0,1),c("blue","white","red"))
# circos.genomicHeatmap(F1_CNV,col = color_CNV_log,side = "inside")

#subset the generations and groups
F1_obese=CNV_annotated%>%filter(patient == "F1", sex == "Overnutrition")
chr_order <- c(paste0("chr", 1:19), "chrX", "chrY")
F1_obese_log_cnv=F1_obese[,c(2,3,4,16)]
F1_obese_log_cnv=distinct(F1_obese_log_cnv)
F1_obese_log_cnv=F1_obese_log_cnv %>%
  mutate(chrom = factor(chrom, levels = chr_order)) %>%
  arrange(chrom, start, end)
F1_obese_type_cnv=F1_obese[,c(2,8,9,15)]
#F1_obese_type_cnv=distinct(F1_obese_type_cnv)
F1_obese_type_cnv=F1_obese_type_cnv %>%
  mutate(chrom = factor(chrom, levels = chr_order)) %>%
  arrange(chrom, start_RE, end_RE)
levels(factor(F1_obese_type_cnv$general_type))

#this is to do everything inside the genomic track
# circos.par(
#   track.margin  = c(0,0),  # defaults are often ~0.05 or so
#   track.height  = 0.05,
#   cell.padding = c(0, 0, 0, 0)
# )
# circos.initializeWithIdeogram(species = "mm39",ideogram.height  = 0.04,
#                               plotType= c("ideogram", "axis"))
# color_CNV_log=colorRamp2(c(min(F1_obese_log_cnv$mean_log2),0,
#                            max(F1_obese_log_cnv$mean_log2)),
#                          c("blue","white","red"))
# circos.genomicHeatmap(F1_obese_log_cnv,col = color_CNV_log,side = "inside",
#                       heatmap_height  = 0.05)
# bed_list = split(F1_obese_type_cnv[, c("chrom","start_RE","end_RE")],
#                  F1_obese_type_cnv$general_type)
# type_levels = names(bed_list)
# color_CNV_type = c("#332288","#88CCEE","#117733","#DDCC77",
#                    "#CC6677","#AA4499","#999933","#882255")
# type_colors = setNames(color_CNV_type[seq_along(type_levels)],
#                        type_levels)
# for(i in seq_along(bed_list)) {
#   circos.genomicDensity(
#     bed_list[[i]],
#     col          = type_colors[i],
#     track.height = 0.08  # tweak height as you like
#   )
# }
# circos.clear()

# this is to do the heatmap outside, and the rest inside
circos.initializeWithIdeogram(species = "mm39",plotType = NULL)
color_CNV_log=colorRamp2(c(min(F1_obese_log_cnv$mean_log2),0,
                           max(F1_obese_log_cnv$mean_log2)),
                         c("blue","white","red"))
circos.genomicHeatmap(F1_obese_log_cnv,col = color_CNV_log,side = "outside",
                      heatmap_height  = 0.05)
circos.genomicIdeogram(species = "mm39")
bed_list = split(F1_obese_type_cnv[, c("chrom","start_RE","end_RE")],
                 F1_obese_type_cnv$general_type)
type_levels = names(bed_list)
color_CNV_type = c("#332288","#88CCEE","#117733","#DDCC77",
                   "#CC6677","#AA4499","#999933","#882255")
type_colors = setNames(color_CNV_type[seq_along(type_levels)],
                       type_levels)
for(i in seq_along(bed_list)) {
  circos.genomicDensity(
    bed_list[[i]],
    col          = type_colors[i],
    track.height = 0.05  # tweak height as you like
  )
}
circos.clear()


tiff("CNV_F1_obese_RE_circus_plot.tiff",width = 2000, height = 1500,units = "px", res = 150)
circos.initializeWithIdeogram(species = "mm39",plotType = NULL)
color_CNV_log=colorRamp2(c(min(F1_obese_log_cnv$mean_log2),0,
                           max(F1_obese_log_cnv$mean_log2)),
                         c("blue","white","red"))
circos.genomicHeatmap(F1_obese_log_cnv,col = color_CNV_log,side = "outside",
                      heatmap_height  = 0.05,connection_height = mm_h(1))
circos.genomicIdeogram(species = "mm39")
bed_list = split(F1_obese_type_cnv[, c("chrom","start_RE","end_RE")],
                 F1_obese_type_cnv$general_type)
type_levels = names(bed_list)
color_CNV_type = c("#332288","#88CCEE","#117733","#DDCC77",
                   "#CC6677","#AA4499","#999933","#882255")
type_colors = setNames(color_CNV_type[seq_along(type_levels)],
                       type_levels)
for(i in seq_along(bed_list)) {
  circos.genomicDensity(
    bed_list[[i]],
    col          = type_colors[i],
    track.height = 0.05  # tweak height as you like
  )
}
circos.clear()
lgd_heatmap <- Legend(
  title    = "CNV Variation",
  col_fun  = color_CNV_log,
  # You can pick suitable breaks/labels, here is an example:
  at       = c(
    min(F1_obese_log_cnv$mean_log2),
    0,
    max(F1_obese_log_cnv$mean_log2)
  ),
  labels   = c("Delection", "0", "Duplication")
)
lgd_types <- Legend(
  title     = "Repeated Element\nTypes",
  labels    = type_levels,
  legend_gp = gpar(fill = color_CNV_type[seq_along(type_levels)])
)
combined_legend <- packLegend(lgd_heatmap, lgd_types)
draw(
  combined_legend,
  x    = unit(1, "npc") - unit(1, "mm"),  # Right side of the device
  just = "right"
)
dev.off()

F1_control=CNV_annotated%>%filter(patient == "F1", sex == "Control")
F1_control_log_cnv=F1_control[,c(2,3,4,16)]
F1_control_log_cnv=distinct(F1_control_log_cnv)
F1_control_log_cnv=F1_control_log_cnv %>%
  mutate(chrom = factor(chrom, levels = chr_order)) %>%
  arrange(chrom, start, end)
F1_control_type_cnv=F1_control[,c(2,8,9,15)]
#F1_control_type_cnv=distinct(F1_control_type_cnv)
F1_control_type_cnv=F1_control_type_cnv %>%
  mutate(chrom = factor(chrom, levels = chr_order)) %>%
  arrange(chrom, start_RE, end_RE)

tiff("CNV_F1_control_RE_circus_plot.tiff",width = 2000, height = 1500,units = "px", res = 150)
circos.initializeWithIdeogram(species = "mm39",plotType = NULL)
color_CNV_log=colorRamp2(c(min(F1_control_log_cnv$mean_log2),0,
                           max(F1_control_log_cnv$mean_log2)),
                         c("blue","white","red"))
circos.genomicHeatmap(F1_control_log_cnv,col = color_CNV_log,side = "outside",
                      heatmap_height  = 0.05,connection_height = mm_h(1))
circos.genomicIdeogram(species = "mm39")
bed_list = split(F1_control_type_cnv[, c("chrom","start_RE","end_RE")],
                 F1_control_type_cnv$general_type)
type_levels = names(bed_list)
color_CNV_type = c("#332288","#88CCEE","#117733","#DDCC77",
                   "#CC6677","#AA4499","#999933","#882255")
type_colors = setNames(color_CNV_type[seq_along(type_levels)],
                       type_levels)
for(i in seq_along(bed_list)) {
  circos.genomicDensity(
    bed_list[[i]],
    col          = type_colors[i],
    track.height = 0.05  # tweak height as you like
  )
}
circos.clear()
lgd_heatmap <- Legend(
  title    = "CNV Variation",
  col_fun  = color_CNV_log,
  # You can pick suitable breaks/labels, here is an example:
  at       = c(
    min(F1_control_log_cnv$mean_log2),
    0,
    max(F1_control_log_cnv$mean_log2)
  ),
  labels   = c("Delection", "0", "Duplication")
)
lgd_types <- Legend(
  title     = "Repeated Element\nTypes",
  labels    = type_levels,
  legend_gp = gpar(fill = color_CNV_type[seq_along(type_levels)])
)
combined_legend <- packLegend(lgd_heatmap, lgd_types)
draw(
  combined_legend,
  x    = unit(1, "npc") - unit(1, "mm"),  # Right side of the device
  just = "right"
)
dev.off()

F2_obese=CNV_annotated%>%filter(patient == "F2", sex == "Overnutrition")
F2_obese_log_cnv=F2_obese[,c(2,3,4,16)]
F2_obese_log_cnv=distinct(F2_obese_log_cnv)
F2_obese_log_cnv=F2_obese_log_cnv %>%
  mutate(chrom = factor(chrom, levels = chr_order)) %>%
  arrange(chrom, start, end)
F2_obese_type_cnv=F2_obese[,c(2,8,9,15)]
#F2_obese_type_cnv=distinct(F2_obese_type_cnv)
F2_obese_type_cnv=F2_obese_type_cnv %>%
  mutate(chrom = factor(chrom, levels = chr_order)) %>%
  arrange(chrom, start_RE, end_RE)

tiff("CNV_F2_obese_RE_circus_plot.tiff",width = 2000, height = 1500,units = "px", res = 150)
circos.initializeWithIdeogram(species = "mm39",plotType = NULL)
color_CNV_log=colorRamp2(c(min(F2_obese_log_cnv$mean_log2),0,
                           max(F2_obese_log_cnv$mean_log2)),
                         c("blue","white","red"))
circos.genomicHeatmap(F2_obese_log_cnv,col = color_CNV_log,side = "outside",
                      heatmap_height  = 0.05,connection_height = mm_h(1))
circos.genomicIdeogram(species = "mm39")
bed_list = split(F2_obese_type_cnv[, c("chrom","start_RE","end_RE")],
                 F2_obese_type_cnv$general_type)
type_levels = names(bed_list)
color_CNV_type = c("#332288","#88CCEE","#117733","#DDCC77",
                   "#CC6677","#AA4499","#999933","#882255")
type_colors = setNames(color_CNV_type[seq_along(type_levels)],
                       type_levels)
for(i in seq_along(bed_list)) {
  circos.genomicDensity(
    bed_list[[i]],
    col          = type_colors[i],
    track.height = 0.05  # tweak height as you like
  )
}
circos.clear()
lgd_heatmap <- Legend(
  title    = "CNV Variation",
  col_fun  = color_CNV_log,
  # You can pick suitable breaks/labels, here is an example:
  at       = c(
    min(F2_obese_log_cnv$mean_log2),
    0,
    max(F2_obese_log_cnv$mean_log2)
  ),
  labels   = c("Delection", "0", "Duplication")
)
lgd_types <- Legend(
  title     = "Repeated Element\nTypes",
  labels    = type_levels,
  legend_gp = gpar(fill = color_CNV_type[seq_along(type_levels)])
)
combined_legend <- packLegend(lgd_heatmap, lgd_types)
draw(
  combined_legend,
  x    = unit(1, "npc") - unit(1, "mm"),  # Right side of the device
  just = "right"
)
dev.off()

F2_control=CNV_annotated%>%filter(patient == "F2", sex == "Control")
F2_control_log_cnv=F2_control[,c(2,3,4,16)]
F2_control_log_cnv=distinct(F2_control_log_cnv)
F2_control_log_cnv=F2_control_log_cnv %>%
  mutate(chrom = factor(chrom, levels = chr_order)) %>%
  arrange(chrom, start, end)
F2_control_type_cnv=F2_control[,c(2,8,9,15)]
#F2_control_type_cnv=distinct(F2_control_type_cnv)
F2_control_type_cnv=F2_control_type_cnv %>%
  mutate(chrom = factor(chrom, levels = chr_order)) %>%
  arrange(chrom, start_RE, end_RE)

tiff("CNV_F2_control_RE_circus_plot.tiff",width = 2000, height = 1500,units = "px", res = 150)
circos.initializeWithIdeogram(species = "mm39",plotType = NULL)
color_CNV_log=colorRamp2(c(min(F2_control_log_cnv$mean_log2),0,
                           max(F2_control_log_cnv$mean_log2)),
                         c("blue","white","red"))
circos.genomicHeatmap(F2_control_log_cnv,col = color_CNV_log,side = "outside",
                      heatmap_height  = 0.05,connection_height = mm_h(1))
circos.genomicIdeogram(species = "mm39")
bed_list = split(F2_control_type_cnv[, c("chrom","start_RE","end_RE")],
                 F2_control_type_cnv$general_type)
type_levels = names(bed_list)
color_CNV_type = c("#332288","#88CCEE","#117733","#DDCC77",
                   "#CC6677","#AA4499","#999933","#882255")
type_colors = setNames(color_CNV_type[seq_along(type_levels)],
                       type_levels)
for(i in seq_along(bed_list)) {
  circos.genomicDensity(
    bed_list[[i]],
    col          = type_colors[i],
    track.height = 0.05  # tweak height as you like
  )
}
circos.clear()
lgd_heatmap <- Legend(
  title    = "CNV Variation",
  col_fun  = color_CNV_log,
  # You can pick suitable breaks/labels, here is an example:
  at       = c(
    min(F2_control_log_cnv$mean_log2),
    0,
    max(F2_control_log_cnv$mean_log2)
  ),
  labels   = c("Delection", "0", "Duplication")
)
lgd_types <- Legend(
  title     = "Repeated Element\nTypes",
  labels    = type_levels,
  legend_gp = gpar(fill = color_CNV_type[seq_along(type_levels)])
)
combined_legend <- packLegend(lgd_heatmap, lgd_types)
draw(
  combined_legend,
  x    = unit(1, "npc") - unit(1, "mm"),  # Right side of the device
  just = "right"
)
dev.off()

