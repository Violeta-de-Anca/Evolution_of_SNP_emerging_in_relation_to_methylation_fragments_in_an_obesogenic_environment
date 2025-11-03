setwd("/proj/naiss2024-23-57/ICR_male_lineage/GBS_violeta/cnv_detection/GATK/")
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

#let's do the probability of the permutation test, we already have the numbers for the distribution
random_RE=fread("cnv_repeat_counts.tsv")
names(random_RE)=c("coordinate","type_RE","count")
random_RE=random_RE%>%pivot_wider(names_from = type_RE,values_from = count,values_fill = 0)

#we need to unmerge the randomization by coordinates:
bed_files <- sprintf("cnv_coverage_check/cnv_%d.bed", 1:7)
read_bed=function(f){
  fread(f,sep = "\t",header = F)
}

cnv_list=lapply(bed_files,read_bed)

for(i in seq_along(cnv_list)){
  cnv_list[[i]][,coordinate:=sprintf("%s:%d-%d",V1,as.integer(V2),as.integer(V3))]
}

cnv_list=lapply(cnv_list,function(dt){
  as_tibble(dt)%>%left_join(random_RE,by="coordinate")
}
  )

#now get the average of each of the cnvs on the OG group
F2_cnv=fread("F2_cnvs_GATK_annotated_UCSC_RE.bed")
F3_cnv=fread("F3_cnvs_GATK_annotated_UCSC_RE.bed")

F2_cnv=F2_cnv[,c(1:3,5,27)]
F3_cnv=F3_cnv[,c(1:3,5,27)]

F2_cnv$coordinate=sprintf("%s:%d-%d",F2_cnv$V1,as.integer(F2_cnv$V2),as.integer(F2_cnv$V3))
F3_cnv$coordinate=sprintf("%s:%d-%d",F3_cnv$V1,as.integer(F3_cnv$V2),as.integer(F3_cnv$V3))

F2_uniq_cnv=unique(sort(F2_cnv$coordinate))
F3_uniq_cnv=unique(sort(F3_cnv$coordinate))

cnv_merged_coord=fread("cnvs_all_gen_merged.bed")
cnv_merged_coord$coordinate=sprintf("%s:%d-%d",cnv_merged_coord$V1,as.integer(cnv_merged_coord$V2),
                                    as.integer(cnv_merged_coord$V3))

#i need to merge the coordinates of both F2 and F3 in the order we have in cnv_merged_coord
setkey(cnv_merged_coord,V1,V2,V3)
setkey(F2_cnv,V1,V2,V3)
setkey(F3_cnv,V1,V2,V3)
F2_good_coordinates=foverlaps(F2_cnv,cnv_merged_coord,type = "within",nomatch = 0L)
F3_good_coordinates=foverlaps(F2_cnv,cnv_merged_coord,type = "within",nomatch = 0L)

#now let's count per type per individual with the good coordinates
F2_good_coordinates=F2_good_coordinates[,c(4,7,8)]
F3_good_coordinates=F3_good_coordinates[,c(4,7,8)]
F2_good_coordinates[, V27 := trimws(V27)]
F3_good_coordinates[, V27 := trimws(V27)]

F2_cnv_counts=F2_good_coordinates[, .N,by=.(coordinate,V5,V27)][order(coordinate,V5,V27)]
F3_cnv_counts=F3_good_coordinates[, .N,by=.(coordinate,V5,V27)][order(coordinate,V5,V27)]

F2_cnv_counts=F2_cnv_counts%>%pivot_wider(names_from = V27,values_from = N,values_fill = 0)
F3_cnv_counts=F3_cnv_counts%>%pivot_wider(names_from = V27,values_from = N,values_fill = 0)

#now let's do the average to be able to do the permutation test
re_types=c("LTR",  "SINE",  "LINE",   "DNA", "Retroposon", "Simple_repeat")
F2_mean_re_cnv=F2_cnv_counts%>%group_by(coordinate)%>%summarise(across(all_of(re_types), 
                                                                       ~mean(.x,na.rm = T)),n_V5=n(),.groups = "drop")
F3_mean_re_cnv=F3_cnv_counts%>%group_by(coordinate)%>%summarise(across(all_of(re_types), 
                                                                       ~mean(.x,na.rm = T)),n_V5=n(),.groups = "drop")

#now that we have the mean, we can compare it with each of the dataframes in the cnv_list
F2_mean_re_cnv=F2_mean_re_cnv%>%arrange(match(coordinate,cnv_merged_coord$coordinate))
F3_mean_re_cnv=F3_mean_re_cnv%>%arrange(match(coordinate,cnv_merged_coord$coordinate))

#now do the test, two sided basially cause we want both enrichment and depletion
perm_test=function(obs_row,perm_dist){
  sapply(re_types,function(re){
    obs=obs_row[[re]]
    perm=perm_dist[[re]]
    n=length(perm)
    mu=mean(perm,na.rm=T)
    (sum(abs(perm-mu)>=abs(obs-mu))+1)/(n+1)
  })
}

F2_pvals=map2_dfr(seq_len(nrow(F2_mean_re_cnv)),cnv_list,function(i,perm_dist){
  obs_row=F2_mean_re_cnv[i,]
  p=perm_test(obs_row,perm_dist)
  enframe(p,name = "RE",value = "p")%>%mutate(coordinate=obs_row$coordinate)%>%select(coordinate,RE,p)
  
})










