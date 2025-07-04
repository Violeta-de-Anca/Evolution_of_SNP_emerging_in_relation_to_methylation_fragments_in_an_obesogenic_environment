setwd("/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/pathway_enrichment")

require(devtools)
  #install_bitbucket("sonnhammergroup/anubix",lib="/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/bin")
#library(vctrs, lib.loc = "/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/bin/")
library(data.table)
#install_version("vctrs", lib = "/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/bin/", version = "0.6.4")
library(ANUBIX, lib.loc = "/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/bin/")
#install.packages("vctrs", lib = "/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/bin/")
#source("https://bioconductor.org/biocLite.R")
library(KEGGREST)
library(reactome.db)
library(neat)
library(msigdb)
library(ExperimentHub)
library(GSEABase)
# Set all the networks and gene sets #####
#reactome_set=reactome.db
#list_organism=keggList("organism")
#list_organism[grep("Mus musculus",list_organism)]
#list_organism[19692]
#kegg_set=keggList("mmu")
##############################################
kegg_pathway=keggLink("pathway","mmu")
kegg_pathway_genes= keggLink("mmu", "pathway")
kegg_all_pathways_mus_musculus=data.frame(genes=kegg_pathway_genes,pathway=names(kegg_pathway_genes))
##############################################
#gsea_genes=
eh = ExperimentHub()
query(eh , 'msigdb')
##############################################
metabolic_disease_set=fread("/proj/naiss2024-23-57/reference_genomes/neural_network/C0025517_disease_gda_summary.tsv")
load("/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/pathway_enrichment/funcoup5.rda")
#network=fread("/proj/naiss2024-23-57/reference_genomes/mus_musculus/funcoup/FC5.0_M.musculus_full")
#save(network,file="/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/pathway_enrichment/funcoup5.rda")

# First get how many link does it have the metabolic disease set in the network ####
links_matrix_disgenet=anubix_links(network,metabolic_disease_set[,c(1,5)],cutoff = 0.75,network_type = "weighted")
links_matrix_kegg=anubix_links(network,kegg_all_pathways_mus_musculus,cutoff = 0.75,network_type = "weighted")
#save(links_matrix, file="/proj/naiss2024-23-57/ICR_male_lineage/GBS-MEDIP/pathway_enrichment/link_matrix.rda")

# Load all individuals ####


# let's try neat as anubix doesn't work ####
neat(network = network, alist = , blist = )

