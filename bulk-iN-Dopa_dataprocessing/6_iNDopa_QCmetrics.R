## Load packages
if(!require(rstudioapi)){
  install.packages("rstudioapi")
  library(rstudioapi)
}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(ggpointdensity)){
  install.packages("ggpointdensity")
  library(ggpointdensity)
}
if(!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
if(!require(biomaRt)){
  install.packages("biomaRt")
  library(biomaRt)
}

## Specify data and results directories
script_path <- rstudioapi::getActiveDocumentContext()$path
script_dir <- dirname(script_path)  
setwd(script_dir)
maindir <- dirname(script_dir)
datadir <- paste0(maindir, "/8_salmon_nospikes_tximport_transcriptlevel_to_genelevel")
resultsdir <- paste0(maindir, "/9_sampleQC_nospikes")
dir.create(resultsdir)

########################################
###### Load the unfiltered counts ######
########################################

load(paste0(datadir, "/Shong_iN_salmon_counts_unfiltered.RData"))

####################################################
###### Convert Ensembl gene IDs to gene names ######
####################################################

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "http://aug2020.archive.ensembl.org")
genes <- row.names(salmon.counts)
genes_list <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
  values = genes,
  mart = mart)
save(genes_list, file = paste0(resultsdir, "/biomart_aug2020.RData"))
load(paste0(resultsdir, "/biomart_aug2020.RData"))
salmon.counts$ensembl_gene_id <- row.names(salmon.counts)
salmon.counts <- merge(salmon.counts, genes_list[, c("ensembl_gene_id", "external_gene_name", "chromosome_name")], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
length(unique(salmon.counts$ensembl_gene_id)) ## 60237
length(unique(salmon.counts$external_gene_name)) ## 59276 -- some gene names are duplicated
unique(salmon.counts[duplicated(salmon.counts$external_gene_name),]$external_gene_name) ## these duplicated gene names are mostly pseudogenes and small (nucleolar) RNAs
## For these genes, we will keep the original Ensembl gene IDs:
salmon.counts$unique_gene_name <- salmon.counts$external_gene_name
salmon.counts[duplicated(salmon.counts$external_gene_name),]$unique_gene_name <- salmon.counts[duplicated(salmon.counts$external_gene_name),]$ensembl_gene_id
length(unique(salmon.counts$unique_gene_name)) ## 60237
row.names(salmon.counts) <- salmon.counts$unique_gene_name

##################################
###### Calculate QC metrics ######
##################################

hkgenes <- read.table("housekeepers.txt") ## load the the list of house keeping genes
hkgenes <- as.vector(hkgenes$V1)
hkgenes.found <- which(toupper(rownames(salmon.counts)) %in% hkgenes) ## identify hkgenes in the data
length(hkgenes.found) ## 95/98 HKGs found
n.expressed.hkgenes <- Matrix::colSums(salmon.counts[hkgenes.found, -which(names(salmon.counts) %in% c("ensembl_gene_id", "external_gene_name", "chromosome_name", "unique_gene_name"))] > 0) # sum the number of detected house keeping genes for each cell, then add this to the meta data

genes_per_sample <- Matrix::colSums(salmon.counts[, -which(names(salmon.counts) %in% c("ensembl_gene_id", "external_gene_name", "chromosome_name", "unique_gene_name"))] > 0) # count gene only if it has non-zero reads mapped.

write.csv(cbind(n.expressed.hkgenes, genes_per_sample), file = paste0(resultsdir, "/iNDopa_QCmetrics.csv"))