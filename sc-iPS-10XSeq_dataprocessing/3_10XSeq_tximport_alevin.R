## SCRIPT NAME: tximport.R
## This R script summarizes the transcript-level results from alevin to gene-level expression using the tximport() function
## NOTE: This R script needs to be sourced
## NOTE: The path to the alevin output needs to be specified below in the "files.alevin" variable

## Load packages
if(!require(rstudioapi)){
  install.packages("rstudioapi")
  library(rstudioapi)
}
if(!require(tximport)){
  BiocManager::install("tximport")
  library(tximport)
}
if(!require(Seurat)){
  install.packages("Seurat")
  library(Seurat)
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
datadir <- paste0(maindir, "/2_salmon_alevin_output_nospikes")
resultsdir <- paste0(maindir, "/4_SeuratObject_nospikes")
dir.create(resultsdir)

####################################################################
###### [Sample HET_UNT] Get gene-level expression from alevin ######
####################################################################

files.alevin <- paste0(datadir, "/HET_UNT/alevin/quants_mat.gz")
all(file.exists(files.alevin)) ## TRUE
txi.alevin <- tximport(files.alevin, type = "alevin")
row.names(txi.alevin$counts) <- gsub("\\..*", "", row.names(txi.alevin$counts)) ## remove the ENSG version information from the Ensembl gene identifiers

## Export the txi.alevin object
save(txi.alevin, file = paste0(resultsdir, "/HET_UNT_alevin_txim.RData"))

####################################################
###### [Sample HET_UNT] Prepare Seurat object ######
####################################################

## [Sample HET_UNT] Create and save Seurat object from counts matrix
HET_UNT <- CreateSeuratObject(as.data.frame(as.matrix(txi.alevin$counts)), project = "HET_UNT")
HET_UNT@assays$RNA[1:5,1:5] ## look at the first 5 cells and first 5 genes
dim(HET_UNT) ## 60237 genes and 4483 cells
save(HET_UNT, file = paste0(resultsdir, "/HET_UNT_Seurat_counts_unfiltered.RData"))

## [Sample HET_UNT] Create and save Seurat object from TPM matrix
counts <- as.data.frame(HET_UNT@assays$RNA@counts)
tpm <- function(counts) {
  counts / sum(counts) * 1e6
}
tpms <- apply(counts, 2, function(x) tpm(x))
colSums(tpms) ## verify that all TPMs sum up to 1e+06 for each sample
HET_UNT_TPM <- CreateSeuratObject(tpms, project = "HET_UNT")
save(HET_UNT_TPM, file = paste0(resultsdir, "/HET_UNT_Seurat_TPM_unfiltered.RData"))

####################################################################
###### [Sample WT_UNT] Get gene-level expression from alevin ######
####################################################################

rm(files.alevin, txi.alevin, counts, tpms)

files.alevin <- paste0(datadir, "/WT_UNT/alevin/quants_mat.gz")
all(file.exists(files.alevin)) ## TRUE
txi.alevin <- tximport(files.alevin, type = "alevin")
row.names(txi.alevin$counts) <- gsub("\\..*", "", row.names(txi.alevin$counts)) ## remove the ENSG version information from the Ensembl gene identifiers

## Export the txi.alevin object
save(txi.alevin, file = paste0(resultsdir, "/WT_UNT_alevin_txim.RData"))

####################################################
###### [Sample WT_UNT] Prepare Seurat object ######
####################################################

## [Sample WT_UNT] Create and save Seurat object from counts matrix
WT_UNT <- CreateSeuratObject(as.data.frame(as.matrix(txi.alevin$counts)), project = "WT_UNT")
WT_UNT@assays$RNA[1:5,1:5] ## look at the first 5 cells and first 5 genes
dim(WT_UNT) ## 60237 genes and 4466 cells
save(WT_UNT, file = paste0(resultsdir, "/WT_UNT_Seurat_counts_unfiltered.RData"))

## [Sample WT_UNT] Create and save Seurat object from TPM matrix
counts <- as.data.frame(WT_UNT@assays$RNA@counts)
tpm <- function(counts) {
  counts / sum(counts) * 1e6
}
tpms <- apply(counts, 2, function(x) tpm(x))
colSums(tpms) ## verify that all TPMs sum up to 1e+06 for each sample
WT_UNT_TPM <- CreateSeuratObject(tpms, project = "WT_UNT")
save(WT_UNT_TPM, file = paste0(resultsdir, "/WT_UNT_Seurat_TPM_unfiltered.RData"))

##########################################################
###### Merge the Seurat objects for the two samples ######
##########################################################

Fernandes_seurat <- merge(HET_UNT, y = WT_UNT, add.cell.ids = c("HET_UNT", "WT_UNT"))
dim(Fernandes_seurat) ## 60237 genes & 8949 cells
table(Fernandes_seurat$orig.ident)
##
Fernandes_seurat_TPM <- merge(HET_UNT_TPM, y = WT_UNT_TPM, add.cell.ids = c("HET_UNT", "WT_UNT"))
dim(Fernandes_seurat_TPM)
table(Fernandes_seurat_TPM$orig.ident)

## Save
save(Fernandes_seurat, file = paste0(resultsdir, "/Fernandes_Seurat_counts_unfiltered.RData"))
save(Fernandes_seurat_TPM, file = paste0(resultsdir, "/Fernandes_Seurat_TPM_unfiltered.RData"))
##
Fernandes_counts <- as.data.frame(as.matrix(GetAssayData(object = Fernandes_seurat, slot = "counts")))
save(Fernandes_counts, file = paste0(resultsdir, "/Fernandes_alevin_counts_unfiltered.RData"))
##
Fernandes_TPM <- as.data.frame(as.matrix(GetAssayData(object = Fernandes_seurat_TPM, slot = "counts")))
save(Fernandes_TPM, file = paste0(resultsdir, "/Fernandes_alevin_TPM_unfiltered.RData"))
