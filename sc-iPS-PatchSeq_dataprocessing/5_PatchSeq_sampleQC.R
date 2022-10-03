## This R script:
## Performs Quality Control of the counts data (saved in Seurat object) for filtering cells

## Load packages
if(!require(rstudioapi)){
  install.packages("rstudioapi")
  library(rstudioapi)
}
if(!require(Seurat)){
  install.packages("Seurat")
  library(Seurat)
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
library(scater)

## Specify data and results directories
script_path <- rstudioapi::getActiveDocumentContext()$path
script_dir <- dirname(script_path)  
setwd(script_dir)
maindir <- dirname(script_dir)
datadir <- paste0(maindir, "/9_salmon_nospikes_tximport_transcriptlevel_to_genelevel")
resultsdir <- paste0(maindir, "/10_sampleQC")
dir.create(resultsdir)

########################################
###### Load the unfiltered counts ######
########################################

load(file = paste0(datadir, "/Bardy_salmon_counts_unfiltered.RData"))

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
###### Create Seurat object ######
##################################

salmon.counts <- salmon.counts[,which(!(names(salmon.counts) %in% c("ensembl_gene_id", "chromosome_name", "external_gene_name", "unique_gene_name")))]

Bardy_seurat <- CreateSeuratObject(as.data.frame(as.matrix(salmon.counts)), project = "PatchSeq")
Bardy_seurat@assays$RNA[1:5,1:5] ## look at the first 5 cells and first 5 genes
dim(Bardy_seurat) ## 60237 genes and 53 cells

##################################
###### Calculate QC metrics ######
##################################

## Compute the proportion of reads that are of mitochondrial origin for every cell (percent.mt & percent.ribo)
Bardy_seurat <- PercentageFeatureSet(object = Bardy_seurat, pattern = "^MT-", col.name = "percent.mt") ## store mitochondrial proportion in object meta data

## Compute the number of house keeping genes expressed in each cell -- HKGs are abundant and are usually steadily expressed in cells, thus less sensitive to the high dropout
## housekeepers.txt: list of 98 housekeeping genes compiled in Tirosh et al., 2016, to be used in data pre-processing, to remove sources of unwanted variation -- downloaded from https://github.com/Michorlab/tnbc_scrnaseq/blob/master/data/housekeepers.txt
hkgenes <- read.table("housekeepers.txt") ## load the the list of house keeping genes
hkgenes <- as.vector(hkgenes$V1)
hkgenes.found <- which(toupper(rownames(Bardy_seurat@assays[["RNA"]]@counts)) %in% hkgenes) ## identify hkgenes in the data
length(hkgenes.found) ## 95/98 HKGs found
n.expressed.hkgenes <- Matrix::colSums(Bardy_seurat@assays[["RNA"]]@counts[hkgenes.found, ] > 0) # sum the number of detected house keeping genes for each cell, then add this to the meta data
Bardy_seurat <- AddMetaData(object = Bardy_seurat, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")

#####################
###### PLOT QC ######
#####################

## Plot of % mitochondrial reads per cell VS # of genes detected
scatterPlot <- ggplot(Bardy_seurat@meta.data, aes(x = nFeature_RNA, y = percent.mt)) +
  theme_light() +
  geom_point(size = 2) +
  geom_pointdensity(adjust = 0.5) +
  scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'Spectral')))(8)) +
  theme(legend.position = "none", axis.text.x = element_text(size = 22, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 22), axis.title.x = element_text(size = 26), axis.title.y = element_text(size = 26), panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  xlab("# of genes detected") + ylab("% of mitochondrial reads") +
  scale_x_continuous(breaks = seq(0, 15000, by = 1000)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  coord_cartesian(xlim = c(0, 15000), ylim = c(0, 100)) +
  geom_vline(xintercept = 1000, linetype = "dashed", col = "red", size = 1)
scatterPlot
ggsave(scatterPlot, filename = paste(resultsdir, "PatchSeq_genesdetected_versus_percentmito.pdf", sep = "/"), height = 6, width = 6, dpi = 300, useDingbats = FALSE)

## Plot of # of HKG genes expressed VS # of genes detected
scatterPlot <- ggplot(Bardy_seurat@meta.data, aes(x = nFeature_RNA, y = n.exp.hkgenes)) +
  theme_light() +
  geom_point(size = 2) +
  geom_pointdensity(adjust = 0.5) +
  scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'Spectral')))(8)) +
  theme(legend.position = "none", axis.text.x = element_text(size = 22, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 22), axis.title.x = element_text(size = 26), axis.title.y = element_text(size = 26), panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  xlab("# of genes detected") + ylab("# of HK genes expressed\n") +
  scale_x_continuous(breaks = seq(0, 15000, by = 1000)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  coord_cartesian(xlim = c(0, 15000), ylim = c(0, 100))+
  geom_hline(yintercept = 65, linetype = "dashed", col = "red", size = 1) + geom_vline(xintercept = 1000, linetype = "dashed", col = "red", size = 1)
scatterPlot
ggsave(scatterPlot, filename = paste(resultsdir, "PatchSeq_genesdetected_versus_nexpHKgenes.pdf", sep = "/"), height = 6, width = 6, dpi = 300, useDingbats = FALSE)

## Plot of relationship between # of genes detected (y axis) and # of read salmon.counts (x axis) --> to detect doublets
scatterPlot <- ggplot(Bardy_seurat@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  theme_light() +
  geom_point(size = 2) +
  geom_pointdensity(adjust = 0.5) +
  scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'Spectral')))(8)) +
  theme(legend.position = "none", axis.text.x = element_text(size = 22, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 22), axis.title.x = element_text(size = 26), axis.title.y = element_text(size = 26), panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  xlab("# of reads per cell") + ylab("# of genes detected\n") +
  scale_x_continuous(breaks = seq(0, 20000000, by = 2000000)) +
  scale_y_continuous(breaks = seq(0, 15000, by = 1500)) +
  coord_cartesian(xlim = c(0, 20000000), ylim = c(0, 15000))+
  geom_hline(yintercept = 1000, linetype = "dashed", col = "red", size = 1)
scatterPlot
ggsave(scatterPlot, filename = paste(resultsdir, "PatchSeq_reads_versus_genes_per_cell.pdf", sep = "/"), height = 6, width = 6, dpi = 300, useDingbats = FALSE)

## Make histogram of number of cells against # of genes detected
histoPlot <- ggplot(Bardy_seurat@meta.data, aes(nFeature_RNA)) +
  theme_light() +
  geom_histogram(binwidth = 250, fill = "grey", col = "black") +
  xlab("# of genes detected") + ylab("# of cells") +
  theme(legend.position = "none", axis.text.x = element_text(size = 22, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 22), axis.title.x = element_text(size = 22), axis.title.y = element_text(size = 22), panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  scale_x_continuous(breaks = seq(0, 15000, by = 1500)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) +
  coord_cartesian(xlim = c(0, 15000), ylim = c(0, 10)) +
  geom_vline(xintercept = 1000, linetype = "dashed", col = "red", size = 1)
histoPlot
ggsave(histoPlot, filename = paste(resultsdir, "PatchSeq_histogram_genes_per_cell.pdf", sep = "/"), height = 6, width = 6, dpi = 300, useDingbats = FALSE)

##################################
#### SET FILTERING THRESHOLDS #####
##################################

Bardy_seurat$seurat_QC <- c("EMPTY")
Bardy_seurat@meta.data[which(Bardy_seurat@meta.data$n.exp.hkgenes < 65 | Bardy_seurat@meta.data$nFeature_RNA < 1000 | Bardy_seurat@meta.data$nFeature_RNA > 8000),]$seurat_QC <- "FAIL"
Bardy_seurat@meta.data[which(Bardy_seurat@meta.data$seurat_QC != "FAIL"),]$seurat_QC <- "PASS"

## Save the metadata
write.csv(Bardy_seurat@meta.data, file = paste(resultsdir, "PatchSeq_Seurat_QCmetrics.csv", sep = "/"))

############################
#### FILTERING THE DATA ####
############################

dim(Bardy_seurat) ## 60237 genes and 53 cells BEFORE filtering

## Filter the cells based on the quality control metrics: filtering based on (1) # of detected genes and (2) no. of housekeeping genes expressed
Bardy_seurat_filtered <- subset(x = Bardy_seurat, subset = nFeature_RNA >= 1000 & nFeature_RNA <= 8000 & n.exp.hkgenes >= 65)
## We will remove cells with < 1000 genes detected, with less than 65 hKG expressed and with >8,000 genes expressed

dim(Bardy_seurat_filtered) ## 60237 genes and 46 cells AFTER filtering -- 7 cells got removed
table(Bardy_seurat_filtered$orig.ident)

## Data visualization post-filtering
VlnPlot(object = Bardy_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","n.exp.hkgenes"), ncol = 4,  pt.size = 0.1)
dev.copy(pdf, file = paste(resultsdir, "PatchSeq_Seurat_filtered_QCmetrics.pdf", sep = "/"), width = 10, height = 6)
dev.off()

## Save the filtered data
save(Bardy_seurat_filtered, file = paste(resultsdir, "Bardy_Seurat_counts_filtered.RData", sep = "/"))
