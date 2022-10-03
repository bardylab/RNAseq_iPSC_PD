## This R script performs Quality Control of the Fernandes et al. 2020 data for filtering cells using Seurat

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

## Specify data and results directories
script_path <- rstudioapi::getActiveDocumentContext()$path
script_dir <- dirname(script_path)  
setwd(script_dir)
maindir <- dirname(script_dir)
datadir <- paste0(maindir, "/4_SeuratObject_nospikes")
resultsdir <- paste0(maindir, "/5_sampleQC_nospikes")
dir.create(resultsdir)

###############################################
###### Load the unfiltered Seurat object ######
###############################################

## Load the Seurat object containing the Fernandes 10X data
load(paste0(datadir, "/Fernandes_alevin_counts_unfiltered.RData"))

####################################################
###### Convert Ensembl gene IDs to gene names ######
####################################################

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "http://aug2020.archive.ensembl.org")
genes <- row.names(Fernandes_counts)
genes_list <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
  values = genes,
  mart = mart)
save(genes_list, file = paste0(resultsdir, "/biomart_aug2020.RData"))
load(paste0(resultsdir, "/biomart_aug2020.RData"))
Fernandes_counts$ensembl_gene_id <- row.names(Fernandes_counts)
Fernandes_counts <- merge(Fernandes_counts, genes_list[, c("ensembl_gene_id", "external_gene_name", "chromosome_name")], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
length(unique(Fernandes_counts$ensembl_gene_id)) ## 60237
length(unique(Fernandes_counts$external_gene_name)) ## 59276 -- some gene names are duplicated
unique(Fernandes_counts[duplicated(Fernandes_counts$external_gene_name),]$external_gene_name) ## these duplicated gene names are mostly pseudogenes and small (nucleolar) RNAs
## For these genes, we will keep the original Ensembl gene IDs:
Fernandes_counts$unique_gene_name <- Fernandes_counts$external_gene_name
Fernandes_counts[duplicated(Fernandes_counts$external_gene_name),]$unique_gene_name <- Fernandes_counts[duplicated(Fernandes_counts$external_gene_name),]$ensembl_gene_id
length(unique(Fernandes_counts$unique_gene_name)) ## 60237
row.names(Fernandes_counts) <- Fernandes_counts$unique_gene_name

##################################
###### Create Seurat object ######
##################################

Fernandes_counts <- Fernandes_counts[,which(!(names(Fernandes_counts) %in% c("ensembl_gene_id", "chromosome_name", "external_gene_name", "unique_gene_name")))]

Fernandes_seurat <- CreateSeuratObject(as.data.frame(as.matrix(Fernandes_counts)), project = "10XSeq")
Fernandes_seurat@assays$RNA[1:5,1:5] ## look at the first 5 cells and first 5 genes
dim(Fernandes_seurat) ## 60237 genes and 8949 cells

## Check number of cells from each sample (stored in the orig.ident slot)
table(Idents(Fernandes_seurat)) ## HET_UNT: 4483 cells; WT_UNT: 4466 cells

##################################
###### Calculate QC metrics ######
##################################

## Compute the proportion of reads that are of mitochondrial origin for every cell (percent.mt & percent.ribo)
Fernandes_seurat <- PercentageFeatureSet(object = Fernandes_seurat, pattern = "^MT-", col.name = "percent.mt") ## store mitochondrial proportion in object meta data

## Compute the number of house keeping genes expressed in each cell -- HKGs are abundant and are usually steadily expressed in cells, thus less sensitive to the high dropout
## housekeepers.txt: list of 98 housekeeping genes compiled in Tirosh et al., 2016, to be used in data pre-processing, to remove sources of unwanted variation -- downloaded from https://github.com/Michorlab/tnbc_scrnaseq/blob/master/data/housekeepers.txt
hkgenes <- read.table("housekeepers.txt") ## load the the list of house keeping genes
hkgenes <- as.vector(hkgenes$V1)
hkgenes.found <- which(toupper(rownames(Fernandes_seurat@assays[["RNA"]]@counts)) %in% hkgenes) ## identify hkgenes in the data
length(hkgenes.found) ## 95/98 HKGs found
n.expressed.hkgenes <- Matrix::colSums(Fernandes_seurat@assays[["RNA"]]@counts[hkgenes.found, ] > 0) # sum the number of detected house keeping genes for each cell, then add this to the meta data
Fernandes_seurat <- AddMetaData(object = Fernandes_seurat, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")

#####################
###### PLOT QC ######
#####################

## Plot of % mitochondrial reads per cell VS # of genes detected
scatterPlot <- ggplot(Fernandes_seurat@meta.data, aes(x = nFeature_RNA, y = percent.mt)) +
  theme_light() +
  geom_point(size = 2) +
  geom_pointdensity(adjust = 0.5) +
  scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'Spectral')))(8)) +
  theme(legend.position = "none", axis.text.x = element_text(size = 22, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 22), axis.title.x = element_text(size = 26), axis.title.y = element_text(size = 26), panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  xlab("# of genes detected") + ylab("% of mitochondrial reads") +
  scale_x_continuous(breaks = seq(0, 10000, by = 1000)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  coord_cartesian(xlim = c(0, 10000), ylim = c(0, 100)) +
  geom_hline(yintercept = 20, linetype = "dashed", col = "red", size = 1) + geom_vline(xintercept = 1000, linetype = "dashed", col = "red", size = 1)
scatterPlot
ggsave(scatterPlot, filename = paste(resultsdir, "10XSeq_genesdetected_versus_percentmito.pdf", sep = "/"), height = 6, width = 6, dpi = 300, useDingbats = FALSE)

## Plot of # of HKG genes expressed VS # of genes detected
scatterPlot <- ggplot(Fernandes_seurat@meta.data, aes(x = nFeature_RNA, y = n.exp.hkgenes)) +
  theme_light() +
  geom_point(size = 2) +
  geom_pointdensity(adjust = 0.5) +
  scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'Spectral')))(8)) +
  theme(legend.position = "none", axis.text.x = element_text(size = 22, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 22), axis.title.x = element_text(size = 26), axis.title.y = element_text(size = 26), panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  xlab("# of genes detected") + ylab("# of HK genes expressed\n") +
  scale_x_continuous(breaks = seq(0, 10000, by = 1000)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  coord_cartesian(xlim = c(0, 10000), ylim = c(0, 100))+
  geom_hline(yintercept = 65, linetype = "dashed", col = "red", size = 1) + geom_vline(xintercept = 1000, linetype = "dashed", col = "red", size = 1)
scatterPlot
ggsave(scatterPlot, filename = paste(resultsdir, "10XSeq_genesdetected_versus_nexpHKgenes.pdf", sep = "/"), height = 6, width = 6, dpi = 300, useDingbats = FALSE)

## Plot of relationship between # of genes detected (y axis) and # of read counts (x axis) --> to detect doublets
scatterPlot <- ggplot(Fernandes_seurat@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  theme_light() +
  geom_point(size = 2) +
  geom_pointdensity(adjust = 0.5) +
  scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'Spectral')))(8)) +
  theme(legend.position = "none", axis.text.x = element_text(size = 22, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 22), axis.title.x = element_text(size = 26), axis.title.y = element_text(size = 26), panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  xlab("# of reads per cell") + ylab("# of genes detected\n") +
  scale_x_continuous(breaks = seq(0, 150000, by = 25000)) +
  scale_y_continuous(breaks = seq(0, 9000, by = 1000)) +
  coord_cartesian(xlim = c(0, 150000), ylim = c(0, 9000))+
  geom_hline(yintercept = 1000, linetype = "dashed", col = "red", size = 1) + geom_vline(xintercept = 37500, linetype = "dashed", col = "red", size = 1) +
  geom_hline(yintercept = 8000, linetype = "dashed", col = "red", size = 1)
scatterPlot
ggsave(scatterPlot, filename = paste(resultsdir, "10XSeq_reads_versus_genes_per_cell.pdf", sep = "/"), height = 6, width = 6, dpi = 300, useDingbats = FALSE)

## Make histogram of number of cells against # of genes detected
histoPlot <- ggplot(Fernandes_seurat@meta.data, aes(nFeature_RNA)) +
  theme_light() +
  geom_histogram(binwidth = 250, fill = "grey", col = "black") +
  xlab("# of genes detected") + ylab("# of cells") +
  theme(legend.position = "none", axis.text.x = element_text(size = 22, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 22), axis.title.x = element_text(size = 22), axis.title.y = element_text(size = 22), panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  scale_x_continuous(breaks = seq(0, 9000, by = 1000)) +
  scale_y_continuous(breaks = seq(0, 1000, by = 500)) +
  coord_cartesian(xlim = c(0, 9000), ylim = c(0, 1500)) +
  geom_vline(xintercept = 1000, linetype = "dashed", col = "red", size = 1)
histoPlot
ggsave(histoPlot, filename = paste(resultsdir, "10XSeq_histogram_genes_per_cell.pdf", sep = "/"), height = 6, width = 6, dpi = 300, useDingbats = FALSE)

#######################################################
#### ADD MAPPING QC INFO FOR EACH CELLULAR BARCODE ####
#######################################################

featureDump_HET <- read.csv(paste0(maindir, "/2_salmon_alevin_output_nospikes/HET_UNT/alevin/featureDump.txt"), sep = "\t")
featureDump_HET$CB <- paste0("HET_UNT_", featureDump_HET$CB)
featureDump_WT <- read.csv(paste0(maindir, "/2_salmon_alevin_output_nospikes/WT_UNT/alevin/featureDump.txt"), sep = "\t")
featureDump_WT$CB <- paste0("WT_UNT_", featureDump_WT$CB)
featureDump <- rbind(featureDump_HET, featureDump_WT); rm(featureDump_HET, featureDump_WT)

Fernandes_seurat@meta.data$CB <- row.names(Fernandes_seurat@meta.data)
Fernandes_seurat@meta.data <- merge(Fernandes_seurat@meta.data, featureDump[,c("CB", "MappedReads", "MappingRate")], by.x = "CB", by.y = "CB")
row.names(Fernandes_seurat@meta.data) <- Fernandes_seurat@meta.data$CB
Fernandes_seurat@meta.data$CB <- NULL

###################################
#### SET FILTERING THRESHOLDS #####
###################################

Fernandes_seurat$seurat_QC <- c("EMPTY")
Fernandes_seurat@meta.data[which(Fernandes_seurat@meta.data$n.exp.hkgenes < 65 | Fernandes_seurat@meta.data$nFeature_RNA < 1000 | Fernandes_seurat@meta.data$nFeature_RNA > 8000 | Fernandes_seurat@meta.data$nCount_RNA > 37500),]$seurat_QC <- "FAIL"
Fernandes_seurat@meta.data[which(Fernandes_seurat@meta.data$seurat_QC != "FAIL"),]$seurat_QC <- "PASS"

Fernandes_seurat@meta.data$mapping_QC <- c("EMPTY")
Fernandes_seurat@meta.data[which(Fernandes_seurat@meta.data$MappingRate > 0.5),]$mapping_QC <- "PASS"
Fernandes_seurat@meta.data[which(Fernandes_seurat@meta.data$mapping_QC != "PASS"),]$mapping_QC <- "FAIL"

Fernandes_seurat@meta.data$final_QC <- c("EMPTY")
Fernandes_seurat@meta.data[which(Fernandes_seurat@meta.data$seurat_QC == "PASS" & Fernandes_seurat@meta.data$mapping_QC == "PASS"),]$final_QC <- "PASS"
Fernandes_seurat@meta.data[which(Fernandes_seurat@meta.data$final_QC != "PASS"),]$final_QC <- "FAIL"

## Save the metadata
write.csv(Fernandes_seurat@meta.data, file = paste(resultsdir, "10XSeq_Seurat_QCmetrics.csv", sep = "/"))

############################
#### FILTERING THE DATA ####
############################

dim(Fernandes_seurat) ## 60237 genes and 8949 cells BEFORE filtering

## Filter the cells based on the quality control metrics: filtering based on (1) # of detected genes, (2) no. of housekeeping genes expressed, and (3) number of reads <= 37,500/cell
Fernandes_seurat_filtered <- subset(x = Fernandes_seurat, subset = nFeature_RNA >= 1000 & nFeature_RNA <= 8000 & n.exp.hkgenes >= 65 & nCount_RNA <= 37500)
## We will remove cells with < 1000 genes detected, with > 8000 genes detected, with less than 65 hKG expressed and with > 37,500 reads detected
dim(Fernandes_seurat_filtered) ## 60237 genes and 6200 cells AFTER filtering
table(Fernandes_seurat_filtered$orig.ident) ## 3178 x HET_UNT, 3022 x WT_UNT

## Additionally remove genes with low mapping
Fernandes_seurat_filtered <- subset(x = Fernandes_seurat_filtered, subset = mapping_QC == "PASS")
dim(Fernandes_seurat_filtered) ## 60237 genes and 5315 cells AFTER filtering
table(Fernandes_seurat_filtered$orig.ident) ## 2762 x HET_UNT, 2553 x WT_UNT

## Data visualization post-filtering
VlnPlot(object = Fernandes_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","n.exp.hkgenes"), ncol = 4,  pt.size = 0.1)
dev.copy(pdf, file = paste(resultsdir, "10XSeq_Seurat_filtered_QCmetrics.pdf", sep = "/"), width = 10, height = 6)
dev.off()

## Save the filtered data
save(Fernandes_seurat_filtered, file = paste(resultsdir, "Fernandes_Seurat_counts_filtered.RData", sep = "/"))
