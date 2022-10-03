## SCRIPT NAME: tximport.R
## This R script summarizes the transcript-level results from salmon to gene-level expression using the tximport() function
## NOTE: This R script needs to be sourced
## NOTE: This R script needs to be located is a sub-directory of the main directory, along with the Homo_sapiens_GRCh38_txp2gene.tsv file
## NOTE: The salmon output needs to be located in a sub-directory of the main directory

## Load packages
if(!require(rstudioapi)){
  install.packages("rstudioapi")
  library(rstudioapi)
}
if(!require(tximport)){
  BiocManager::install("tximport")
  library(tximport)
}

## Specify data and results directories
script_path <- rstudioapi::getActiveDocumentContext()$path
script_dir <- dirname(script_path) 
setwd(script_dir)
maindir <- dirname(script_dir)
datadir <- paste0(maindir, "/6_salmon_output_nospikes/iN")
resultsdir <- paste0(maindir, "/8_salmon_nospikes_tximport_transcriptlevel_to_genelevel")
dir.create(resultsdir)

## Import tx2gene data frame correlating transcript IDs (first column) to gene IDs (second column) for use with tximport()
tx2gene <- read.table(paste0(script_dir, "/Homo_sapiens_GRCh38_txp2gene.tsv"), sep = "\t")
tx2gene <- tx2gene[,c(1,2)]
names(tx2gene) <- c("TxID", "GeneID")
tx2gene$TxID <- gsub("\\..*", "", tx2gene$TxID); tx2gene$GeneID <- gsub("\\..*", "", tx2gene$GeneID) ## remove version information to facilitate matching with the tx id
length(unique(tx2gene$TxID)) ## 229,518 unique transcript IDs
length(unique(tx2gene$GeneID)) ## 60,710 unique gene IDs

# Get gene-level expression from transcript-level results of salmon
sampleID <- list.files(datadir)
files.salmon <- file.path(datadir, sampleID, "quant.sf")
all(file.exists(files.salmon)) ## TRUE
names(files.salmon) <- sampleID
txi.salmon <- tximport(files.salmon, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
head(txi.salmon$counts)
## txi.salmon$counts will be the count table from the original salmon quantification, but gene-level summarized
## txi.salmon$abundance will be TPM table from the original salmon quantification, but gene-level summarized

## Export the txi.salmon object
save(txi.salmon, file = paste0(resultsdir, "/Shong_iN_salmon_txim.RData"))

## Retrieve and export the counts table
salmon.counts <- txi.salmon$counts
salmon.counts <- as.data.frame(salmon.counts)
write.csv(salmon.counts, file = paste0(resultsdir, "/Shong_iN_salmon_counts_unfiltered.csv"), quote = FALSE, row.names = TRUE)
save(salmon.counts, file = paste0(resultsdir, "/Shong_iN_salmon_counts_unfiltered.RData"))

## Retrieve and export the TPM table
salmon.TPM <- txi.salmon$abundance
salmon.TPM <- as.data.frame(salmon.TPM)
write.csv(salmon.TPM, file = paste0(resultsdir, "/Shong_iN_salmon_TPM_unfiltered.csv"), quote = FALSE, row.names = TRUE)
save(salmon.TPM, file = paste0(resultsdir, "/Shong_iN_salmon_TPM_unfiltered.RData"))
